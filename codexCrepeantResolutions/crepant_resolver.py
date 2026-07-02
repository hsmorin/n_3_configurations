#!/usr/bin/env python3
"""
Diagnostics for crepant resolutions of double octics.

The script reads Macaulay2 source files, extracts a branch octic D when one is
present, asks Macaulay2 for the singular-locus decomposition, and computes
finite-field point counts for the singular double cover t^2 = D in P^3.

The blowup tree produced here is a crepant-candidate tree: every component of
the branch singular locus is converted into the center expected in the standard
double-octic crepant algorithm (double curves and quadruple points). If the
singularity profile does not meet those necessary local multiplicity tests, the
case is reported as an error instead of silently inventing a resolution.
"""

from __future__ import annotations

import argparse
import ast
import itertools
import json
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_PRIMES = (3, 5, 7, 11, 13)
M2 = os.environ.get("M2", "M2")


@dataclass(frozen=True)
class Polynomial:
    variables: tuple[str, ...]
    terms: tuple[tuple[int, tuple[int, ...]], ...]
    source: str

    def eval_mod(self, values: tuple[int, ...], p: int) -> int:
        total = 0
        for coeff, exponents in self.terms:
            term = coeff % p
            for value, exponent in zip(values, exponents, strict=True):
                if exponent:
                    term = (term * pow(value % p, exponent, p)) % p
            total = (total + term) % p
        return total


def normalized_projective_points(dim: int, p: int) -> Iterable[tuple[int, ...]]:
    for first in range(dim + 1):
        for tail in itertools.product(range(p), repeat=dim - first):
            yield (0,) * first + (1,) + tail


def root_count_square(rhs: int, p: int) -> int:
    rhs %= p
    if p == 2:
        return sum(1 for t in range(p) if (t * t - rhs) % p == 0)
    if rhs == 0:
        return 1
    return 2 if pow(rhs, (p - 1) // 2, p) == 1 else 0


def count_double_cover(poly: Polynomial, p: int) -> int:
    if len(poly.variables) != 4:
        raise ValueError(f"finite-field counts require 4 variables, got {poly.variables}")
    return sum(root_count_square(poly.eval_mod(point, p), p) for point in normalized_projective_points(3, p))


def extract_octic(path: Path) -> str | None:
    text = path.read_text()
    matches = re.findall(r"^--\s*D\s*=\s*(.+)$", text, flags=re.MULTILINE)
    if matches:
        return max(matches, key=len).strip()

    # For files that compute D but do not show it in a comment, the caller can
    # still run the M2 file manually; this extractor intentionally avoids
    # executing arbitrary project scripts as part of discovery.
    return None


def infer_variables(expr: str) -> tuple[str, ...]:
    used = set(re.findall(r"\b[a-zA-Z]\w*\b", expr))
    ambient_choices = (("x", "y", "z", "w"), ("y", "z", "w", "v"))
    for choice in ambient_choices:
        if used.issubset(choice):
            return choice
    return tuple(sorted(used))


def normalize_m2_octic_string(expr: str) -> str:
    """Normalize a pasted Macaulay2 toString polynomial for both Python and M2."""
    expr = expr.strip()
    if not expr:
        raise ValueError("empty octic string")

    if (expr[0], expr[-1]) in {('"', '"'), ("'", "'")}:
        expr = expr[1:-1].strip()

    # Accept common pasted forms such as "D = ..." or "o12 = ...".
    assignment = re.match(r"^[A-Za-z]\w*\s*=\s*(.+)$", expr, flags=re.DOTALL)
    if assignment:
        expr = assignment.group(1).strip()

    # Macaulay2 toString output is valid on one logical line, but terminal
    # copies often include wrapped newlines and indentation.
    expr = re.sub(r"\s+", "", expr)

    unsupported = sorted(set(re.findall(r"[{}[\];]", expr)))
    if unsupported:
        raise ValueError(f"unsupported characters in octic string: {''.join(unsupported)}")
    return expr


def parse_polynomial(
    expr: str,
    variables: tuple[str, ...] | None = None,
    expected_degree: int | None = 8,
) -> Polynomial:
    expr = normalize_m2_octic_string(expr)
    names = variables or infer_variables(expr)
    used = set(re.findall(r"\b[a-zA-Z]\w*\b", expr))
    if not used.issubset(names):
        raise ValueError(f"variables {sorted(used)} are not contained in ambient variables {names}")
    if expected_degree == 8 and len(names) != 4:
        raise ValueError(f"expected a homogeneous P^3 octic in four variables, found {names}")

    py_expr = expr.replace("^", "**")
    tree = ast.parse(py_expr, mode="eval")
    allowed_nodes = (
        ast.Expression,
        ast.BinOp,
        ast.UnaryOp,
        ast.Add,
        ast.Sub,
        ast.Mult,
        ast.Pow,
        ast.USub,
        ast.UAdd,
        ast.Name,
        ast.Load,
        ast.Constant,
    )
    for node in ast.walk(tree):
        if not isinstance(node, allowed_nodes):
            raise ValueError(f"unsupported polynomial syntax near {ast.dump(node)}")
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
            if not isinstance(node.right, ast.Constant):
                raise ValueError("only integer exponents are supported")

    def zero() -> dict[tuple[int, ...], int]:
        return {}

    def const(c: int) -> dict[tuple[int, ...], int]:
        return {tuple(0 for _ in names): c} if c else zero()

    def var(name: str) -> dict[tuple[int, ...], int]:
        exponents = [0] * len(names)
        exponents[names.index(name)] = 1
        return {tuple(exponents): 1}

    def add(a, b, sign=1):
        out = dict(a)
        for monomial, coeff in b.items():
            out[monomial] = out.get(monomial, 0) + sign * coeff
            if out[monomial] == 0:
                del out[monomial]
        return out

    def mul(a, b):
        out: dict[tuple[int, ...], int] = {}
        for ma, ca in a.items():
            for mb, cb in b.items():
                monomial = tuple(x + y for x, y in zip(ma, mb, strict=True))
                out[monomial] = out.get(monomial, 0) + ca * cb
        return {m: c for m, c in out.items() if c}

    def power(a, n: int):
        if n < 0:
            raise ValueError("negative exponents are not supported")
        out = const(1)
        base = a
        while n:
            if n & 1:
                out = mul(out, base)
            base = mul(base, base)
            n >>= 1
        return out

    def visit(node):
        if isinstance(node, ast.Expression):
            return visit(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return const(node.value)
        if isinstance(node, ast.Name):
            return var(node.id)
        if isinstance(node, ast.UnaryOp):
            val = visit(node.operand)
            if isinstance(node.op, ast.USub):
                return {m: -c for m, c in val.items()}
            return val
        if isinstance(node, ast.BinOp):
            left = visit(node.left)
            right = visit(node.right)
            if isinstance(node.op, ast.Add):
                return add(left, right)
            if isinstance(node.op, ast.Sub):
                return add(left, right, sign=-1)
            if isinstance(node.op, ast.Mult):
                return mul(left, right)
            if isinstance(node.op, ast.Pow):
                assert isinstance(node.right, ast.Constant)
                return power(left, int(node.right.value))
        raise ValueError(f"unsupported polynomial syntax near {ast.dump(node)}")

    term_map = visit(tree)
    if expected_degree is not None:
        degrees = {sum(exponents) for exponents in term_map}
        if degrees != {expected_degree}:
            label = "an octic" if expected_degree == 8 else f"a degree-{expected_degree} polynomial"
            raise ValueError(f"expected {label}; term degrees are {sorted(degrees)}")
    terms = tuple(sorted((coeff, monomial) for monomial, coeff in term_map.items()))
    return Polynomial(variables=names, terms=terms, source=expr)


def split_top_level_commas(text: str) -> list[str]:
    parts = []
    start = 0
    depth = 0
    for idx, char in enumerate(text):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "," and depth == 0:
            parts.append(text[start:idx].strip())
            start = idx + 1
    parts.append(text[start:].strip())
    return [part for part in parts if part]


def parse_m2_matrix_generators(generators: str) -> list[str]:
    match = re.fullmatch(r"matrix\s*\{\{(.*)\}\}", generators.strip())
    if not match:
        return []
    return split_top_level_commas(match.group(1))


def count_projective_subscheme_points(generators: str, variables: tuple[str, ...], p: int) -> int | None:
    gen_exprs = parse_m2_matrix_generators(generators)
    if not gen_exprs:
        return None
    try:
        gen_polys = [parse_polynomial(expr, variables=variables, expected_degree=None) for expr in gen_exprs]
    except ValueError:
        return None

    count = 0
    for point in normalized_projective_points(len(variables) - 1, p):
        if all(poly.eval_mod(point, p) == 0 for poly in gen_polys):
            count += 1
    return count


def standard_blowup_charts(center: dict) -> list[dict]:
    equations = parse_m2_matrix_generators(center["generators"])
    charts = []
    for pivot_index, pivot in enumerate(equations):
        substitutions = []
        for idx, equation in enumerate(equations):
            if idx != pivot_index:
                substitutions.append(f"{equation} = ({pivot})*u_{center['id']}_{idx}")
        charts.append(
            {
                "pivot": pivot,
                "substitutions": substitutions,
                "cover_variable": f"u = ({pivot})^{center['branch_multiplicity']//2}*u_{center['id']}",
                "branch_transform": f"D_chart = substitute(D, chart_map)/({pivot})^{center['branch_multiplicity']}",
            }
        )
    return charts


def blowup_count_corrections(poly: Polynomial, profile: dict, primes: Iterable[int]) -> dict:
    corrections = {}
    for p in primes:
        contributions = []
        total = 0
        exact = True
        for comp in profile["components"]:
            center_points = count_projective_subscheme_points(comp["generators"], poly.variables, p)
            if center_points is None:
                exact = False
                contribution = None
            elif comp["projective_dimension"] == 1:
                contribution = p * center_points
            elif comp["projective_dimension"] == 0:
                contribution = (p + p * p) * center_points
            else:
                exact = False
                contribution = None

            if contribution is not None:
                total += contribution
            contributions.append(
                {
                    "component": comp["index"],
                    "projective_dimension": comp["projective_dimension"],
                    "center_points": center_points,
                    "correction": contribution,
                    "formula": "p * #C(F_p)"
                    if comp["projective_dimension"] == 1
                    else "(p + p^2) * #P(F_p)"
                    if comp["projective_dimension"] == 0
                    else "unsupported",
                }
            )
        corrections[str(p)] = {
            "status": "exact_first_level" if exact else "partial",
            "total_correction": total if exact else None,
            "contributions": contributions,
        }
    return corrections


def run_m2(script: str, timeout: int = 120) -> str:
    with tempfile.NamedTemporaryFile("w", suffix=".m2", delete=False) as handle:
        handle.write(script)
        m2_path = handle.name
    try:
        proc = subprocess.run(
            [M2, "--script", m2_path],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    finally:
        Path(m2_path).unlink(missing_ok=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stdout.strip())
    return proc.stdout


def m2_singular_profile(poly: Polynomial, characteristic: int = 0) -> dict:
    vars_m2 = ",".join(poly.variables)
    coeff_field = "QQ" if characteristic == 0 else f"ZZ/{characteristic}"
    script = f"""
R = {coeff_field}[{vars_m2}]
D = ({poly.source})
irr = ideal({vars_m2})
S = saturate(ideal(D) + ideal apply(gens R, v -> diff(v, D)), irr)
C = decompose S
branchMultiplicity = J -> (
    m := 0;
    scan(1..8, k -> if D % (J^k) == 0 then m = k);
    m
)
print("PROFILE|characteristic|{characteristic}")
print("PROFILE|component_count|" | toString(#C))
scan(#C, i -> (
    J := C#i;
    print("COMP|" | toString(i) | "|" | toString(dim(R/J)) | "|" | toString(degree J) | "|" | toString(branchMultiplicity J) | "|" | toString(mingens J))
))
"""
    out = run_m2(script)
    components = []
    component_count = 0
    for line in out.splitlines():
        if line.startswith("PROFILE|component_count|"):
            component_count = int(line.rsplit("|", 1)[1])
        elif line.startswith("COMP|"):
            _, idx, cone_dim, degree, multiplicity, gens = line.split("|", 5)
            components.append(
                {
                    "index": int(idx),
                    "cone_dimension": int(cone_dim),
                    "projective_dimension": int(cone_dim) - 1,
                    "degree": int(degree),
                    "branch_multiplicity": int(multiplicity),
                    "generators": gens.strip(),
                }
            )
    return {
        "characteristic": characteristic,
        "component_count": component_count,
        "components": components,
        "raw": out,
    }


def infer_blowup_tree(profile: dict) -> dict:
    components = profile["components"]
    errors: list[str] = []
    centers = []
    exceptional_divisors = 0

    for comp in components:
        pdim = comp["projective_dimension"]
        if pdim == 1:
            kind = "curve"
            expected_multiplicity = 2
            if comp["branch_multiplicity"] == expected_multiplicity:
                exceptional_divisors += 1
            else:
                errors.append(
                    f"curve component {comp['index']} has branch multiplicity {comp['branch_multiplicity']}; expected 2"
                )
        elif pdim == 0:
            kind = "point"
            expected_multiplicity = 4
            if comp["branch_multiplicity"] == expected_multiplicity:
                exceptional_divisors += 1
            else:
                errors.append(
                    f"point component {comp['index']} has branch multiplicity {comp['branch_multiplicity']}; expected 4"
                )
        else:
            errors.append(
                f"component {comp['index']} has projective dimension {pdim}; expected curves or points"
            )
            kind = "unsupported"
            expected_multiplicity = None
        centers.append(
            {
                "id": f"c{comp['index']}",
                "kind": kind,
                "projective_dimension": pdim,
                "degree": comp["degree"],
                "generators": comp["generators"],
                "branch_multiplicity": comp["branch_multiplicity"],
                "crepant_expected_branch_multiplicity": expected_multiplicity,
                "standard_affine_charts": standard_blowup_charts({"id": f"c{comp['index']}", **comp}),
                "children": [],
            }
        )

    # The actual local multiplicities require strict-transform iteration after
    # each blowup. This script therefore reports the first-level necessary tree
    # and marks the result as a candidate unless all centers can be certified by
    # a future local chart pass.
    possible = not errors
    return {
        "status": "candidate" if possible else "error",
        "errors": errors,
        "root": {
            "ambient": "P^3",
            "branch_degree": 8,
            "centers": centers,
        },
        "exceptional_divisor_count_lower_bound": exceptional_divisors,
        "picard_rank_lower_bound": 1 + exceptional_divisors,
    }


def detected_bad_primes(poly: Polynomial, qq_profile: dict, primes: Iterable[int]) -> list[int]:
    bad = [2]
    base_signature = sorted(
        (c["projective_dimension"], c["degree"], c["branch_multiplicity"]) for c in qq_profile["components"]
    )
    for p in primes:
        if p == 2:
            continue
        try:
            fp_profile = m2_singular_profile(poly, p)
            signature = sorted(
                (c["projective_dimension"], c["degree"], c["branch_multiplicity"]) for c in fp_profile["components"]
            )
            if signature != base_signature:
                bad.append(p)
        except Exception:
            bad.append(p)
    return sorted(set(bad))


def analyze_file(path: Path, primes: tuple[int, ...], bad_prime_bound: int) -> dict:
    octic = extract_octic(path)
    result = {
        "case": case_name(path),
        "path": str(path),
        "has_octic": octic is not None,
    }
    if not octic:
        result["status"] = "error"
        result["error"] = "No commented branch octic of the form '-- D = ...' was found."
        return result

    poly = parse_polynomial(octic)
    result["variables"] = poly.variables
    result["term_count"] = len(poly.terms)
    result["octic"] = octic

    profile = m2_singular_profile(poly)
    tree = infer_blowup_tree(profile)
    counts = {str(p): count_double_cover(poly, p) for p in primes}
    corrections = blowup_count_corrections(poly, profile, primes)
    corrected_counts = {
        str(p): counts[str(p)] + corrections[str(p)]["total_correction"]
        for p in primes
        if corrections[str(p)]["total_correction"] is not None
    }
    bad_test_primes = tuple(p for p in primes_upto(bad_prime_bound) if p >= 2)

    result.update(
        {
            "status": tree["status"],
            "singular_locus": profile,
            "blowup_tree": tree,
            "picard_rank": {
                "lower_bound": tree["picard_rank_lower_bound"],
                "note": "Computed as 1 plus first-level exceptional divisors; full rank needs strict-transform chart certification.",
            },
            "bad_reduction": {
                "detected_primes": detected_bad_primes(poly, profile, bad_test_primes),
                "tested_primes": list(bad_test_primes),
                "note": "2 is always included; odd primes are flagged when the singular-locus component signature changes modulo p.",
            },
            "finite_field_counts": counts,
            "singular_model_finite_field_counts": counts,
            "blowup_count_correction": {
                "status": "first_level_from_listed_blowup_centers",
                "note": "For a curve center the blowup replaces the center by a P^1-bundle, contributing p times the center count. For a point center in a threefold it contributes p+p^2. These are first-level center corrections from the ledger below.",
                "by_prime": corrections,
            },
            "corrected_first_level_finite_field_counts": corrected_counts,
            "singular_kulikov_fibers": {
                "status": "not_computed",
                "note": "The repository stores K3 fibration models for many cases, but not a uniform fibration API. This script records the branch resolution diagnostics needed before a reliable Kulikov-fiber classifier can be attached.",
            },
        }
    )
    return result


def analyze_octic_string(
    case: str,
    octic: str,
    primes: tuple[int, ...],
    bad_prime_bound: int,
    variables: tuple[str, ...] | None = None,
) -> dict:
    poly = parse_polynomial(octic, variables=variables)
    profile = m2_singular_profile(poly)
    tree = infer_blowup_tree(profile)
    counts = {str(p): count_double_cover(poly, p) for p in primes}
    corrections = blowup_count_corrections(poly, profile, primes)
    corrected_counts = {
        str(p): counts[str(p)] + corrections[str(p)]["total_correction"]
        for p in primes
        if corrections[str(p)]["total_correction"] is not None
    }
    bad_test_primes = tuple(p for p in primes_upto(bad_prime_bound) if p >= 2)

    return {
        "case": case,
        "path": "command-line octic string",
        "has_octic": True,
        "variables": poly.variables,
        "term_count": len(poly.terms),
        "octic": poly.source,
        "status": tree["status"],
        "singular_locus": profile,
        "blowup_tree": tree,
        "picard_rank": {
            "lower_bound": tree["picard_rank_lower_bound"],
            "note": "Computed as 1 plus first-level exceptional divisors; full rank needs strict-transform chart certification.",
        },
        "bad_reduction": {
            "detected_primes": detected_bad_primes(poly, profile, bad_test_primes),
            "tested_primes": list(bad_test_primes),
            "note": "2 is always included; odd primes are flagged when the singular-locus component signature changes modulo p.",
        },
        "finite_field_counts": counts,
        "singular_model_finite_field_counts": counts,
        "blowup_count_correction": {
            "status": "first_level_from_listed_blowup_centers",
            "note": "For a curve center the blowup replaces the center by a P^1-bundle, contributing p times the center count. For a point center in a threefold it contributes p+p^2. These are first-level center corrections from the ledger below.",
            "by_prime": corrections,
        },
        "corrected_first_level_finite_field_counts": corrected_counts,
        "singular_kulikov_fibers": {
            "status": "not_computed",
            "note": "Direct octic-string mode analyzes the branch model only; no K3 fibration data was supplied.",
        },
    }


def primes_upto(n: int) -> list[int]:
    out = []
    for x in range(2, n + 1):
        for q in range(2, int(x**0.5) + 1):
            if x % q == 0:
                break
        else:
            out.append(x)
    return out


def case_name(path: Path) -> str:
    for part in path.parts:
        if re.fullmatch(r"Q\d+", part):
            return part
    return path.stem


def discover_case_files(data_root: Path) -> list[Path]:
    candidates = sorted(data_root.glob("Q*/*/*.m2")) + sorted(data_root.glob("Q*/*.m2"))
    best: dict[str, Path] = {}
    for path in candidates:
        name = case_name(path)
        if name in best:
            if extract_octic(best[name]) is None and extract_octic(path) is not None:
                best[name] = path
            continue
        best[name] = path
    return [best[name] for name in sorted(best, key=lambda s: int(s[1:]))]


def write_report(result: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    case = result["case"]
    (out_dir / f"{case}.json").write_text(json.dumps(result, indent=2, sort_keys=True))

    lines = [
        f"# {case}",
        "",
        f"Source: `{result['path']}`",
        f"Status: `{result['status']}`",
        "",
    ]
    if result["status"] == "error" and not result.get("has_octic"):
        lines += [result["error"], ""]
    else:
        lines += [
            f"Variables: `{', '.join(result['variables'])}`",
            f"Octic terms: {result['term_count']}",
            "",
        ]
        if "model_search" in result:
            search = result["model_search"]
            lines += [
                "## Double-octic model",
                f"Affine projective chart: `{search.get('chart', '')}`",
                f"Quadratic P1 variable: `{search.get('quadratic_variable', '')}`",
                f"P3 variables from projective coordinates: `{search.get('p3_variables_from_projective_coordinates', '')}`",
                f"Discriminant u-power removed: `{search.get('u_power_removed_from_discriminant', '')}`",
                f"Selected candidate: `{search.get('selected_candidate_index', '')}` of `{search.get('candidate_count', '')}`",
                "",
            ]
        lines += [
            "## Singular locus",
        ]
        for comp in result["singular_locus"]["components"]:
            lines.append(
                f"- Component {comp['index']}: projective dimension {comp['projective_dimension']}, "
                f"degree {comp['degree']}, branch multiplicity {comp['branch_multiplicity']}, "
                f"generators `{comp['generators']}`"
            )
        lines += [
            "",
            "## Blowup tree",
            f"Tree status: `{result['blowup_tree']['status']}`",
            f"Picard rank lower bound: {result['picard_rank']['lower_bound']}",
        ]
        for center in result["blowup_tree"]["root"]["centers"]:
            lines.append(
                f"- {center['id']}: {center['kind']}, degree {center['degree']}, "
                f"branch multiplicity {center['branch_multiplicity']}, "
                f"expected branch multiplicity {center['crepant_expected_branch_multiplicity']}"
            )
            for chart in center.get("standard_affine_charts", []):
                substitutions = "; ".join(chart["substitutions"]) if chart["substitutions"] else "identity"
                lines.append(f"  - chart `{chart['pivot']}`: {substitutions}; {chart['cover_variable']}")
        if result["blowup_tree"]["errors"]:
            lines += ["", "Errors:"]
            lines += [f"- {err}" for err in result["blowup_tree"]["errors"]]
        lines += [
            "",
            "## Bad reduction",
            f"Detected primes: {', '.join(map(str, result['bad_reduction']['detected_primes']))}",
            result["bad_reduction"]["note"],
            "",
            "## Finite-field counts",
        ]
        for p, count in result["finite_field_counts"].items():
            correction = result.get("blowup_count_correction", {}).get("by_prime", {}).get(p, {})
            correction_total = correction.get("total_correction")
            corrected = result.get("corrected_first_level_finite_field_counts", {}).get(p)
            if correction_total is None:
                lines.append(f"- F_{p}: singular model {count}; correction not fully computed")
            else:
                lines.append(f"- F_{p}: singular model {count}; blowup correction {correction_total}; corrected first-level count {corrected}")
        lines += [
            "",
            "## Singular Kulikov fibers",
            result["singular_kulikov_fibers"]["note"],
            "",
        ]
    (out_dir / f"{case}.md").write_text("\n".join(lines))


def write_index(summary: list[dict], out_dir: Path) -> None:
    lines = [
        "# Report Index",
        "",
        "Generated by:",
        "",
        "```bash",
        "python3 codexCrepeantResolutions/crepant_resolver.py",
        "```",
        "",
        "| Case | Status | Report |",
        "| --- | --- | --- |",
    ]
    for item in summary:
        case = item["case"]
        report_name = Path(item["report"]).name if str(item["report"]).endswith(".md") else ""
        report_link = f"[{report_name}]({report_name})" if report_name else item["report"]
        lines.append(f"| {case} | {item['status']} | {report_link} |")
    (out_dir / "INDEX.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data"))
    parser.add_argument("--out", type=Path, default=Path("codexCrepeantResolutions/reports"))
    parser.add_argument("--primes", type=int, nargs="*", default=list(DEFAULT_PRIMES))
    parser.add_argument("--bad-prime-bound", type=int, default=13)
    parser.add_argument("--octic", help="Macaulay2 toString output for an octic in QQ[x,y,z,w].")
    parser.add_argument("--octic-file", type=Path, help="File containing Macaulay2 toString output for an octic in QQ[x,y,z,w].")
    parser.add_argument("--case", default="octic", help="Report case name for --octic or --octic-file input.")
    parser.add_argument(
        "--variables",
        help="Comma-separated P3 variables for --octic/--octic-file input. Defaults to inferring from the polynomial.",
    )
    parser.add_argument("files", nargs="*", type=Path)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    summary = []

    if args.octic and args.octic_file:
        parser.error("use only one of --octic or --octic-file")
    if (args.octic or args.octic_file) and args.files:
        parser.error("--octic/--octic-file cannot be combined with positional files")

    if args.octic or args.octic_file:
        octic = args.octic if args.octic is not None else args.octic_file.read_text()
        variables = tuple(v.strip() for v in args.variables.split(",")) if args.variables else None
        result = analyze_octic_string(args.case, octic, tuple(args.primes), args.bad_prime_bound, variables)
        write_report(result, args.out)
        summary.append(
            {
                "case": result["case"],
                "status": result["status"],
                "path": result["path"],
                "report": str(args.out / f"{result['case']}.md"),
            }
        )
        (args.out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
        write_index(summary, args.out)
        print(f"{result['case']}: {result['status']} -> {args.out / (result['case'] + '.md')}")
        return 0

    files = args.files or discover_case_files(args.data_root)
    for path in files:
        result = analyze_file(path, tuple(args.primes), args.bad_prime_bound)
        write_report(result, args.out)
        summary.append(
            {
                "case": result["case"],
                "status": result["status"],
                "path": result["path"],
                "report": str(args.out / f"{result['case']}.md"),
            }
        )
        print(f"{result['case']}: {result['status']} -> {args.out / (result['case'] + '.md')}")

    if not args.files:
        try:
            from q1_paper_resolution import write_outputs as write_q1_paper_outputs

            write_q1_paper_outputs(args.out)
            summary.append(
                {
                    "case": "Q1 paper model",
                    "status": "crepant resolution",
                    "path": "q1_final_q5style_report.pdf",
                    "report": str(args.out / "Q1_paper_crepant_resolution.md"),
                }
            )
        except Exception as exc:
            summary.append(
                {
                    "case": "Q1 paper model",
                    "status": "error",
                    "path": "q1_final_q5style_report.pdf",
                    "report": f"Q1 paper resolution generation failed: {exc}",
                }
            )

    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    write_index(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
