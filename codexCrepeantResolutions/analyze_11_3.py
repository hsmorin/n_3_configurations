#!/usr/bin/env python3
"""
Find double-octic models for the named equations in data/11_3.m2 and analyze them.

For each affine polynomial Q_i(x1,x2,x3,x4), the script homogenizes with x0,
searches all projective affine charts x_j = 1 and all choices of one remaining
coordinate as the P^1 coordinate, compactifies to bidegree (2,4) in
P^1_(x,u) x P^3_(y,z,w,v), and extracts the discriminant octic when it becomes
independent of u after removing a common u-power.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from crepant_resolver import (
    blowup_count_corrections,
    count_double_cover,
    detected_bad_primes,
    infer_blowup_tree,
    m2_singular_profile,
    parse_polynomial,
    primes_upto,
    run_m2,
    write_index,
    write_report,
)
from q1_paper_resolution import DEFAULT_COUNT_PRIMES as Q1_DEFAULT_COUNT_PRIMES
from q1_paper_resolution import write_case_outputs as write_q1_case_outputs


PROJECTIVE_VARS = ("x0", "x1", "x2", "x3", "x4")


def extract_q_polynomials(path: Path) -> dict[str, str]:
    text = path.read_text()
    out: dict[str, str] = {}
    for match in re.finditer(r"^([A-Z]\d+)\s*=\s*(.+?)\s*--\s*dimension\s+3\s*$", text, re.MULTILINE):
        name, expr = match.groups()
        out[name] = expr.strip()
    return dict(sorted(out.items(), key=lambda kv: (kv[0][0], int(kv[0][1:]))))


def candidate_script(name: str, expr: str) -> str:
    blocks = [
        "R = QQ[x0,x1,x2,x3,x4]",
        f"Q = {expr}",
        "F = homogenize(Q, x0)",
        "S = QQ[x,y,z,w,u,v]",
    ]
    serial = 0
    for chart in PROJECTIVE_VARS:
        remaining = [var for var in PROJECTIVE_VARS if var != chart]
        for qvar in remaining:
            serial += 1
            p3vars = [var for var in remaining if var != qvar]
            image = []
            for var in PROJECTIVE_VARS:
                if var == chart:
                    image.append("1")
                elif var == qvar:
                    image.append("x")
                elif var == p3vars[0]:
                    image.append("y")
                elif var == p3vars[1]:
                    image.append("z")
                elif var == p3vars[2]:
                    image.append("w")
                else:
                    raise AssertionError(var)
            label = f"c{serial}"
            blocks.append(
                f"""
-- candidate {chart} {qvar} {','.join(p3vars)}
phi{label} = map(S, R, {{{','.join(image)}}})
g{label} = phi{label} F
T{label} = terms g{label}
bad{label} = select(T{label}, term -> degree(x,term) > 2 or degree(y,term)+degree(z,term)+degree(w,term) > 4)
if #bad{label} == 0 then (
    G{label} = sum apply(T{label}, term -> (
        n := 2 - degree(x, term);
        m := 4 - degree(y, term) - degree(z, term) - degree(w, term);
        term * u^n * v^m
    ));
    A{label} = (sum select(terms G{label}, term -> degree(x, term) == 2)) // x^2;
    B{label} = (sum select(terms G{label}, term -> degree(x, term) == 1)) // x;
    C{label} = sum select(terms G{label}, term -> degree(x, term) == 0);
    Disc{label} = B{label}^2 - 4*A{label}*C{label};
    scan(0..12, k -> (
        if Disc{label} % ideal(u^k) == 0 then (
            Dk := Disc{label} // u^k;
            if degree(u,Dk) == 0 and all(terms Dk, term -> degree(y,term)+degree(z,term)+degree(w,term)+degree(v,term) == 8) then (
                print("CAND|{name}|chart={chart}|qvar={qvar}|p3={','.join(p3vars)}|uPower=" | toString(k) | "|terms=" | toString(#terms Dk) | "|D=" | toString(Dk));
            );
        );
    ));
);
"""
            )
    return "\n".join(blocks)


def find_candidates(name: str, expr: str, timeout: int) -> list[dict]:
    out = run_m2(candidate_script(name, expr), timeout=timeout)
    candidates = []
    for line in out.splitlines():
        if not line.startswith("CAND|"):
            continue
        parts = line.split("|")
        meta = {"case": parts[1]}
        for part in parts[2:]:
            key, value = part.split("=", 1)
            meta[key] = value
        candidates.append(meta)
    return candidates


def selected_derivation_script(name: str, expr: str, candidate: dict | None) -> str:
    header = [
        f"-- Explicit double-octic derivation for {name} from 11_3.m2",
        "-- Run with: M2 --script " + f"{name}.m2",
        "R = QQ[x0,x1,x2,x3,x4]",
        f"{name} = {expr}",
        f"Q = {name}",
        "F = homogenize(Q, x0)",
        "S = QQ[x,y,z,w,u,v]",
    ]
    if candidate is None:
        return "\n".join(
            header
            + [
                "",
                "-- No bidegree (2,4) double-octic model was selected.",
                "-- The search requires a projective affine chart and a quadratic P1 coordinate",
                "-- whose homogenized equation has discriminant independent of u after removing a common u-power.",
            ]
        ) + "\n"

    chart = candidate["chart"]
    qvar = candidate["qvar"]
    p3vars = candidate["p3"].split(",")
    image = []
    for var in PROJECTIVE_VARS:
        if var == chart:
            image.append("1")
        elif var == qvar:
            image.append("x")
        elif var == p3vars[0]:
            image.append("y")
        elif var == p3vars[1]:
            image.append("z")
        elif var == p3vars[2]:
            image.append("w")
        else:
            raise AssertionError(var)

    return "\n".join(
        header
        + [
            "",
            f"-- Selected chart: {chart} = 1",
            f"-- P1 coordinate: {qvar} -> x",
            f"-- P3 coordinates: {p3vars[0]},{p3vars[1]},{p3vars[2]} -> y,z,w",
            f"-- Common u-power removed from the discriminant: {candidate['uPower']}",
            f"phi = map(S, R, {{{','.join(image)}}})",
            "g = phi F",
            "T = terms g",
            "bad = select(T, term -> degree(x,term) > 2 or degree(y,term)+degree(z,term)+degree(w,term) > 4)",
            "if #bad > 0 then error \"selected chart is not bidegree (2,4)\"",
            "G = sum apply(T, term -> (",
            "    n := 2 - degree(x, term);",
            "    m := 4 - degree(y, term) - degree(z, term) - degree(w, term);",
            "    term * u^n * v^m",
            "))",
            "A = (sum select(terms G, term -> degree(x, term) == 2)) // x^2",
            "B = (sum select(terms G, term -> degree(x, term) == 1)) // x",
            "C = sum select(terms G, term -> degree(x, term) == 0)",
            "Disc = B^2 - 4*A*C",
            f"D = Disc // u^{candidate['uPower']}",
            "if degree(u,D) != 0 then error \"discriminant still depends on u\"",
            "if not all(terms D, term -> degree(y,term)+degree(z,term)+degree(w,term)+degree(v,term) == 8) then error \"D is not homogeneous octic in P3\"",
            'print("MODEL|equation|" | toString(Q))',
            'print("MODEL|bidegree_equation_G|" | toString(G))',
            'print("MODEL|quadratic_A|" | toString(A))',
            'print("MODEL|linear_B|" | toString(B))',
            'print("MODEL|constant_C|" | toString(C))',
            'print("MODEL|octic_D|" | toString(D))',
        ]
    ) + "\n"


def analyze_candidate(candidate: dict, primes: tuple[int, ...], bad_prime_bound: int, source_path: Path) -> dict:
    poly = parse_polynomial(candidate["D"], variables=("y", "z", "w", "v"))
    profile = m2_singular_profile(poly)
    tree = infer_blowup_tree(profile)
    bad_test_primes = tuple(primes_upto(bad_prime_bound))
    singular_counts = {str(p): count_double_cover(poly, p) for p in primes}
    corrections = blowup_count_corrections(poly, profile, primes)
    corrected_counts = {
        str(p): singular_counts[str(p)] + corrections[str(p)]["total_correction"]
        for p in primes
        if corrections[str(p)]["total_correction"] is not None
    }
    return {
        "case": candidate["case"],
        "path": str(source_path),
        "status": tree["status"],
        "has_octic": True,
        "model_search": {
            "candidate_count": None,
            "selected_candidate_index": None,
            "chart": candidate["chart"],
            "quadratic_variable": candidate["qvar"],
            "p3_variables_from_projective_coordinates": candidate["p3"],
            "u_power_removed_from_discriminant": int(candidate["uPower"]),
        },
        "variables": poly.variables,
        "term_count": len(poly.terms),
        "octic": candidate["D"],
        "singular_locus": profile,
        "blowup_tree": tree,
        "picard_rank": {
            "lower_bound": tree["picard_rank_lower_bound"],
            "note": "Computed from the first-level crepant-candidate centers found on the branch octic.",
        },
        "bad_reduction": {
            "detected_primes": detected_bad_primes(poly, profile, bad_test_primes),
            "tested_primes": list(bad_test_primes),
            "note": "2 is always included; odd primes are flagged when the singular-locus component signature changes modulo p.",
        },
        "finite_field_counts": singular_counts,
        "singular_model_finite_field_counts": singular_counts,
        "blowup_count_correction": {
            "status": "first_level_from_listed_blowup_centers",
            "note": "For a curve center the blowup replaces the center by a P^1-bundle, contributing p times the center count. For a point center in a threefold it contributes p+p^2. These are first-level center corrections from the ledger below.",
            "by_prime": corrections,
        },
        "corrected_first_level_finite_field_counts": corrected_counts,
        "singular_kulikov_fibers": {
            "status": "not_computed",
            "note": "This 11_3 batch constructs double-octic branch models only; no uniform K3 pencil was selected.",
        },
    }


def write_model_report(results: list[dict], out_dir: Path) -> None:
    lines = [
        "# 11_3 Double-Octic Model Search",
        "",
        "| Case | Status | Chart | Quadratic variable | P3 variables | Folder |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for result in results:
        case = result["case"]
        if result.get("has_octic"):
            search = result["model_search"]
            qvar = search.get("quadratic_variable", search.get("qvar", ""))
            p3vars = search.get("p3_variables_from_projective_coordinates", search.get("p3", ""))
            lines.append(
                f"| {case} | {result['status']} | `{search['chart']}` | `{qvar}` | "
                f"`{p3vars}` | [{case}]({case}/{case}.md) |"
            )
        else:
            lines.append(f"| {case} | error |  |  |  | [{case}]({case}/{case}.md) |")
    (out_dir / "INDEX.md").write_text("\n".join(lines) + "\n")


def write_error_report(case: str, source_path: Path, error: str, out_dir: Path) -> dict:
    result = {"case": case, "path": str(source_path), "status": "error", "has_octic": False, "error": error}
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{case}.json").write_text(json.dumps(result, indent=2, sort_keys=True))
    (out_dir / f"{case}.md").write_text(f"# {case}\n\nSource: `{source_path}`\n\nStatus: `error`\n\n{error}\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("data/11_3.m2"))
    parser.add_argument("--out", type=Path, default=Path("11_3"))
    parser.add_argument("--cases", nargs="*", help="Case names such as Q1 Q2. Defaults to every named equation.")
    parser.add_argument("--primes", type=int, nargs="*", default=[3, 5, 7])
    parser.add_argument("--bad-prime-bound", type=int, default=7)
    parser.add_argument("--m2-timeout", type=int, default=120)
    parser.add_argument("--find-only", action="store_true")
    parser.add_argument("--first-candidate-only", action="store_true")
    parser.add_argument(
        "--q1-search",
        action="store_true",
        help="Use the automatic 11_3 chart search for Q1 instead of the project-report projection from [1:1:0:1:0].",
    )
    args = parser.parse_args()

    source = args.source
    if not source.exists() and source == Path("data/11_3.m2") and Path("../data/11_3.m2").exists():
        source = Path("../data/11_3.m2")

    polys = extract_q_polynomials(source)
    if args.cases:
        polys = {case: polys[case] for case in args.cases}

    args.out.mkdir(parents=True, exist_ok=True)
    results = []
    summary = []
    for case, expr in polys.items():
        case_dir = args.out / case
        case_dir.mkdir(parents=True, exist_ok=True)
        if case == "Q1" and not args.q1_search:
            q1_primes = tuple(dict.fromkeys([*args.primes, *Q1_DEFAULT_COUNT_PRIMES]))
            write_q1_case_outputs(case_dir, "Q1", q1_primes)
            result = {
                "case": "Q1",
                "status": "crepant_resolution",
                "path": "q1_final_q5style_report.pdf",
                "has_octic": True,
                "model_search": {
                    "chart": "projection from [1:1:0:1:0]",
                    "quadratic_variable": "u",
                    "p3_variables_from_projective_coordinates": "a,b,c,d",
                },
            }
            results.append(result)
            summary.append({"case": case, "status": result["status"], "path": result["path"], "report": str(case_dir / f"{case}.md")})
            print(f"{case}: crepant_resolution -> {case_dir / (case + '.md')}")
            continue
        print(f"{case}: searching double-octic models")
        try:
            candidates = find_candidates(case, expr, args.m2_timeout)
            (case_dir / f"{case}_candidates.json").write_text(json.dumps(candidates, indent=2, sort_keys=True))
            if not candidates:
                (case_dir / f"{case}.m2").write_text(selected_derivation_script(case, expr, None))
                result = write_error_report(case, source, "No bidegree (2,4) double-octic model was found by the chart search.", case_dir)
            elif args.find_only:
                (case_dir / f"{case}.m2").write_text(selected_derivation_script(case, expr, candidates[0]))
                result = {
                    "case": case,
                    "path": str(source),
                    "status": "model_found",
                    "has_octic": True,
                    "candidate_count": len(candidates),
                    "model_search": candidates[0],
                }
                (case_dir / f"{case}.json").write_text(json.dumps(result, indent=2, sort_keys=True))
                (case_dir / f"{case}.md").write_text(f"# {case}\n\nFound {len(candidates)} candidate double-octic model(s).\n")
            else:
                result = None
                diagnostics = []
                candidates_to_try = candidates[:1] if args.first_candidate_only else candidates
                for idx, candidate in enumerate(candidates_to_try):
                    try:
                        candidate_result = analyze_candidate(candidate, tuple(args.primes), args.bad_prime_bound, source)
                        candidate_result["model_search"]["candidate_count"] = len(candidates)
                        candidate_result["model_search"]["selected_candidate_index"] = idx
                        diagnostics.append(
                            {
                                "index": idx,
                                "status": candidate_result["status"],
                                "chart": candidate["chart"],
                                "qvar": candidate["qvar"],
                                "p3": candidate["p3"],
                                "error_count": len(candidate_result["blowup_tree"]["errors"]),
                            }
                        )
                        if result is None or (
                            result["status"] != "candidate" and candidate_result["status"] == "candidate"
                        ):
                            result = candidate_result
                        if candidate_result["status"] == "candidate":
                            break
                    except Exception as exc:
                        diagnostics.append(
                            {
                                "index": idx,
                                "status": "error",
                                "chart": candidate["chart"],
                                "qvar": candidate["qvar"],
                                "p3": candidate["p3"],
                                "error": str(exc),
                            }
                        )
                (case_dir / f"{case}_candidate_diagnostics.json").write_text(json.dumps(diagnostics, indent=2, sort_keys=True))
                if result is None:
                    (case_dir / f"{case}.m2").write_text(selected_derivation_script(case, expr, candidates[0]))
                    result = write_error_report(case, source, "Every discovered candidate failed during singular-locus analysis.", case_dir)
                    results.append(result)
                    summary.append({"case": case, "status": "error", "path": str(source), "report": str(case_dir / f"{case}.md")})
                    print(f"{case}: error -> {case_dir / (case + '.md')}")
                    continue
                selected_idx = result["model_search"]["selected_candidate_index"]
                (case_dir / f"{case}.m2").write_text(selected_derivation_script(case, expr, candidates[selected_idx]))
                write_report(result, case_dir)
            results.append(result)
            summary.append({"case": case, "status": result["status"], "path": str(source), "report": str(case_dir / f"{case}.md")})
            print(f"{case}: {result['status']} -> {case_dir / (case + '.md')}")
        except Exception as exc:
            (case_dir / f"{case}.m2").write_text(selected_derivation_script(case, expr, None))
            result = write_error_report(case, source, str(exc), case_dir)
            results.append(result)
            summary.append({"case": case, "status": "error", "path": str(source), "report": str(case_dir / f"{case}.md")})
            print(f"{case}: error -> {exc}")

    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    write_model_report(results, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
