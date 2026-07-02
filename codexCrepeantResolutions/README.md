# Crepant Resolution Diagnostics

Run:

```bash
python3 codexCrepeantResolutions/crepant_resolver.py
```

The script scans `data/Q*/*`, extracts branch octics written as `-- D = ...`,
uses Macaulay2 to decompose the singular locus, builds a first-level crepant
blowup candidate tree, and computes point counts for `t^2 = D` over the default
fields `F_3`, `F_5`, `F_7`, `F_11`, and `F_13`.

Outputs are written to `codexCrepeantResolutions/reports` as Markdown and JSON.

To analyze one octic pasted from Macaulay2 `toString` output in variables
`x,y,z,w`, run:

```bash
python3 codexCrepeantResolutions/crepant_resolver.py --case my_octic --octic 'x^8+y^8+z^8+w^8'
```

For longer octics, put the `toString` output in a text file and use:

```bash
python3 codexCrepeantResolutions/crepant_resolver.py --case my_octic --octic-file octic.txt
```

If the octic uses variables other than `x,y,z,w`, either let the script infer
them or pass them explicitly:

```bash
python3 codexCrepeantResolutions/crepant_resolver.py --case my_octic --variables a,b,c,d --octic-file octic.txt
```

For the Q1 crepant resolution from `q1_final_q5style_report.pdf`, run:

```bash
M2 --script codexCrepeantResolutions/q1_compare_models.m2
python3 codexCrepeantResolutions/q1_paper_resolution.py --out codexCrepeantResolutions/11_3/Q1 --case Q1
```

This writes `11_3/Q1/Q1.md`, `.json`, and `.m2` from the projection away from
`[1:1:0:1:0]`. The JSON contains the expanded affine chart tree, the projective
small-resolution data, and the resolved point-count correction
`#X1(F_p) = #Y0(F_p) + 58p^2 + 82p`. The comparison script shows that the
paper's branch octic in `P^3_{a,b,c,d}` is not projectively isomorphic to the
older octic extracted from `data/Q1/K3_Fibration.m2`: their singular-locus
component signatures differ.

For every named equation in `data/11_3.m2`, run:

```bash
cd codexCrepeantResolutions
python3 analyze_11_3.py
```

This uses the project-report projection for Q1 by default. For the other named
equations it searches projective affine charts and P1 variables for bidegree
`(2,4)` models in `P^1 x P^3`, extracts discriminant octics, and reuses the
crepant diagnostics. Pass `--q1-search` to force the old automatic Q1 chart
search. Outputs are written to `11_3/<equation-name>/`. Each folder
contains:

- `<name>.m2`: the explicit Macaulay2 derivation of the selected double-octic
  model.
- `<name>.md` and `<name>.json`: the singular-locus analysis, blowup-center
  ledger, singular-model point counts for `u^2 = D`, and first-level blowup
  count correction.
- `<name>_candidates.json` and, when analysis was attempted,
  `<name>_candidate_diagnostics.json`.

Important limitations:

- The tree is a crepant-candidate tree from the singular locus of the branch
  octic. Full certification needs strict-transform chart iteration after each
  blowup.
- The Picard rank is reported as a lower bound, `1 +` the first-level
  exceptional-divisor count.
- Bad reduction is detected by comparing singular-locus component signatures
  modulo small primes; this is a practical test, not a proof of the complete
  bad-prime set.
- Singular Kulikov fiber classification is not automatic yet because the
  repository does not expose one uniform K3 fibration interface for all cases.
