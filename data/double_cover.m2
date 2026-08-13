
doubleCover = f -> (
    S := QQ[y0,y1,y2,y3,y4];
    R := ring f;
    phi := map(S, R, vars S);
    F := homogenize(phi f, y0);
    degTwoVars := select(flatten entries vars S, v -> degree(v, F) === 2);
    models = {};
    for var in degTwoVars do (
        A := (sum select(terms F, m -> degree(var, m) == 2)) / var^2;
        B := (sum select(terms F, m -> degree(var, m) == 1)) / var;
        C := sum select(terms F, m -> degree(var, m) == 0);
        D := B^2 - 4*A*C;
        D = sub(D, S);
        for f in (toList factor D) do (
            m = f#1;
            f = f#0;
            if (m % 2) == 0 then (
                D = D / (f^m);
                );
            );
        models = models | {sub(D, S)};
        );
    psi := map(R, S, vars R);
    apply(models, D -> psi D)
)
