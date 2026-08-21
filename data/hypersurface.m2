
hypersurface = I -> (
if not isPolynomialRing ring I then return I;
gens := flatten entries mingens I;
if (length gens == 1) or (length gens == 1) then return I;
usedVars := support I;
indSet := support first independentSets I;
depVars := toList(set usedVars - set indSet);
elimVars := drop(depVars, 1);
J := eliminate(elimVars, I);
J)
