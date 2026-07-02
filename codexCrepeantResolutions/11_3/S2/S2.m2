-- Explicit double-octic derivation for S2 from 11_3.m2
-- Run with: M2 --script S2.m2
R = QQ[x0,x1,x2,x3,x4]
S2 = x1^2*x2*x3^2 + 2*x1^2*x2*x3*x4 - 2*x1^2*x2*x3 + x1^2*x2*x4^2 - 2*x1^2*x2*x4 + x1^2*x2 - x1^2*x3^2*x4 - x1^2*x3*x4^2 + 2*x1^2*x3*x4 + x1^2*x4^2 - x1^2*x4 + x1*x2^2*x3*x4^2 - x1*x2^2*x3*x4 - 2*x1*x2^2*x4^2 + x1*x2^2*x4 - x1*x2*x3^2*x4^2 + x1*x2*x3*x4 - x1*x2*x4^3 + 2*x1*x2*x4^2 - x1*x2*x4 + x1*x3^2*x4^2 + x1*x3*x4^3 - 3*x1*x3*x4^2 - 2*x1*x4^3 + 2*x1*x4^2 + x2^2*x4^3 - x2*x3*x4^3 + x3*x4^3 + x4^4 - x4^3
Q = S2
F = homogenize(Q, x0)
S = QQ[x,y,z,w,u,v]

-- Selected chart: x4 = 1
-- P1 coordinate: x1 -> x
-- P3 coordinates: x0,x2,x3 -> y,z,w
-- Common u-power removed from the discriminant: 2
phi = map(S, R, {y,x,z,w,1})
g = phi F
T = terms g
bad = select(T, term -> degree(x,term) > 2 or degree(y,term)+degree(z,term)+degree(w,term) > 4)
if #bad > 0 then error "selected chart is not bidegree (2,4)"
G = sum apply(T, term -> (
    n := 2 - degree(x, term);
    m := 4 - degree(y, term) - degree(z, term) - degree(w, term);
    term * u^n * v^m
))
A = (sum select(terms G, term -> degree(x, term) == 2)) // x^2
B = (sum select(terms G, term -> degree(x, term) == 1)) // x
C = sum select(terms G, term -> degree(x, term) == 0)
Disc = B^2 - 4*A*C
D = Disc // u^2
if degree(u,D) != 0 then error "discriminant still depends on u"
if not all(terms D, term -> degree(y,term)+degree(z,term)+degree(w,term)+degree(v,term) == 8) then error "D is not homogeneous octic in P3"
print("MODEL|equation|" | toString(Q))
print("MODEL|bidegree_equation_G|" | toString(G))
print("MODEL|quadratic_A|" | toString(A))
print("MODEL|linear_B|" | toString(B))
print("MODEL|constant_C|" | toString(C))
print("MODEL|octic_D|" | toString(D))
