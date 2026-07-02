-- Explicit double-octic derivation for V2 from 11_3.m2
-- Run with: M2 --script V2.m2
R = QQ[x0,x1,x2,x3,x4]
V2 = x1^2*x2^2*x3^2 - 2*x1^2*x2*x3^2*x4 + x1^2*x2*x3*x4^2 - x1^2*x2*x3*x4 + x1^2*x3^2*x4^2 - x1^2*x3*x4^3 + x1^2*x3*x4^2 + x1*x2^2*x3^2*x4^2 - 3*x1*x2^2*x3^2*x4 - 2*x1*x2^2*x3*x4^2 + 2*x1*x2^2*x3*x4 + 2*x1*x2*x3^2*x4^2 + 2*x1*x2*x3^2*x4 + x1*x2*x3*x4^4 - 2*x1*x2*x3*x4^3 + 3*x1*x2*x3*x4^2 - 2*x1*x2*x3*x4 - x1*x2*x4^4 + 2*x1*x2*x4^3 - x1*x2*x4^2 - x1*x3^2*x4^3 - x1*x3^2*x4^2 + x1*x4^4 - 2*x1*x4^3 + x1*x4^2 - x2^2*x3^2*x4^3 + 2*x2^2*x3^2*x4^2 - x2^2*x3*x4^4 + 4*x2^2*x3*x4^3 - 3*x2^2*x3*x4^2 + x2^2*x4^4 - 2*x2^2*x4^3 + x2^2*x4^2 - 2*x2*x3^2*x4^2 - 3*x2*x3*x4^3 + 3*x2*x3*x4^2 - x2*x4^4 + 2*x2*x4^3 - x2*x4^2 + x3^2*x4^3 + x3*x4^4 - x3*x4^3
Q = V2
F = homogenize(Q, x0)
S = QQ[x,y,z,w,u,v]

-- No bidegree (2,4) double-octic model was selected.
-- The search requires a projective affine chart and a quadratic P1 coordinate
-- whose homogenized equation has discriminant independent of u after removing a common u-power.
