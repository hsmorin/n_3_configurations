-- Explicit double-octic derivation for F2 from 11_3.m2
-- Run with: M2 --script F2.m2
R = QQ[x0,x1,x2,x3,x4]
F2 = 0
Q = F2
F = homogenize(Q, x0)
S = QQ[x,y,z,w,u,v]

-- No bidegree (2,4) double-octic model was selected.
-- The search requires a projective affine chart and a quadratic P1 coordinate
-- whose homogenized equation has discriminant independent of u after removing a common u-power.
