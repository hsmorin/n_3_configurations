R = QQ[x1,x2,x3,x4]

Q13 = x1^2*x2^2*x3-2*x1^2*x2*x3*x4+x1^2*x3*x4^2-x1^2*x2^2-x1^2*x2*x3+x1^2*x3^2-x1*x2*x3^2+x1^2*x2*x4+x1*x2^2*x4+2*x1^2*x3*x4-x1*x3^2*x4+x2*x3^2*x4-x1*x2*x4^2-2*x1*x3*x4^2+x2*x3*x4^2+x1^2*x2-x1^2*x3+x1*x2*x3-x1^2*x4-x1*x2*x4+x1*x3*x4-x2*x3*x4+x1*x4^2

M = matrix {{0, 1, 0, 1}, {1, 0, 0, 1}, {0, 0, 1, 0}, {0, 0, 0, 1}} -- T^{-1} change of variables

phi = map(R, R, flatten entries (M * transpose vars R))

f = phi Q13

-- f = x1^2*x2^2*x3+2*x1^2*x2*x3*x4+x1^2*x3*x4^2-x1^2*x2^2-x1*x2^2*x3-x1*x2*x3^2+x2^2*x3^2-x1^2*x2*x4-x1*x2^2*x4-2*x1*x2*x3*x4+x2^2*x3*x4-x1*x2*x4^2+x1*x2^2+x1*x2*x3-x2^2*x3+x1*x2*x4

-- Possible homogenizations: 
-- PP^1 \times PP^1 \times PP^1 \times PP^1
-- PP^1 \times PP^3 with either x1 or x3  sent to PP^1
