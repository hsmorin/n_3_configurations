R = QQ[x1,x2,x3,x4]

Q12 = x1^2*x2^2*x3+x1^2*x2*x3^2-2*x1*x2^2*x3^2+x1^2*x2^2*x4-4*x1*x2^2*x3*x4+2*x2^2*x3^2*x4-x1^2*x2*x4^2-2*x1*x2^2*x4^2+2*x1*x2*x3*x4^2+3*x2^2*x3*x4^2-x2*x3^2*x4^2+2*x1*x2*x4^3+x2^2*x4^3-2*x2*x3*x4^3-x2*x4^4-x1^2*x2^2-x1^2*x2*x3+3*x1*x2^2*x3-x1^2*x3^2+2*x1*x2*x3^2-x2^2*x3^2+x1^2*x2*x4+3*x1*x2^2*x4-x1^2*x3*x4-3*x2^2*x3*x4-x1^2*x4^2-2*x1*x2*x4^2-2*x2^2*x4^2+x1*x3*x4^2+x2*x3*x4^2+x1*x4^3+x2*x4^3-x1*x2^2+x1^2*x3-2*x1*x2*x3+x2^2*x3+x1^2*x4+x2^2*x4-x1*x4^2

M = matrix{{0,1,0,1},{1,0,0,1},{0,0,1,0},{0,0,0,1}} -- Linear Change of Variables (M^{-1})

phi = map(R, R, flatten entries (M * transpose vars R))

f = phi Q12

degree f(x1, 1, 1, 1)

-- o: 1 => rational by projection away from x1
-- f = x1^2*x2^2*x3+x1^2*x2*x3^2-2*x1*x2^2*x3^2+x1^2*x2^2*x4-4*x1*x2^2*x3*x4+2*x2^2*x3^2*x4-x1^2*x2*x4^2-2*x1*x2^2*x4^2+2*x1*x2*x3*x4^2+3*x2^2*x3*x4^2-x2*x3^2*x4^2+2*x1*x2*x4^3+x2^2*x4^3-2*x2*x3*x4^3-x2*x4^4-x1^2*x2^2-x1^2*x2*x3+3*x1*x2^2*x3-x1^2*x3^2+2*x1*x2*x3^2-x2^2*x3^2+x1^2*x2*x4+3*x1*x2^2*x4-x1^2*x3*x4-3*x2^2*x3*x4-x1^2*x4^2-2*x1*x2*x4^2-2*x2^2*x4^2+x1*x3*x4^2+x2*x3*x4^2+x1*x4^3+x2*x4^3-x1*x2^2+x1^2*x3-2*x1*x2*x3+x2^2*x3+x1^2*x4+x2^2*x4-x1*x4^2 

-- Motivates homogenization in \P^3 \times \P^1
-- (x2, v) and (x1, x3, x4, u)

