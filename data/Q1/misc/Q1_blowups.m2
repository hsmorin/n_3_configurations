R = QQ[x0,x1,x2,x3,x4]
S = QQ[x, y, z, w]

Q1 = x2^2*x3^3+x1*x2^2*x3*x4-2*x1*x2*x3^2*x4+x2^2*x3^2*x4-x1^2*x2*x4^2+x1^2*x3*x4^2-2*x1*x2*x3*x4^2+x1^2*x4^3+x1*x2^2*x3+x1*x2*x3^2-3*x2^2*x3^2-x2*x3^3-x1^2*x3*x4+3*x1*x2*x3*x4-2*x2^2*x3*x4+x1*x3^2*x4-x1^2*x4^2+2*x1*x2*x4^2+x2*x3*x4^2-x1*x4^3-x1*x2^2-2*x1*x2*x3+2*x2^2*x3+2*x2*x3^2+x1^2*x4-x1*x2*x4-x1*x3*x4+x1*x2-x2*x3

F = homogenize(Q1, x0)

-- Chart 1: q1, q2, q3, q8

f1 = F(x,y,z,w)

-- move q3 = (1,0,1,0) to 0

g = sub(f1, {x => x - 1, z => z - 1})

min(apply(terms g, i -> degree(i))) -- o: 3 so multiplicity 3

-- Blowup! \A^4 \times \P^3 has 5 charts?
