S = QQ[x0,x1,x2,x3,x4]

F = x0^3*x1*x2-x0^2*x1*x2^2-x0^3*x2*x3-2*x0^2*x1*x2*x3+2*x0^2*x2^2*x3+x0*x1*x2^2*x3+2*x0^2*x2*x3^2+x0*x1*x2*x3^2-3*x0*x2^2*x3^2-x0*x2*x3^3+x2^2*x3^3+x0^2*x1^2*x4-x0^2*x1*x2*x4-x0^2*x1*x3*x4-x0*x1^2*x3*x4+3*x0*x1*x2*x3*x4-2*x0*x2^2*x3*x4+x1*x2^2*x3*x4+x0*x1*x3^2*x4-2*x1*x2*x3^2*x4+x2^2*x3^2*x4-x0*x1^2*x4^2+2*x0*x1*x2*x4^2-x1^2*x2*x4^2+x1^2*x3*x4^2+x0*x2*x3*x4^2-2*x1*x2*x3*x4^2-x0*x1*x4^3+x1^2*x4^3

-- target P^3 coords y0,y1,y2,y3  <->  linear forms x0-x1, x0-x3, x2, x4  (vanish at P)
A = QQ[y0,y1,y2,y3]
B = A[t]

d0 = y0; d1 = 0; d2 = y2; d3 = y0-y1; d4 = y3;   -- representative direction Q(y)

phi = map(B,S,{1+t*d0, 1+t*d1, t*d2, 1+t*d3, t*d4})
Ft = phi F

-- confirm multiplicity 3 at P
print apply(0..2, i -> coefficient(t^i,Ft) == 0)   -- {true,true,true}

A3 = coefficient(t^3,Ft)
A4 = coefficient(t^4,Ft)
A5 = coefficient(t^5,Ft)

Disc = A4^2 - 4*A3*A5     -- discriminant of the quadratic g = A3 + A4 t + A5 t^2
factor Disc
