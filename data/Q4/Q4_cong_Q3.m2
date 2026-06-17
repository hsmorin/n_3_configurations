R = QQ[x0,x1,x2,x3,x4]

Q3 = x1*x2^2*x3^2-2*x1^2*x2*x3*x4+x1^3*x4^2+x1*x2*x3*x4^2-x1^2*x4^3+x1^2*x2*x3-2*x1*x2^2*x3-x1*x2*x3^2-x1^3*x4+2*x1^2*x2*x4+x1^2*x3*x4+2*x1*x2*x3*x4-x2*x3^2*x4-x1^2*x4^2-x2*x3*x4^2+x1*x4^3-x1^2*x2+x1*x2^2+x2*x3^2+x1^2*x4-3*x1*x2*x4-x1*x3*x4+2*x2*x3*x4+x1*x2-x2*x3

Q4 = x1^2*x2*x3*x4+x2^2*x3^2*x4-x1^3*x4^2-2*x1*x2*x3*x4^2+x1^2*x4^3+x1*x2^2*x3-x2^2*x3^2+x1^3*x4-x1^2*x2*x4-x1^2*x3*x4-x2^2*x3*x4+x1^2*x4^2+x1*x2*x4^2+x2*x3*x4^2-x1*x4^3-x1*x2^2+x2^2*x3-x1^2*x4+x1*x2*x4+x1*x3*x4-x2*x3*x4

F = homogenize(Q3, x0)
G = homogenize(Q4, x0)

T = matrix {{1, 0, 0, 0, 0}, {1, 0, 0, 0, -1}, {1, 0, 0, -1, 0}, {1, 0, -1, 0, 0}, {1, -1, 0, 0, 0}} -- Linear Change of Variables

M = inverse T

phi = map(R, R, flatten entries (M * transpose vars R))

phi F == - G -- o : True therefore Z(F) \cong Z(G) in PP^4
