R = QQ[x0,x1,x2,x3,x4]

Q3 = x1*x2^2*x3^2-2*x1^2*x2*x3*x4+x1^3*x4^2+x1*x2*x3*x4^2-x1^2*x4^3+x1^2*x2*x3-2*x1*x2^2*x3-x1*x2*x3^2-x1^3*x4+2*x1^2*x2*x4+x1^2*x3*x4+2*x1*x2*x3*x4-x2*x3^2*x4-x1^2*x4^2-x2*x3*x4^2+x1*x4^3-x1^2*x2+x1*x2^2+x2*x3^2+x1^2*x4-3*x1*x2*x4-x1*x3*x4+2*x2*x3*x4+x1*x2-x2*x3

F = homogenize(Q3, x0)

f1 = F(x0,1,x2,x3,x4) -- Takes F to the first chart

f = sub(f1, {x3 => x0, x1 => x4, x0 => x3, x4 => x1}) -- Convient Change of Variables

-- f = -x0^2*x1*x2*x3-x0*x1^2*x2*x3+x0^2*x2*x3^2+2*x0*x1*x2*x3^2-x0*x2*x3^3+x0*x1^2*x2+x0^2*x2^2+x1^3*x3-x0^2*x2*x3+2*x0*x1*x2*x3-2*x0*x2^2*x3-x0*x1*x3^2-3*x1*x2*x3^2+x2^2*x3^2+x2*x3^3-x1^3-2*x0*x1*x2+x0*x1*x3-x1^2*x3+x0*x2*x3+2*x1*x2*x3+x1*x3^2-x2*x3^2+x1^2-x1*x3

-- Now f admits a compactification G in P_(x,u) \times P_(y,z,w,v)

-- Avoid reusing variables
S = QQ[x,y,z,w,u,v]
phi = map(S, R, {x0 => x, x1 => y, x2 => z, x3 =>w, x4 => 1})
g = phi f

T = terms g
TBar = apply(T, i -> (n = 2 - degree(x, i); m = 4 - degree(y, i) - degree(z, i) - degree(w,i); i * u^n * v^m))
G = sum TBar

G = -x*y^2*z*w*u+2*x*y*z*w^2*u-x*z*w^3*u+y^3*w*u^2-3*y*z*w^2*u^2+z^2*w^2*u^2+z*w^3*u^2-x^2*y*z*w*v+x^2*z*w^2*v+x*y^2*z*u*v+2*x*y*z*w*u*v-2*x*z^2*w*u*v-x*y*w^2*u*v-y^3*u^2*v-y^2*w*u^2*v+2*y*z*w*u^2*v+y*w^2*u^2*v-z*w^2*u^2*v+x^2*z^2*v^2-x^2*z*w*v^2-2*x*y*z*u*v^2+x*y*w*u*v^2+x*z*w*u*v^2+y^2*u^2*v^2-y*w*u^2*v^2

