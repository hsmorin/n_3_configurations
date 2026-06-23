R = QQ[x0,x1,x2,x3,x4]

Q2 = -x1^3*x4^2 + x1^3*x4 + x1^2*x2*x3*x4 - x1^2*x2*x3 + x1^2*x2*x4 + x1^2*x4^3 - x1^2*x4^2 - x1^2*x4 - x1*x2^2*x3*x4 - x1*x2^2*x3 + x1*x2^2*x4 - x1*x2*x3*x4^2 + x1*x2*x3*x4 + 2*x1*x2*x3 - 2*x1*x2*x4 + x1*x3*x4^2 - x1*x3*x4 - x1*x4^3 + 2*x1*x4^2 + x2^3*x3^2 - x2^3*x3 - x2^2*x3^2 + 2*x2^2*x3 + x2*x3*x4^2 - x2*x3*x4 - x2*x3 - x3*x4^2 + x3*x4

F = homogenize(Q2, x0)

f4 = F(x0,x1,x2,x3,1) -- Takes F to the fourth chart

-- Now f admits a compactification G in P_(x,u) \times P_(y,z,w,v)

-- Avoid reusing variables
S = QQ[x,y,z,w,u,v]
phi = map(S, R, {x3 => x, x0 => y, x1 => z, x2 => w, x4 => 1})
g = phi f4

-- Compactify g with bidegree(2,4) in P^1_(x,u) \times P^3_(y,z,w,v)
T = terms g
TBar = apply(T, i -> (n = 2 - degree(x, i); m = 4 - degree(y, i) - degree(z, i) - degree(w,i); i * u^n * v^m))
G = sum TBar

-- G = -x*y^3*w*u+2*x*y^2*z*w*u-x*y*z^2*w*u+2*x*y^2*w^2*u-x*y*z*w^2*u-x*y*w^3*u-y^2*z^2*u^2+y*z^3*u^2-2*y^2*z*w*u^2+y*z^2*w*u^2+y*z*w^2*u^2-x^2*y*w^2*v+x^2*w^3*v+x*y^3*u*v-x*y^2*z*u*v-x*y^2*w*u*v+x*y*z*w*u*v+x*z^2*w*u*v-x*z*w^2*u*v+2*y^2*z*u^2*v-y*z^2*u^2*v-z^3*u^2*v-x*y^2*u*v^2+x*y*z*u*v^2+x*y*w*u*v^2-x*z*w*u*v^2-y*z*u^2*v^2+z^2*u^2*v^2

-- Compute octic

a = (sum select(terms G, m -> degree(x, m) == 2)) / x^2
b = (sum select(terms G, m -> degree(x, m) == 1)) / x
c = sum select(terms G, m -> degree(x, m) == 0)

D = b^2 - 4*a*c
D = D / u^2
D = sub(D, S)

-- D = y^6*w^2-4*y^5*z*w^2+6*y^4*z^2*w^2-4*y^3*z^3*w^2+y^2*z^4*w^2-4*y^5*w^3+10*y^4*z*w^3-8*y^3*z^2*w^3+2*y^2*z^3*w^3+6*y^4*w^4-8*y^3*z*w^4+3*y^2*z^2*w^4-4*y^3*w^5+2*y^2*z*w^5+y^2*w^6-2*y^6*w*v+6*y^5*z*w*v-6*y^4*z^2*w*v+2*y^3*z^3*w*v+6*y^5*w^2*v-12*y^4*z*w^2*v+2*y^3*z^2*w^2*v+6*y^2*z^3*w^2*v-2*y*z^4*w^2*v-6*y^4*w^3*v+2*y^3*z*w^3*v+6*y^2*z^2*w^3*v-4*y*z^3*w^3*v+2*y^3*w^4*v+6*y^2*z*w^4*v-4*y*z^2*w^4*v-2*y*z*w^5*v+y^6*v^2-2*y^5*z*v^2+y^4*z^2*v^2-2*y^4*z*w*v^2+6*y^3*z^2*w*v^2-4*y^2*z^3*w*v^2-5*y^4*w^2*v^2+16*y^3*z*w^2*v^2-11*y^2*z^2*w^2*v^2+z^4*w^2*v^2+6*y^3*w^3*v^2-14*y^2*z*w^3*v^2+4*y*z^2*w^3*v^2+2*z^3*w^3*v^2-2*y^2*w^4*v^2+2*y*z*w^4*v^2+z^2*w^4*v^2-2*y^5*v^3+4*y^4*z*v^3-2*y^3*z^2*v^3+4*y^4*w*v^3-8*y^3*z*w*v^3+2*y^2*z^2*w*v^3+2*y*z^3*w*v^3-2*y^3*w^2*v^3+2*y^2*z*w^2*v^3+2*y*z^2*w^2*v^3-2*z^3*w^2*v^3+2*y*z*w^3*v^3-2*z^2*w^3*v^3+y^4*v^4-2*y^3*z*v^4+y^2*z^2*v^4-2*y^3*w*v^4+4*y^2*z*w*v^4-2*y*z^2*w*v^4+y^2*w^2*v^4-2*y*z*w^2*v^4+z^2*w^2*v^4

-- Prime

