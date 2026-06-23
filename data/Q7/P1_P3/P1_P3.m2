R = QQ[x0,x1,x2,x3,x4]

Q7 = x1^3*x4^2 - x1^3*x4 - 2*x1^2*x2*x3*x4 + x1^2*x2*x3 - x1^2*x3*x4^2 + 2*x1^2*x3*x4 - x1^2*x4^3 + x1^2*x4^2 + x1*x2^2*x3^2 + x1*x2^2*x3*x4 - x1*x2^2*x3 - x1*x2^2*x4 + x1*x2^2 + x1*x2*x3^2*x4 - 2*x1*x2*x3^2 + x1*x2*x3*x4^2 - x1*x2*x3*x4 + x1*x2*x3 + x1*x2*x4 - x1*x2 + x1*x3*x4^2 - 2*x1*x3*x4 + x1*x4^3 - 2*x1*x4^2 + x1*x4 - x2^3*x3^2 + x2^3*x3 + x2^2*x3^2 - 2*x2^2*x3 - x2*x3^2*x4 + x2*x3^2 - x2*x3*x4^2 + 2*x2*x3*x4 -- dimension 3


F = homogenize(Q7, x0)

f2 = F(x0,x1,1,x3,x4) -- Takes F to the second chart

-- Now f admits a compactification G in P_(x,u) \times P_(y,z,w,v)

-- Avoid reusing variables
S = QQ[x,y,z,w,u,v]
phi = map(S, R, {x3 => x, x0 => y, x2 => z, x4 => w, x1 => 1})
g = phi f2

-- Compactify g with bidegree(2,4) in P^1_(x,u) \times P^3_(y,z,w,v)
T = terms g
TBar = apply(T, i -> (n = 2 - degree(x, i); m = 4 - degree(y, i) - degree(z, i) - degree(w,i); i * u^n * v^m))
G = sum TBar

-- G = y^3*w*u^2-2*y^2*w^2*u^2+y*w^3*u^2-y^3*u^2*v+y^2*w*u^2*v+y*w^2*u^2*v-w^3*u^2*v+x^2*y^2*v^2-x^2*y*w*v^2-x*y^2*u*v^2+x*y*w*u*v^2+y^2*u^2*v^2-2*y*w*u^2*v^2+w^2*u^2*v^2-x^2*y*v^3+x^2*w*v^3+x*y*u*v^3-x*w*u*v^3

-- Compute octic

a = (sum select(terms G, m -> degree(x, m) == 2)) / x^2
b = (sum select(terms G, m -> degree(x, m) == 1)) / x
c = sum select(terms G, m -> degree(x, m) == 0)

D = b^2 - 4*a*c
D = D / u^2
D = sub(D, S)

-- D = -4*y^5*w*v^2+12*y^4*w^2*v^2-12*y^3*w^3*v^2+4*y^2*w^4*v^2+4*y^5*v^3-4*y^4*w*v^3-12*y^3*w^2*v^3+20*y^2*w^3*v^3-8*y*w^4*v^3-7*y^4*v^4+18*y^3*w*v^4-11*y^2*w^2*v^4-4*y*w^3*v^4+4*w^4*v^4+2*y^3*v^5-8*y^2*w*v^5+10*y*w^2*v^5-4*w^3*v^5+y^2*v^6-2*y*w*v^6+w^2*v^6
-- Three planes (with mult 2) and a quadratic

