R = QQ[x1,x2,x3,x4]

Q5 = x1*x2^2*x3^2-x1^2*x2*x3*x4+x1^2*x2*x3-2*x1*x2^2*x3-x1*x2*x3^2+x2^2*x3^2+2*x1^2*x2*x4+x2^2*x3*x4-x2*x3^2*x4-x1^2*x4^2-x1*x2*x4^2+x1*x3*x4^2-x2*x3*x4^2+x1*x4^3-x1^2*x2+x1*x2^2+x1*x2*x3-x2^2*x3-x1*x2*x4+x2*x3*x4

S = QQ[x,y,z,w,u,v]

phi = map(S, R, {x2 => x, x1 => y, x3 => z, x4 => w})

f = phi Q5

-- Homogenize f into P^1 \times P^3
T = terms f
TBar = apply(T, i -> (n = 2 - degree(x, i); m = 4 - degree(y, i) - degree(z, i) - degree(w, i); i * u^n * v^m))
F = sum TBar

-- Bidegree (2,4)
-- F = -x*y^2*z*w*u-y^2*w^2*u^2+y*z*w^2*u^2+y*w^3*u^2+x^2*y*z^2*v+x*y^2*z*u*v-x*y*z^2*u*v+2*x*y^2*w*u*v-x*z^2*w*u*v-x*y*w^2*u*v-x*z*w^2*u*v-2*x^2*y*z*v^2+x^2*z^2*v^2+x^2*z*w*v^2-x*y^2*u*v^2+x*y*z*u*v^2-x*y*w*u*v^2+x*z*w*u*v^2+x^2*y*v^3-x^2*z*v^3 

-- Analyze Singular Fibres


-- Compute Octic

a = (sum select(terms F, m -> degree(x, m) == 2)) / x^2
b = (sum select(terms F, m -> degree(x, m) == 1)) / x
c = sum select(terms F, m -> degree(x, m) == 0)

D = (b^2 - 4 * a * c) / u^2
D = sub(D, S)

-- D = y^4*z^2*w^2-2*y^4*z^2*w*v+2*y^3*z^3*w*v-4*y^4*z*w^2*v+4*y^3*z^2*w^2*v-2*y^2*z^3*w^2*v+2*y^3*z*w^3*v-2*y^2*z^2*w^3*v+y^4*z^2*v^2-2*y^3*z^3*v^2+y^2*z^4*v^2+6*y^4*z*w*v^2-6*y^3*z^2*w*v^2-2*y^2*z^3*w*v^2+2*y*z^4*w*v^2+4*y^4*w^2*v^2-8*y^3*z*w^2*v^2+6*y^2*z^2*w^2*v^2-2*y*z^3*w^2*v^2+z^4*w^2*v^2-4*y^3*w^3*v^2+8*y^2*z*w^3*v^2-6*y*z^2*w^3*v^2+2*z^3*w^3*v^2+y^2*w^4*v^2-2*y*z*w^4*v^2+z^2*w^4*v^2-2*y^4*z*v^3+4*y^3*z^2*v^3-2*y^2*z^3*v^3-4*y^4*w*v^3+2*y^3*z*w*v^3+6*y^2*z^2*w*v^3-4*y*z^3*w*v^3+2*y^3*w^2*v^3-4*y^2*z*w^2*v^3+4*y*z^2*w^2*v^3-2*z^3*w^2*v^3-2*y^2*w^3*v^3+4*y*z*w^3*v^3-2*z^2*w^3*v^3+y^4*v^4-2*y^3*z*v^4+y^2*z^2*v^4+2*y^3*w*v^4-4*y^2*z*w*v^4+2*y*z^2*w*v^4+y^2*w^2*v^4-2*y*z*w^2*v^4+z^2*w^2*v^4
