R = QQ[x1,x2,x3,x4]

Q5 = x1*x2^2*x3^2-x1^2*x2*x3*x4+x1^2*x2*x3-2*x1*x2^2*x3-x1*x2*x3^2+x2^2*x3^2+2*x1^2*x2*x4+x2^2*x3*x4-x2*x3^2*x4-x1^2*x4^2-x1*x2*x4^2+x1*x3*x4^2-x2*x3*x4^2+x1*x4^3-x1^2*x2+x1*x2^2+x1*x2*x3-x2^2*x3-x1*x2*x4+x2*x3*x4

S = QQ[x,y,z,w,u,v,r,s,t]

phi = map(S, R, {x1 => x, x2 => y, x3 => z, x4 => w})

f = phi Q5

-- Homogenize f into P^1 \times P^3
T = terms f
TBar = apply(T, i -> (n = 2 - degree(x, i); m = 4 - degree(y, i) - degree(z, i) - degree(w, i); i * u^n * v^m))
F = sum TBar

-- Bidegree (2,4)
-- F = x*y^2*z^2*u+y^2*z^2*u^2+y^2*z*w*u^2-y*z^2*w*u^2-y*z*w^2*u^2-x^2*y*z*w*v-2*x*y^2*z*u*v-x*y*z^2*u*v-x*y*w^2*u*v+x*z*w^2*u*v+x*w^3*u*v-y^2*z*u^2*v+y*z*w*u^2*v+x^2*y*z*v^2+2*x^2*y*w*v^2-x^2*w^2*v^2+x*y^2*u*v^2+x*y*z*u*v^2-x*y*w*u*v^2-x^2*y*v^3

-- Analyze Singular Fibres


-- Compute Octic

a = (sum select(terms F, m -> degree(x, m) == 2)) / x^2
b = (sum select(terms F, m -> degree(x, m) == 1)) / x
c = sum select(terms F, m -> degree(x, m) == 0)

D = (b^2 - 4 * a * c) / u^2
D = sub(D, S)
D = D(1,x,y,z,1,w,1,1,1)
-- D = x^4*y^4-4*x^4*y^3*w-2*x^3*y^4*w+4*x^3*y^3*z*w+2*x^3*y^2*z^2*w-2*x^2*y^3*z^2*w-2*x^2*y^2*z^3*w+6*x^4*y^2*w^2+2*x^3*y^3*w^2+x^2*y^4*w^2-18*x^3*y^2*z*w^2+4*x^2*y^3*z*w^2-4*x^3*y*z^2*w^2+18*x^2*y^2*z^2*w^2-2*x*y^3*z^2*w^2+8*x^2*y*z^3*w^2-6*x*y^2*z^3*w^2+x^2*z^4*w^2-6*x*y*z^4*w^2+y^2*z^4*w^2-2*x*z^5*w^2+2*y*z^5*w^2+z^6*w^2-4*x^4*y*w^3+2*x^3*y^2*w^3-2*x^2*y^3*w^3+16*x^3*y*z*w^3-6*x^2*y^2*z*w^3-2*x^3*z^2*w^3-16*x^2*y*z^2*w^3+2*x*y^2*z^2*w^3+4*x^2*z^3*w^3+4*x*y*z^3*w^3-2*x*z^4*w^3+x^4*w^4-2*x^3*y*w^4+x^2*y^2*w^4-2*x^3*z*w^4+2*x^2*y*z*w^4+x^2*z^2*w^4 
-- Prime

-- Singular Locus
-- L = {ideal(w,y), ideal(w,x), ideal(z,y-w), ideal(z,x), ideal(z-w,x-w), ideal(y,x-z), ideal(y-w,x-z), ideal(y+z-w,x-w)}

-- L1 = ideal(w,y)
-- L2 = ideal(w,x)
-- L3 = ideal(z,y-w)
-- L4 = ideal(z,x)
-- L5 = ideal(z-w,x-w)
-- L6 = ideal(y,x-z)
-- L7 = ideal(y-w,x-z)
-- L8 = ideal(y+z-w,x-w)

-- Intersection Points
-- P1 = [0 : 0 : 1 : 0] -- L1, L2 mult 2
-- P2 = [1 : 0 : 0 : 0] -- L1, L3 mult 4
-- P3 = [0 : 1 : 0 : 0] -- L2, L3, L5 mult 4 - chart 2
-- P4 = [0 : 1 : 0 : 1] -- L3, L4, L7 mult 4 - chart 4
-- P5 = [1 : 0 : 1 : 0] -- L1, L6, L7 mult 4 - chart 3
-- P6 = [0 : 0 : 0 : 1] -- L4, L6 mult 4 - chart 4
-- P7 = [1 : 0 : 1 : 1] -- L5, L6, L8 mult 4 - chart 4
-- P8 = [1 : 1 : 1 : 1] -- L5, L7 mult 3 - chart 4
-- P9 = [0 : 1 : -1 : 0] -- L2, L8 mult 4 - chart 2
-- P10 = [1 : 1 : 0 : 1] -- L3, L8 mult 2


-- Blowup order
-- P6

-- Chart 1
d1 = D(1,y,z,w,1,1)


-- Chart 2
d2 = D(x,1,z,w,1,1)


-- Chart 3
d3 = D(x,y,1,w,1,1)


-- Chart 4
d4 = D(x,y,z,1,1,1)

-- P4 = (0, 1, 0) -- L3, L4, L7 mult 4
-- P6 = (0, 0, 0) -- L4, L6 mult 4
-- P7 = (1, 0, 1) -- L5, L6, L8 mult 4
-- P8 = (1, 1, 1) -- L5, L7 mult 3
-- P10 = (1, 1, 0) -- L3, L8 mult 2

-- 0 has multiplicity 4, blowup P6 (x,y,z)[r,s,t]

d41 = (mingens (decompose eliminate({y,z}, ideal(d4, x*s - y, x*t - z, y*t - z*s)))#1)_(0,0)
d41 = d41(x,y,z,w,u,v,x,y,z)
-- P7 = (1,0,1) mult 4
-- P8 = (1,1,1) mult 3
-- P10 = (1, 1, 0)  mult 2

-- Move P7 to 0
d41 = d41(x+1,y,z+1,w,1,1,1,1,1)

-- Blowup P7 = 0
d411 = (mingens (decompose eliminate({y,z}, ideal(d41, x*s - y, x*t - z, y*t - z*s)))#1)_(0,0)
d411 = d411(x,y,z,w,u,v,x,y,z)

-- Move P8 to 0
d411 = d411(x+1,y+1,z+1,w,1,1,1,1,1)


d412 = (mingens (decompose eliminate({y,z}, ideal(d41, x*s - y, x*t - z, y*t - z*s)))#1)_(0,0)
d412 = d412(x,y,z,w,u,v,x,y,z)

d413 = (mingens (decompose eliminate({y,z}, ideal(d41, x*s - y, x*t - z, y*t - z*s)))#1)_(0,0)
d413 = d413(x,y,z,w,u,v,x,y,z)


d42 = (mingens (decompose eliminate({x,z}, ideal(d41, x - y*r, x*t - z*r, y*t - z)))#1)_(0,0)
-- P6 = (0,0,0) mult 2
-- P4 = (0,1,0)
-- P8 = (1,1,1)
-- P10 = (1,1,0)

d43 = (mingens (decompose eliminate({x,y}, ideal(d4, x*s - y*r, x - z*r, y - z*s)))#1)_(0,0)
-- P7 = (1,0,1) mult 4
-- P8 = (1,1,1) mult 3


