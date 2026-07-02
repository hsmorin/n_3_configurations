-- Q1 branch octic from q1_final_q5style_report.pdf.
-- Projection center: [1:1:0:1:0].
-- Projection substitution: X0=u+a, X1=u+b, X2=c, X3=u, X4=d.
-- Coordinates are P^3_{a,b,c,d}; the double cover is R^2 = Delta.

R = QQ[a,b,c,d]

q3 = a*d*(b + c - d)

q4 = (
     a^2*b*c + a^2*b*d + a^2*c^2
     - a^2*c*d + a*b^2*d - a*b*c^2
     + a*b*c*d - 2*a*b*d^2 - 2*a*c^2*d
     + 3*a*c*d^2 - a*d^3 + b*c^2*d
     - 2*b*c*d^2 + b*d^3
)

q5 = (
     a^3*b*c + a^2*b^2*d - a^2*b*c^2
     - a^2*b*c*d - a*b^2*d^2 + 2*a*b*c*d^2
     - a*b*d^3 - b^2*c*d^2 + b^2*d^3
)

Delta = q4^2 - 4*q3*q5
-- D = (a^2*b*c+a^2*b*d+a^2*c^2-a^2*c*d+a*b^2*d-a*b*c^2+a*b*c*d-2*a*b*d^2-2*a*c^2*d+3*a*c*d^2-a*d^3+b*c^2*d-2*b*c*d^2+b*d^3)^2-4*(a*d*(b+c-d))*(a^3*b*c+a^2*b^2*d-a^2*b*c^2-a^2*b*c*d-a*b^2*d^2+2*a*b*c*d^2-a*b*d^3-b^2*c*d^2+b^2*d^3)

singularIdeal = () -> (
    irr := ideal(a,b,c,d);
    saturate(ideal(Delta) + ideal(diff(a,Delta), diff(b,Delta), diff(c,Delta), diff(d,Delta)), irr)
)

singularComponents = () -> decompose singularIdeal()
