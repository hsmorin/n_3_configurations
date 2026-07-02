-- Compare the existing repository Q1 branch octic with the paper's Q1 octic.
-- The output gives projective singular-locus component signatures.

load "q1_paper_model.m2"

paperComponents = singularComponents()
print("PAPER_COMPONENTS|" | toString(#paperComponents))
scan(#paperComponents, i -> (
    J := paperComponents#i;
    print("PAPER|" | toString(i) | "|" | toString(dim(R/J)-1) | "|" | toString(degree J) | "|" | toString(mingens J))
))

S = QQ[y,z,w,v]
Drepo = y^4*z^2*w^2-2*y^2*z^4*w^2+z^6*w^2+4*y^2*z^3*w^3-4*z^5*w^3-2*y^2*z^2*w^4+6*z^4*w^4-4*z^3*w^5+z^2*w^6-4*y^4*z^2*w*v+6*y^3*z^3*w*v-2*y^2*z^4*w*v+2*y*z^5*w*v-2*z^6*w*v-6*y^3*z^2*w^2*v-2*y*z^4*w^2*v+8*z^5*w^2*v+2*y^2*z^2*w^3*v-2*y*z^3*w^3*v-12*z^4*w^3*v+2*y*z^2*w^4*v+8*z^3*w^4*v-2*z^2*w^5*v+4*y^4*z^2*v^2-8*y^3*z^3*v^2+5*y^2*z^4*v^2-2*y*z^5*v^2+z^6*v^2+2*y^4*z*w*v^2+4*y^3*z^2*w*v^2-2*y*z^4*w*v^2-4*z^5*w*v^2+4*y^3*z*w^2*v^2-7*y^2*z^2*w^2*v^2+10*y*z^3*w^2*v^2+6*z^4*w^2*v^2+2*y^2*z*w^3*v^2-6*y*z^2*w^3*v^2-4*z^3*w^3*v^2+z^2*w^4*v^2-4*y^4*z*v^3+6*y^3*z^2*v^3-6*y^2*z^3*v^3+4*y*z^4*v^3-6*y^3*z*w*v^3+8*y^2*z^2*w*v^3-8*y*z^3*w*v^3-2*y^2*z*w^2*v^3+4*y*z^2*w^2*v^3+y^4*v^4
repoSing = saturate(ideal(Drepo) + ideal(diff(y,Drepo), diff(z,Drepo), diff(w,Drepo), diff(v,Drepo)), ideal(y,z,w,v))
repoComponents = decompose repoSing
print("REPO_COMPONENTS|" | toString(#repoComponents))
scan(#repoComponents, i -> (
    J := repoComponents#i;
    print("REPO|" | toString(i) | "|" | toString(dim(S/J)-1) | "|" | toString(degree J) | "|" | toString(mingens J))
))

paperSignature = sort apply(paperComponents, J -> (dim(R/J)-1, degree J))
repoSignature = sort apply(repoComponents, J -> (dim(S/J)-1, degree J))
print("SIGNATURE_EQUAL|" | toString(paperSignature == repoSignature))
