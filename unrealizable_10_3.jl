import Combinatorics as comb
import Polymake
using Oscar


matroid = matroid_from_revlex_basis_encoding(0******0******0***0*******************************0*********************0*******0********************0*******0*******0**)
R = realization_space(matroid, ground_ring = QQ)
        display(R)
        I = gens(defining_ideal(R))
        L = inequations(R)
        I_str = "f$(i) = "
        L_str = "L$(i) = {"
        for f in I
            I_str *= string(f) * ", "
        end
        for g in L
            L_str *= string(g) * ", "
        end
        I_str = chop(I_str, tail=2)
        I_str *= " -- dimension "
        if is_realizable(R)
            I_str *= string(dim(R)) * "\n"
        else
            I_str *= "-∞\n"
        end
        I_str *= "-- revlex_basis_encoding: " * tag
        I_str *= "\n\n"
        L_str = chop(L_str, tail=2)
        L_str *= "}\n\n"
        write(n_3, I_str)
        write(n_3_L, L_str)
        print("\n")
end

close(n_3)
close(n_3_L)

