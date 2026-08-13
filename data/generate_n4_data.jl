import Combinatorics as comb
import Polymake
using Oscar

db = Polymake.Polydb.get_db()
collection = db["Matroids.Small"]


if length(ARGS) == 1
    n = parse(Int, ARGS[1])
else
    exit("Usage: julia generate_n4_data.jl n")
end

#Gets polymake data
results = Polymake.Polydb.find(collection, Dict("RANK" => 4, "N_ELEMENTS" => n, "SIMPLE" => true, "PAVING" => true, "N_MATROID_HYPERPLANES" => (binomial(n, 3) - 3 * n)))

db = nothing
collection = nothing
results = nothing
GC.gc()

csv_str = ""

#= Fill Code
for i in 0:10
    for j in 0:i
        con = true
        for l in H 
            if i ∈ l && j ∈ l
                con = false
            end
        end
        if con
            push!(H, [i,j])
        end
    end
end=#

#Finds 9_4 configurations from polymake database
lst = []
for (i, matroid) in enumerate(results)
    println("Evaluating matroid ${i}")
    tag = matroid.REVLEX_BASIS_ENCODING
    matroid = Matroid(matroid)
    H = hyperplanes(matroid)
    four_point_planes = filter(h -> length(h) == 4, Hs)
    three_point_planes = filter(h -> length(h) == 3, Hs)

    length(four_point_panes) == n || continue
    length(three_point_planes) == (binomial(n, 3) - 4 * n)

    E = matroid_groundset(M)
    degs = [count(h -> p in h, four_point_planes) for p in E]
    all(==(4), degs) || continue

    push!(lst, (matroid, tag))
    print("Matroid ${i} is valid")
end


n_4 = open("./data/$(n)_4.m2", "w")
n_4_L = open("./data/$n_4_avoids.m2", "w")

#Writes data to a file
for (i, (matroid, tag)) in enumerate(lst)
    R = realization_space(matroid, ground_ring=QQ)
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
    write(n_4, I_str)
    write(n_4_L, L_str)
    print("\n")
end

close(n_4)
close(n_4_L)

