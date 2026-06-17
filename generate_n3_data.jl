import Combinatorics as comb
import Polymake
using Oscar

db = Polymake.Polydb.get_db()
collection = db["Matroids.Small"]


if length(ARGS) == 1
    n = parse(Int, ARGS[1])
else
    exit("Usage: julia get_n_3_data.jl n")
end

#Gets polymake data
results = Polymake.Polydb.find(collection, Dict("RANK" => 3, "N_ELEMENTS" => n, "SIMPLE" => true))

matroid_list = collect(results)
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

#Finds n_3 configurations from polymake database
lst = []
for matroid in matroid_list
    tag = matroid.REVLEX_BASIS_ENCODING
    matroid = Matroid(matroid)
    H = hyperplanes(matroid)
    three = 0
    two = 0 
    cond = true
    for l in H
        if length(l) == 2
            two += 1
        elseif length(l) == 3
            three += 1
        else
            cond = false
            break
        end
    end

    if cond == false
        continue
    end

    if three == n && two == binomial(n,2) - 3 * n
        points = Dict()
        for l in H
            if length(l) == 3
                for p in l
                    points[p] = get(points, p, 0) + 1
                end
            end
        end
        con = true
        count = 0
        for item in points
            p = item[1]
            c = item[2]
            if c == 3
                count += 1
            else
                con = false
            end
        end
        
        if con && count == n
            push!(lst, (matroid, tag))
        end
    end
end


n_3 = open("./data/$(n)_3.m2", "w")
n_3_L = open("./data/$(n)_3_avoids.m2", "w")

#Writes data to a file
for (i, (matroid, tag)) in enumerate(lst)
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

