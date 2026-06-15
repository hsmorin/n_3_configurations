import Combinatorics as comb
import Polymake
using Oscar

db = Polymake.Polydb.get_db()
collection = db["Matroids.Small"]

n = 10
results = Polymake.Polydb.find(collection, Dict("RANK" => 3, "N_ELEMENTS" => n))

matroid_list = results

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
lst = []
dims = []
for matroid in matroid_list
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
            push!(lst, matroid)
        end
    end
end


n_3 = open("$(n)_3_auto", "w")
n_3_L = open("$(n)_3_L_auto", "w")

for (i, matroid) in enumerate(lst)
    try
        R = realization_space(matroid, ground_ring = QQ)
        I = gens(defining_ideal(R))
        L = inequations(R)
        I_str = "I$(i) = ideal("
        L_str = "L$(i) = {"
        i += 1
        for f in I
            I_str *= string(f) * ", "
        end
        for g in L
            L_str *= string(g) * ", "
        end
        I_str = chop(I_str, tail=2)
        I_str *= ")\n\n"
        L_str = chop(L_str, tail=2)
        L_str *= "}\n\n"
        write(n_3, I_str)
        write(n_3_L, L_str)
        push!(dims, dim(R))
        print(hyperplanes(matroid))
    catch
        display(realization_space(matroid))
    end
end

close(n_3)
close(n_3_L)

println(length(lst))
println(dims)
