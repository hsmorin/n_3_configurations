#Computes all n_3 configurations as 
import Combinatorics as comb
import Polymake
using Oscar

function n_3(n::Int)
    L1 = [1, 2, 3]
    degArray = zeros(Int, n)
    degArray[1:3] .= (1,1,1)
    configs = [[degArray, [L1]]]

    function getNextConfigurations(config)

        degArray = config[1]
        newConfigs = []
        avaliablePoints = []

        for i in 1:n
            if degArray[i] != 3
                push!(avaliablePoints, i)
            end
        end
        
        m = length(avaliablePoints)
        L1 = avaliablePoints[1]
        #compute n choose 2 of avaliablePoints[2:n]
        
        for (i, p) in enumerate(avaliablePoints[2:(m - 1)])
            for (j,q) in enumerate(avaliablePoints[3:m])
                if [L1, p, q] ∉ config[2]
                    newLine = zeros(Int, 3)
                    newLine[1] += L1
                    newLine[2] += p
                    newLine[3] += q
                    newDegreeArray = copy(degArray)
                    newDegreeArray[L1] += 1
                    newDegreeArray[p] += 1
                    newDegreeArray[q] += 1
                    newConfig = copy(config)
                    newConfig[1] = newDegreeArray
                    newConfig[2] = vcat(newConfig[2], [newLine])
                    push!(newConfigs, newConfig) 
                end
            end
        end
        return(newConfigs)
    end

    function recursiveStep(configs) 

        if configs[1][1] >= fill(3, n)
            return [conf for conf in configs if conf[1] == fill(3,n)]
        end

        allNewConfigs = []
        for config in configs
            newConfigs = getNextConfigurations(config)
            append!(allNewConfigs, newConfigs)
        end
        allNewConfigs = [config for config in unique(allNewConfigs) if config[1] <= fill(3,n)]
        if length(allNewConfigs) >= 1000
            allNewConfigs = allNewConfigs[1:1000]
        end
        recursiveStep(allNewConfigs)
    end

    return recursiveStep(configs)
end
#r3n11 = find_one(db["Combinatorics.SelfProjectingMatroids"], Dict(["data.rank"=>"3", "data.length_groundset"=>"11", "data.dim_s"=>nothing]))
M1 = matroid_from_revlex_basis_encoding("0000************0**********0***************0***********0***********0*****0**************************************0**************************************0*************", 3, 11)

#Not expected dimension
M2 = matroid_from_revlex_basis_encoding("0000************0**********0***************0***********0***********0*****0**************************************0**************************************0***********0*", 3, 11)
#L5 = matroid_from_circuits([[1,2,4], [1,3,8], [1,7,9], [2,3,7], [2,5,9], [3,5,0], [4,5,6], [4,8,0], [6,7,8], [6,9,0]], [10, 0:9])
M3 = matroid_from_revlex_basis_encoding("0000************0**********0***************0***********0***********0*****0**************************************0**************************************0********0****", 3, 11)

#D = dim(realization_space(matroid, ground_ring = QQ)) 
H = [[1,2,4],[1,3,8],[1,7,9],[2,3,7],[2,5,9],[3,5,0],[4,5,6],[4, 8, 0], [6, 7, 8], [6, 9, 0], [3,4,9]];
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
end
println(H)
M = matroid_from_hyperplanes(H, 0:10);
H = hyperplanes(M)
println(length(H), " 3-planes")
println(hyperplanes(M))
println(M)
display(realization_space(M, ground_ring = QQ)) 
M = matroid_from_hyperplanes([[1,2,3],[1,6,0],[1,7,9],[2,4,0],[2,5,6],[3,4,9],[3,5,7],[4,5,8],[6,7,8],[8,9,0,]], 0:9)
H = [[1,2,4],[2,3,5],[1,3,6],[3,4,7],[1,5,7],[2,6,7],[4,6,5]];
M = matroid_from_hyperplanes(H,7)
