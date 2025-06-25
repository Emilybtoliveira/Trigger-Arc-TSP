module InstanceModule
export InstanceModule


struct Instance 
    graph::Matrix{Float64}
    nodes::Int
    edges::Int
    relationships::Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int, Float64}}}

    function Instance(filepath::String) 
        n, m, r, graph, relationships = read_instance(filepath)
        new(graph, n, m, relationships)
    end
end

function read_instance(filepath::String)
    lines = readlines(filepath)
    # println(lines)
    
    n, m, r = split(lines[1], ' ')
    println("n: $n, m: $m, r: $r")

    n = parse(Int, n)
    m = parse(Int, m)
    r = parse(Int, r)
    matrix = zeros(Float64, n, n)

    # Cria o grafo
    for i = 2:m+1
        index, node1, node2, cost = split(lines[i], ' ')
        # println(lines[i])

        node1 = parse(Int, node1) + 1
        node2 = parse(Int, node2) + 1
        matrix[node1, node2] = parse(Float64, cost)
    end

    # println("Graph: ", matrix)

    relationships = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int, Float64}}}()

    # Cria a matriz de relacionamentos
    for i = (m+2):length(lines)
        index1, index2, n1, n2, index3, n3, n4, cost = split(lines[i], ' ')
        # println(lines[i])

        n1 = parse(Int, n1) + 1
        n2 = parse(Int, n2) + 1
        n3 = parse(Int, n3) + 1
        n4 = parse(Int, n4) + 1

        key = (n1, n2)
        value = (n3, n4, parse(Float64, cost)) #provavelmente vai ter que mudar a key pra ser o target e a lista pra ser os triggers
        if haskey(relationships, key)
            push!(relationships[key], value)
        else
            relationships[key] = [value]
        end
    end

    # println("Relationships: ", relationships)

    return n, m, r, matrix, relationships
end

end