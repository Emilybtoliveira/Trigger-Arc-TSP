module InstanceModule
export InstanceModule


struct Instance 
    graph::Matrix{Int}
    num_nodes::Int
    num_edges::Int

    function Instance(filepath::String) 
        read_instance(filepath)
        new(zeros(Int, 0, 0), 0, 0)
    end
end

function read_instance(filepath::String)
    lines = readlines(filepath)
    println(lines)
    return lines
end

end