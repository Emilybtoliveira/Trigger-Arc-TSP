include("./instance.jl")

import .InstanceModule

if length(ARGS) < 1
    println("É necessário passar como argumento o caminho do arquivo de instância.")
    exit(1)
end

instance_filepath = ARGS[1]
instance_module::InstanceModule.Instance = InstanceModule.Instance(instance_filepath)

println("Number of nodes: ", instance_module.num_nodes)