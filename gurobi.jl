module GurobiModule
export GurobiModule

import Pkg
Pkg.add("Gurobi")

using JuMP, Gurobi


function solve_model()
    model = Model(Gurobi.Optimizer)
    println(model)
end

end