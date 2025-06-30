#
# Este deve ser o unico arquivo que deve ser atualizado pelo aluno
#
# Abaixo voce encontra as rotinas vazias das seguintes funcoes:
#                 TriggerArcTSP_lb_lp(T)
#                 TriggerArcTSP_lb_rlxlag(T)
#                 TriggerArcTSP_lb_colgen(T) - Opcional
#                 TriggerArcTSP_ub_lp(T)
#                 TriggerArcTSP_ub_rlxlag(T)
#                 TriggerArcTSP_ub_colgen(T) - Opcional
#                 TriggerArcTSP_ilp(T)
#
using JuMP
using Gurobi

# --------------------------------------------------------------
function TriggerArcTSP_lb_lp(T::TriggerArcTSP)
# This routine only changes the fields
# time_lb_lp
# lb_lp
    println("The function TriggerArcTSP_lb_lp is not implemented yet.")
end

# --------------------------------------------------------------
function TriggerArcTSP_ub_lp(T::TriggerArcTSP)
# This routine only changes the fields
#
# time_ub_lp
# ub_lp
# ub_lp_arcs
    println("The function  TriggerArcTSP_ub_lp is not implemented yet.")
end

# --------------------------------------------------------------
function TriggerArcTSP_lb_rlxlag(T::TriggerArcTSP)
# This routine only changes the fields
# time_lb_rlxlag
# lb_rlxlag
    println("The function TriggerArcTSP_lb_rlxlag is not implemented yet.")
end

# --------------------------------------------------------------
function TriggerArcTSP_ub_rlxlag(T::TriggerArcTSP)
# This routine only changes the fields
# time_ub_rlxlag
# ub_rlxlag
# ub_rlxlag_arcs
    println("The function TriggerArcTSP_ub_rlxlag is not implemented yet.")
end
# --------------------------------------------------------------
function TriggerArcTSP_lb_colgen(T::TriggerArcTSP)
# This routine only changes the fields
# time_lb_colgen
# lb_colgen
    println("The function TriggerArcTSP_lb_colgen is not implemented yet.")
end

# --------------------------------------------------------------
function TriggerArcTSP_ub_colgen(T::TriggerArcTSP)
# This routine only changes the fields
# time_ub_colgen
# ub_colgen
# ub_colgen_arcs
    println("The function TriggerArcTSP_ub_colgen is not implemented yet.")
end

# --------------------------------------------------------------
function TriggerArcTSP_ilp(T::TriggerArcTSP)
# This routine only changes the fields
# time_ilp
# lb_ilp
# ub_ilp
# ub_ilp_arcs
# nn_ilp
    model = Model(Gurobi.Optimizer)

    c = [ 1; 2; 5]
    A = [-1  1  3;
        1  3 -7]
    b = [-5; 10]

    @variable(model, x[1:3] >= 0)
    @objective(model, Max, sum( c[i]*x[i] for i in 1:3) )

    @constraint(model, constraint[j in 1:2], sum( A[j,i]*x[i] for i in 1:3 ) <= b[j] )
    @constraint(model, bound, x[1] <= 10)


    println("Model:", model)

    JuMP.optimize!(model)


    println("Optimal Solutions:")
    for i in 1:3
        println("x[$i] = ", JuMP.value(x[i]))
    end

        println("Dual Variables:")
    for j in 1:2
        println("dual[$j] = ", JuMP.shadow_price(constraint[j]))
    end
end
# --------------------------------------------------------------
