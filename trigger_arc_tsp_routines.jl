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
    println("The function TriggerArcTSP_ub_lp is not implemented yet.")
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

    (InArcs, OutArcs) = Build_In_Out_Arcs(T)
    model = Model(Gurobi.Optimizer)

    # x_a = 1 se arco a está na solução
    @variable(model, x_a[1:T.NArcs], Bin)

    # y_a = 1 se x_a = 1 e não há triggers ativos pro arco a
    @variable(model, y_a[1:T.NArcs], Bin)

    # y_r = 1 se trigger r está ativo
    @variable(model, y_r[1:T.NTriggers], Bin)

    # y_hat_r = 1 se trigger r não está ativo porque o arco target precede o trigger na solução
    @variable(model, y_hat_r[1:T.NTriggers], Bin)

    # u_i = ordem do nó i na solução
    @variable(model, 1 <= u[1:T.NNodes] <= T.NNodes, Int)


    # soma dos custos dos arcos sem trigger e com trigger
    @objective(model, Min,  sum(T.Trigger[i].cost * y_r[i] for i in 1:T.NTriggers)
                          + sum(T.Arc[i].cost * y_a[i] for i in 1:T.NArcs ))

    # R1: Quantidade de arcos na solução deve ser igual ao número de nós
    @constraint(model, R1, sum(x_a[j] for j in 1:T.NArcs) == T.NNodes ) 

    
    # R2 (MTZ): Se um arco (u,v) está na solução, entao u tem que vir antes de v na rota
    @constraint(model, R2, u[1] == 1) # O nó 1 deve ser o primeiro na rota
    @constraint(model, R2_[id in filter(i -> T.Arc[i].v != 1, 1:T.NArcs)], u[T.Arc[id].u] + 1 <= u[T.Arc[id].v] + (T.NNodes-1) * (1 - x_a[id])) #(T.NNodes-1) pra restrição ficar mais justa


    #R3: Conservação de fluxo    
    @constraint(model, R3_in[i in 1:T.NNodes], sum(x_a[ij] for ij in InArcs[i]) == 1) # Soma dos arcos que entram no nó i deve ser igual a 1
    @constraint(model, R3_out[i in 1:T.NNodes], sum(x_a[ij] for ij in OutArcs[i]) == 1) # Soma dos arcos que saem do nó i deve ser igual a 1
    
    #R4: ou y_a ou y_r devem ser 1 se o arco a está na solução
    for arc in 1:T.NArcs
        arc_triggers = filter(t -> T.Trigger[t].target_arc_id == arc, 1:T.NTriggers)
        @constraint(model, y_a[arc] + sum(y_r[trigger] for trigger in arc_triggers) == x_a[arc])    
    end

    
    for t in 1:T.NTriggers
        trigger = T.Trigger[t]
        id_trigger = trigger.trigger_arc_id
        id_target = trigger.target_arc_id
        
        u_trigger = T.Arc[id_trigger].u
        u_target = T.Arc[id_target].u
        
        #R5: Se o trigger está ativo, então o arco trigger deve estar na solução
        @constraint(model, y_r[t] <= x_a[id_trigger]) 
        
        #R6: Se um trigger está ativo, então o arco trigger deve vir antes do arco target na solução
        @constraint(model, u[u_trigger] + 1 <= u[u_target] + T.NNodes * (1 - y_r[t]))

        #R7: Se um trigger não está ativo, então o arco target deve vir antes do arco trigger na solução
        @constraint(model, u[u_target] + 1 <= u[u_trigger] + T.NNodes * (1 - y_hat_r[t]))
        
        #R8: Se o arco target está na solução, 
        @constraint(model, x_a[id_trigger] <= (1 - x_a[id_target]) + (1 - y_a[id_target]) + y_hat_r[t])

        for t2 in 1:T.NTriggers
            if t != t2
                trigger_2 = T.Trigger[t2]
                id_trigger_2 = trigger_2.trigger_arc_id                
                u_trigger_2 = T.Arc[id_trigger_2].u

                #R9: 
                # @constraint(model, u[u_trigger_2] - T.NNodes * y_hat_r[t2] <= u[u_trigger] + T.NNodes * (2 - y_r[t] - x_a[id_trigger_2]) - 1)
            end
        end
    end

    println("Model:", model)    

    optimize!(model)

    obj = objective_value(model) 
    println("Objective value: ", obj)
    if obj != nothing
        T.lb_ilp = objective_value(model)
        T.time_ilp = solve_time(model)
        T.ub_ilp = T.lb_ilp  # Assuming the model is feasible and optimal
        T.ub_ilp_arcs = [value(x_a[i]) for i in 1:T.NArcs]
        T.nn_ilp = sum(T.Arc[i].cost * value(x_a[i]) for i in 1:T.NArcs)
    else
        println("Model is not optimal or feasible.")
        T.time_ilp = 0.0
        T.lb_ilp = Inf
        T.ub_ilp = -Inf
        T.ub_ilp_arcs = zeros(T.NArcs)
        T.nn_ilp = 0.0
    end
    
    println("Time ILP: ", T.time_ilp)
    println("LB ILP: ", T.lb_ilp)
    println("UB ILP: ", T.ub_ilp)  
    println("UB ILP arcs: ", T.ub_ilp_arcs)
    println("NN ILP: ", T.nn_ilp)

end
# --------------------------------------------------------------
