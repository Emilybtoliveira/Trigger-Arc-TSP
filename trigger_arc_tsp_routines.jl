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
using DataStructures
using Plots

function Build_TATSP_Base_Model(model::JuMP.Model, T::TriggerArcTSP, relaxed_vars::Bool=false, active_constraints::Dict{String,Bool}=Dict{String,Bool}(), default_objective::Bool=true)  
    if !relaxed_vars 
        # Variáveis inteiras
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
    else 
        #Variáveis continuas
        @variable(model, 0 <= x_a[1:T.NArcs] <= 1)
        @variable(model, 0 <= y_a[1:T.NArcs] <= 1)
        @variable(model, 0 <= y_r[1:T.NTriggers] <= 1)
        @variable(model, 0 <= y_hat_r[1:T.NTriggers] <= 1)
        @variable(model, 1 <= u[1:T.NNodes] <= T.NNodes, Int)
    end
    

    # soma dos custos dos arcos sem trigger e com trigger
    if default_objective
        @objective(model, Min,  sum(T.Trigger[i].cost * y_r[i] for i in 1:T.NTriggers)
                            + sum(T.Arc[i].cost * y_a[i] for i in 1:T.NArcs ))
    end

    # R1: Quantidade de arcos na solução deve ser igual ao número de nós
    if get(active_constraints, "R1", true)
        @constraint(model, R1, sum(x_a[j] for j in 1:T.NArcs) == T.NNodes ) 
    end

    
    # R2 (MTZ): Se um arco (u,v) está na solução, entao u tem que vir antes de v na rota
    if get(active_constraints, "R2", true)
        @constraint(model, R2, u[1] == 1) # O nó 1 deve ser o primeiro na rota
        @constraint(model, R2_[id in filter(i -> T.Arc[i].v != 1, 1:T.NArcs)], u[T.Arc[id].u] + 1 <= u[T.Arc[id].v] + (T.NNodes-1) * (1 - x_a[id])) #(T.NNodes-1) pra restrição ficar mais justa
    end
    
    #R3: Conservação de fluxo    
    (InArcs, OutArcs) = Build_In_Out_Arcs(T)
    if get(active_constraints, "R3", true)
        @constraint(model, R3_in[i in 1:T.NNodes], sum(x_a[ij] for ij in InArcs[i]) == 1) # Soma dos arcos que entram no nó i deve ser igual a 1
        @constraint(model, R3_out[i in 1:T.NNodes], sum(x_a[ij] for ij in OutArcs[i]) == 1) # Soma dos arcos que saem do nó i deve ser igual a 1
    end

    #R4: ou y_a ou y_r devem ser 1 se o arco a está na solução
    if get(active_constraints, "R4", true)        
        for arc in 1:T.NArcs
            arc_triggers = filter(t -> T.Trigger[t].target_arc_id == arc, 1:T.NTriggers)
            @constraint(model, y_a[arc] + sum(y_r[trigger] for trigger in arc_triggers) == x_a[arc], base_name="R4_$(arc)")    
        end
    end
    
    for t in 1:T.NTriggers
        trigger = T.Trigger[t]
        id_trigger = trigger.trigger_arc_id
        id_target = trigger.target_arc_id
        
        u_trigger = T.Arc[id_trigger].u
        u_target = T.Arc[id_target].u
        
        #R5: Se o trigger está ativo, então o arco trigger deve estar na solução
        if get(active_constraints, "R5", true)
            @constraint(model, y_r[t] <= x_a[id_trigger], base_name="R5_$(t)") 
        end
        
        #R6: Se um trigger está ativo, então o arco trigger deve vir antes do arco target na solução
        if get(active_constraints, "R6", true) 
            @constraint(model, u[u_trigger] + 1 <= u[u_target] + T.NNodes * (1 - y_r[t]), base_name="R6_$(t)")  
        end
        
        #R7: Se um trigger não está ativo, então o arco target deve vir antes do arco trigger na solução
        if get(active_constraints, "R7", true)
            @constraint(model, u[u_target] + 1 <= u[u_trigger] + T.NNodes * (1 - y_hat_r[t]),  base_name="R7_$(t)")
        end

        #R8
        if get(active_constraints, "R8", true)
            @constraint(model, x_a[id_trigger] <= (1 - x_a[id_target]) + (1 - y_a[id_target]) + y_hat_r[t], base_name="R8_$(t)")
        end

        if get(active_constraints, "R9", true)            
            for t2 in 1:T.NTriggers
                if t != t2 && T.Trigger[t2].target_arc_id == id_target && t < t2
                    trigger_2 = T.Trigger[t2]
                    id_trigger_2 = trigger_2.trigger_arc_id                
                    u_trigger_2 = T.Arc[id_trigger_2].u

                    #R9: Se dois triggers estão para o mesmo target, o último que aparece na rota é que deve estar ativo
                    @constraint(model, u[u_trigger_2] - (T.NNodes * y_hat_r[t2]) <= u[u_trigger] + (T.NNodes * (2 - y_r[t] - x_a[id_trigger_2]) - 1), base_name="R9_$(t)_$(t2)")
                end
            end
        end
    end  
end

function Build_lagrangian_objective(model::JuMP.Model, T::TriggerArcTSP, lambda::Vector{FloatType}, gama::Dict{Tuple{IntType,IntType},FloatType}, delta::Vector{FloatType}, mu::Vector{FloatType})
    x_a = model[:x_a]
    y_a = model[:y_a]
    y_r = model[:y_r]
    y_hat_r = model[:y_hat_r]
    u = model[:u]

    relax_objective = AffExpr()

    for i in 1:T.NTriggers
        add_to_expression!(relax_objective, T.Trigger[i].cost, y_r[i])
    end
    for i in 1:T.NArcs
        add_to_expression!(relax_objective, T.Arc[i].cost, y_a[i])
    end

    for t in 1:T.NTriggers
        trigger  = T.Trigger[t]
        trigger_id = trigger.trigger_arc_id
        target_id  = trigger.target_arc_id

        expr_R6 = u[T.Arc[trigger_id].u] + 1 - u[T.Arc[target_id].u] - (T.NNodes * (1 - y_r[t]))
        add_to_expression!(relax_objective, delta[t], expr_R6)

        expr_R7 = u[T.Arc[trigger_id].u] + 1 - u[T.Arc[target_id].u] - (T.NNodes * (1 - y_hat_r[t]))
        add_to_expression!(relax_objective, mu[t], expr_R7)

        expr_R8 = x_a[trigger_id] + x_a[target_id] + y_a[target_id] - y_hat_r[t] - 2
        add_to_expression!(relax_objective, lambda[t], expr_R8)
    end

    for t1 in 1:T.NTriggers
        for t2 in t1+1:T.NTriggers
            if T.Trigger[t1].target_arc_id == T.Trigger[t2].target_arc_id
                trigger_1_id = T.Trigger[t1].trigger_arc_id
                trigger_2_id = T.Trigger[t2].trigger_arc_id

                expr_R9 = u[T.Arc[trigger_2_id].u] - u[T.Arc[trigger_1_id].u] -  (T.NNodes * y_hat_r[t2]) - (T.NNodes * (2 - y_r[t1] - x_a[trigger_2_id])) + 1

                add_to_expression!(relax_objective, gama[(t1,t2)], expr_R9)
            end
        end
    end

    return relax_objective
end

function Get_Active_Constraints_For_lb_rlxlag()
    active_constraints = Dict(
        "R1" => true,
        "R2" => true,
        "R3" => true,
        "R4" => true,
        "R5" => true,
        "R6" => false,  # Relaxa a R6
        "R7" => false,  # Relaxa a R7
        "R8" => false,  # Relaxa a R8
        "R9" => false   # Relaxa a R9
    )
    return active_constraints
end 

function Apply_2_OPT_Heuristic(T::TriggerArcTSP, lb_sol::Vector{FloatType})
    current_route = BuildGreedyRouteFromFractionalSolution(T, lb_sol)  
    
    if length(current_route) == 0
        return Inf, nothing, nothing, nothing, nothing, nothing
    end
    
    node_seq = [T.Arc[a].u for a in current_route] 
    push!(node_seq, node_seq[1]) 
    
    best_route = copy(node_seq)
    best_cost = GetRouteCost(T, current_route)
    max_iterations = 100
    improved = true
    loops = 1
    
    while improved && loops <= max_iterations
        improved = false
        for i in 2:(T.NNodes - 1)
            for j in (i + 1):T.NNodes
                new_route = vcat(best_route[1:i-1], reverse(best_route[i:j]), best_route[j+1:end])

                if length(unique(new_route[1:end-1])) != T.NNodes
                    continue
                end

                valid = true
                arc_list = IntType[]
                for k in 1:T.NNodes
                    u = new_route[k]
                    v = new_route[k + 1]
                    arc = findfirst(a -> T.Arc[a].u == u && T.Arc[a].v == v, 1:T.NArcs)
                    if arc === nothing
                        valid = false
                        break
                    end
                    push!(arc_list, arc)
                end

                if !valid
                    continue
                end

                cost = GetRouteCost(T, arc_list)
                if cost < best_cost
                    best_cost = cost
                    best_route = new_route
                    improved = true
                    break
                end
            end
        end
        loops += 1
    end

    arc_route = [findfirst(a -> T.Arc[a].u == best_route[i] && T.Arc[a].v == best_route[i+1], 1:T.NArcs) for i in 1:T.NNodes]
    
    x_a = zeros(IntType, T.NArcs)
    for a in arc_route
        x_a[a] = 1
    end

    u = zeros(IntType, T.NNodes)
    for i in 1:T.NNodes
        u[best_route[i]] = i
    end

    y_a = zeros(IntType, T.NArcs)
    y_r = zeros(IntType, T.NTriggers)
    y_hat_r = zeros(IntType, T.NTriggers)

    for t in 1:T.NTriggers
        id_trigger = T.Trigger[t].trigger_arc_id
        id_target  = T.Trigger[t].target_arc_id

        if u[T.Arc[id_target].u] < u[T.Arc[id_trigger].u]
            y_hat_r[t] = 1
        end
    end

    for arc_target in 1:T.NArcs
        if x_a[arc_target] != 1
            continue  # arco target não está na solução
        end

        triggers_for_target = filter(t -> T.Trigger[t].target_arc_id == arc_target, 1:T.NTriggers)

        best_t = nothing
        best_pos = -1
        for t in triggers_for_target
            id_trigger = T.Trigger[t].trigger_arc_id

            if x_a[id_trigger] == 1 && y_hat_r[t] == 0
                u_pos = u[T.Arc[id_trigger].u]
                if u_pos > best_pos
                    best_pos = u_pos
                    best_t = t
                end
            end
        end

        if best_t !== nothing
            y_r[best_t] = 1
        end
    end

    for a in 1:T.NArcs
        if x_a[a] == 1 && all(t -> !(T.Trigger[t].target_arc_id == a && y_r[t] == 1), 1:T.NTriggers)
            y_a[a] = 1
        end
    end

    return best_cost, x_a, y_a, y_r, y_hat_r, u
end

function Plot_Bounds(T::TriggerArcTSP, lb_values::Vector{FloatType}, ub_values::Vector{FloatType}, iter_values::Vector{IntType})
    plot(iter_values, lb_values;
        label = "Lower Bound",
        xlabel = "Iteração",
        ylabel = "Bound",
        # title = "Evolução dos Limitantes na Heurística Lagrangiana",
        lw = 2,
        color = :blue,
        marker = :circle,        
        markersize = 4,          
        legend = :bottomright    
    )

    plot!(iter_values, ub_values;
        label = "Upper Bound",
        lw = 2,
        color = :red,
        marker = :square,          
        markersize = 4
    )

    instance_name = splitext(basename(T.inputfilename))[1]
    savefig("./figs/rlxlag_$(instance_name).png")
end

# --------------------------------------------------------------
function TriggerArcTSP_lb_lp(T::TriggerArcTSP)
    println("===== Running TriggerArcTSP_lb_lp...")

    model = Model(Gurobi.Optimizer)    
    
    Build_TATSP_Base_Model(model, T, true)
    if verbose_mode
        println("Model:", model)    
    end
    
    println("Built the model with $(num_variables(model)) variables and $(num_constraints(model, count_variable_in_set_constraints = false)) constraints.")
    
    set_time_limit_sec(model, T.maxtime_lb_lp)
    # set_objective_sense(model, FEASIBILITY_SENSE)
    optimize!(model)
    
    # println(primal_status(model))    
    
    if has_values(model)    
        x_a = model[:x_a]
        u = model[:u]
        
        T.lb_lp = objective_bound(model) 
        T.time_lb_lp = solve_time(model)
        T.ub_lp = objective_value(model)
        T.ub_lp_arcs = [value(x_a[i]) for i in 1:T.NArcs]
        gap = relative_gap(model)
        u = [value(u[i]) for i in 1:T.NNodes]

        # is_feasible = VerifyIfSolutionIsFeasible(T, T.ub_lp, T.ub_lp_arcs)
        # println("Result of feasibility evaluation: $is_feasible")
    else
        println("No solution was found within time limit.")
        T.time_lb_lp = solve_time(model)
        gap = Inf
    end    
 
    
    println("Optimality Gap: ", gap)
    println("Time LB LP: ", T.time_lb_lp)
    println("LB LP: ", T.lb_lp)
    println("UB LP: ", T.ub_lp)  
    println("UB LP arcs: ", T.ub_lp_arcs)

    println("Finished TriggerArcTSP_lb_lp ======== \n\n")
end

function TriggerArcTSP_ub_lp(T::TriggerArcTSP) 
    println("===== Running TriggerArcTSP_ub_lp...")

    start_time = time()    
    model = Model(Gurobi.Optimizer)
    Build_TATSP_Base_Model(model, T)
    relax_integrality(model)
    # println("Built the model with $(num_variables(model)) variables and $(num_constraints(model, count_variable_in_set_constraints = false)) constraints.")

    # println(model)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_time_limit_sec(model, 10) #10s para resolver cada modelo
    optimize!(model)
    
    if !has_values(model)  
        println("Could not solve relaxation within time limit.")
    end
    
    println("Relaxation solution: $(objective_value(model))")
    x_lp = value.(model[:x_a])

    TriggerArcTSP_ub_rlxlag(T)
    x_rlx_lag = T.ub_rlxlag_arcs 

    fixed_vars = 0
    for a in 1:T.NArcs
        if x_lp[a] > 0.9999 && x_rlx_lag[a] > 0.9999
            # println("Fixing x_a[$a] = 1")
            fix(model[:x_a][a], 1.0; force=true)
            fixed_vars += 1
        elseif 1 - x_lp[a] > 0.9999 && 1 - x_rlx_lag[a] > 0.9999
            # println("Fixing x_a[$a] = 0")
            fix(model[:x_a][a], 0.0; force=true)
            fixed_vars += 1
        else
            set_binary(model[:x_a][a])
        end
    end

    T.ub_lp_fixed_vars = fixed_vars
    for v in model[:y_a] set_binary(v) end
    for v in model[:y_r] set_binary(v) end
    for v in model[:y_hat_r] set_binary(v) end
    for v in model[:u] set_integer(v) end

    set_time_limit_sec(model, 30) 
    # set_optimizer_attribute(model, "OutputFlag", 1)
    optimize!(model)
    # println(primal_status(model))

    if has_values(model)
        T.ub_lp = objective_value(model)
        T.ub_lp_arcs = value.(model[:x_a])
        VerifyIfSolutionIsFeasible(T, T.ub_lp, T.ub_lp_arcs)    
        
        new_cost, new_route, _,_,_,_ =  Apply_2_OPT_Heuristic(T, T.ub_lp_arcs)
        # println("2-opt cost: $new_cost")
        
        T.ub_lp = new_cost
        T.ub_lp_arcs = new_route        
        T.time_ub_lp = time() - start_time
    else
        println("No feasible solution found in RINS search.")
    end        

    println("Time UB LP: ", T.time_ub_lp)
    println("UB LP: ", T.ub_lp)  
    println("UB LP arcs: ", T.ub_lp_arcs)
    println("Fixed variables: ", T.ub_lp_fixed_vars)

    println("Finished TriggerArcTSP_ub_lp ======== \n\n")
end
# --------------------------------------------------------------

function TriggerArcTSP_lb_rlxlag(T::TriggerArcTSP; 
                                 model = nothing,
                                 lambda::Vector{FloatType} = zeros(T.NTriggers), 
                                 gama::Dict{Tuple{IntType,IntType}, FloatType} = Dict{Tuple{IntType,IntType}, FloatType}(), 
                                 delta::Vector{FloatType} = zeros(T.NTriggers), 
                                 mu::Vector{FloatType} = zeros(T.NTriggers),
                                 print_running::Bool = false)
    
    if print_running println("===== Running TriggerArcTSP_lb_rlxlag...") end

    if model == nothing
        model = Model(Gurobi.Optimizer)    
        active_constraints = Get_Active_Constraints_For_lb_rlxlag()
        Build_TATSP_Base_Model(model, T, false, active_constraints, false)
    end 
    
    # Init gama with existing triggers if empty
    if isempty(gama)
        for t1 in 1:T.NTriggers, t2 in t1+1:T.NTriggers
            if T.Trigger[t1].target_arc_id == T.Trigger[t2].target_arc_id
                gama[(t1,t2)] = 0.0           
            end
        end
    end

    set_optimizer_attribute(model, "OutputFlag", 0)
    objective_function = Build_lagrangian_objective(model, T, lambda, gama, delta, mu)
    set_objective(model, MIN_SENSE, objective_function)

    if verbose_mode
        println("Model:", model)    
    end
    
    println("Built the model with $(num_variables(model)) variables and $(num_constraints(model, count_variable_in_set_constraints = false)) constraints.")
    set_time_limit_sec(model, T.maxtime_lb_rlxlag)

    optimize!(model)
    
    # println(primal_status(model))  

    if has_values(model)    
        T.lb_rlxlag = objective_value(model)
        T.time_lb_rlxlag = solve_time(model)
    else
        println("No solution was found within time limit.")
        T.lb_rlxlag = -Inf
        T.time_lb_rlxlag = solve_time(model)
    end 

    println("Time UB rlxlab: $(T.time_lb_rlxlag )")
    println("LB rlxlag: $(T.lb_rlxlag)")
    # VerifyIfSolutionIsFeasible(T, T.lb_rlxlag, [value(model[:x_a][i]) for i in 1:T.NArcs])

    if print_running println("Finished TriggerArcTSP_lb_rlxlag! ======== \n\n") end
    return T.lb_rlxlag
end

function TriggerArcTSP_ub_rlxlag(T::TriggerArcTSP, print_running::Bool = false)
    if print_running println("===== Running TriggerArcTSP_ub_rlxlag...") end

    step_size = 0.001
    theta = step_size
    k = 1

    lb_values = FloatType[]
    ub_values = FloatType[]
    iter_values = IntType[]
    
    #Initialize data structures
    UB = Inf 
    LB = -Inf 
    arcs_best_UB = zeros(IntType, T.NArcs)
    arcs_best_LB = zeros(IntType, T.NArcs)
    best_UB = Inf
    best_LB = -Inf
    arcs_UB = zeros(IntType, T.NArcs)
    arcs_LB = zeros(IntType, T.NArcs)
    
    #Initialize Lagrange multipliers
    delta = zeros(T.NTriggers)
    mu = zeros(T.NTriggers)
    lambda = zeros(T.NTriggers)
    gama = Dict{Tuple{IntType,IntType}, FloatType}()  
    for t1 in 1:T.NTriggers, t2 in t1+1:T.NTriggers
        if T.Trigger[t1].target_arc_id == T.Trigger[t2].target_arc_id
            gama[(t1,t2)] = 0.0           
        end
    end

    model = Model(Gurobi.Optimizer)    
    active_constraints = Get_Active_Constraints_For_lb_rlxlag()
    Build_TATSP_Base_Model(model, T, false, active_constraints, false)
    
    start_time = time()

    while (time() - start_time) < T.maxtime_ub_rlxlag
        LB = TriggerArcTSP_lb_rlxlag(T; model=model, lambda=lambda, gama=gama, delta=delta, mu=mu)

        x_a = model[:x_a]
        y_a = model[:y_a]
        y_r = model[:y_r]
        y_hat_r = model[:y_hat_r]
        u = model[:u]

        arcs_LB = [value(x_a[i]) for i in 1:T.NArcs]
        
        if LB > best_LB && LB <= best_UB
            best_LB = LB
            arcs_best_LB = arcs_LB
        end  
        
        UB, arcs_UB, _, _, _, _ = Apply_2_OPT_Heuristic(T, arcs_LB)
        if UB < best_UB
            best_UB = UB
            arcs_best_UB = arcs_UB
        end

        if verbose_mode
            println("Iteration: $k, LB: $LB, Best LB: $best_LB, UB:$UB, Best UB: $best_UB, Theta: $theta")
        end

        push!(lb_values, LB)
        push!(ub_values, UB)
        push!(iter_values, k)

        violations_R6 = zeros(T.NTriggers)
        violations_R7 = zeros(T.NTriggers)
        for t in 1:T.NTriggers
            trigger_id = T.Trigger[t].trigger_arc_id
            target_id  = T.Trigger[t].target_arc_id

            violations_R6[t] = value(u[T.Arc[trigger_id].u]) + 1 - value(u[T.Arc[target_id].u]) - T.NNodes * (1 - value(y_r[t]))
            violations_R7[t] = value(u[T.Arc[target_id].u]) + 1 - value(u[T.Arc[trigger_id].u]) - T.NNodes * (1 - value(y_hat_r[t]))
        end

        violations_R8 = similar(lambda)
        for t in 1:T.NTriggers
            trigger = T.Trigger[t].trigger_arc_id
            target = T.Trigger[t].target_arc_id

            violations_R8[t] = value(x_a[trigger]) + value(x_a[target]) + value(y_a[target]) - value(y_hat_r[t]) - 2 
        end

        violations_R9 = Dict{Tuple{IntType,IntType}, FloatType}()
        for (t1,t2) in keys(gama)
            trigger1 = T.Trigger[t1].trigger_arc_id
            trigger2 = T.Trigger[t2].trigger_arc_id

            violations_R9[(t1,t2)] = value(u[T.Arc[trigger2].u]) - value(u[T.Arc[trigger1].u]) - (T.NNodes * value(y_hat_r[t2])) - (T.NNodes * (2 - value(y_r[t1]) - value(x_a[trigger2]))) + 1
        end
        
        # Updates the Lagrange multipliers
        delta .= max.(0, delta .+ theta .* violations_R6)
        mu .= max.(0, mu .+ theta .* violations_R7)
        lambda .= max.(0, lambda .+ theta .* violations_R8)
        for p in keys(gama)
            gama[p] = max(0, gama[p] + theta * violations_R9[p])
        end

        # ========= Este bloco está comentado, porque acabou não sendo usado nos experimentos do relatório ========
        # Updates the step size using Polyak's rule.       
        # violations = [violations_R6...;violations_R7...;violations_R8...;values(violations_R9)...]
        # norm_sq = sum(v^2 for v in violations)
        
        # if isfinite(LB) && isfinite(best_UB) && norm_sq > 1e-8
        #     theta = (step_size * (best_UB - LB)) / norm_sq
        # end
        # =========


        k += 1
        # println("Elapsed time $(time() - start_time)")
    end

    T.time_ub_rlxlag =  time() - start_time
    T.ub_rlxlag = best_UB
    T.ub_rlxlag_arcs = arcs_best_UB

    println("Time UB rlxlab: ", T.time_ub_rlxlag)
    println("UB rlxlab: ", T.ub_rlxlag)
    println("LB rlxlab: ", best_LB)

    # Plot_Bounds(T, lb_values, ub_values, iter_values)

    if print_running println("Finished TriggerArcTSP_ub_rlxlag! ======== \n\n") end
end
# --------------------------------------------------------------

function TriggerArcTSP_lb_colgen(T::TriggerArcTSP)
    println("===== Running TriggerArcTSP_lb_colgen...")    
    println("The function TriggerArcTSP_lb_colgen is not implemented yet.")
    println("Finished TriggerArcTSP_lb_colgen ======== \n\n")
end

function TriggerArcTSP_ub_colgen(T::TriggerArcTSP)
    println("===== Running TriggerArcTSP_ub_colgen...")    
    println("The function TriggerArcTSP_ub_colgen is not implemented yet.")
    println("Finished TriggerArcTSP_ub_colgen ======== \n\n")
end
# --------------------------------------------------------------

function TriggerArcTSP_ilp(T::TriggerArcTSP)
    println("===== Running TriggerArcTSP_ilp...")

    active_constraints = Dict(
        "R1" => true,
        "R2" => true,
        "R3" => true, # lazy constraint
        "R4" => true, # lazy constraint
        "R5" => true,
        "R6" => true,
        "R7" => true,
        "R8" => false, # lazy constraint
        "R9" => false  # lazy constraint
    )
        
    model = Model(Gurobi.Optimizer)    
    Build_TATSP_Base_Model(model, T, false, active_constraints)
    if verbose_mode
        println("Model:", model)    
    end
    println("Built the model with $(num_variables(model)) variables and $(num_constraints(model, count_variable_in_set_constraints = false)) constraints.")


    function Cutting_Planes_Callback(cb_data)
        status = callback_node_status(cb_data, model)
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            # println("not integer")
            return  # Only run at integer solutions
        end    

        # println("Solution is: ", callback_value.(cb_data, model[:x_a]))

        x_a = callback_value.(cb_data, model[:x_a])
        y_a = callback_value.(cb_data, model[:y_a])
        y_r = callback_value.(cb_data, model[:y_r])
        y_hat_r = callback_value.(cb_data, model[:y_hat_r])
        u = callback_value.(cb_data, model[:u])
        epsilon = 1e-5
        cuts_added = 0
        
        if (!active_constraints["R4"])
            for arc in 1:T.NArcs
                if x_a[arc] >= (1 - epsilon) #arco selecionado
                    arc_triggers = filter(r -> T.Trigger[r].target_arc_id == arc, 1:T.NTriggers)
                    y_sum = y_a[arc] + sum(y_r[r] for r in arc_triggers; init = 0.0)
                    if y_sum < 1 - epsilon
                        #R4: ou y_a ou sum(y_r) devem ser 1 se o arco a está na solução
                        R4 = @build_constraint(model[:y_a][arc] + sum(model[:y_r][r] for r in arc_triggers) >= model[:x_a][arc])
                        MOI.submit(model, MOI.LazyConstraint(cb_data), R4)
                        cuts_added += 1
                    end
                end
            end
        end

        if (!active_constraints["R8"] && !active_constraints["R9"])
            for t1 in 1:T.NTriggers
                #R8
                trigger1_id = T.Trigger[t1].trigger_arc_id
                target_id  = T.Trigger[t1].target_arc_id

                rhs = (1 - x_a[target_id]) + (1 - y_a[target_id]) + y_hat_r[t1]

                if x_a[trigger1_id] > rhs + epsilon # violação da restrição
                    R8 = @build_constraint(model[:x_a][trigger1_id] <= (1 - model[:x_a][target_id]) + (1 - model[:y_a][target_id]) + model[:y_hat_r][t1])
                    MOI.submit(model, MOI.LazyConstraint(cb_data), R8)
                    cuts_added += 1
                end

                for t2 in (t1+1):T.NTriggers
                    #R9
                    if T.Trigger[t1].target_arc_id == T.Trigger[t2].target_arc_id
                        u1 = u[T.Arc[T.Trigger[t1].trigger_arc_id].u]
                        u2 = u[T.Arc[T.Trigger[t2].trigger_arc_id].u]

                        if u2 - (T.NNodes * y_hat_r[t2]) > u1 + (T.NNodes * (2 - y_r[t1] - x_a[T.Trigger[t2].trigger_arc_id])) - 1 + epsilon  # violação da restrição
                            R9 = @build_constraint(model[:u][T.Arc[T.Trigger[t2].trigger_arc_id].u] - (T.NNodes * model[:y_hat_r][t2]) <= 
                                                    model[:u][T.Arc[T.Trigger[t1].trigger_arc_id].u] + T.NNodes * (2 - model[:y_r][t1] - model[:x_a][T.Trigger[t2].trigger_arc_id]) - 1)
                            MOI.submit(model, MOI.LazyConstraint(cb_data), R9)
                            cuts_added += 1
                        end
                    end
                end
            end
        end

        if (!active_constraints["R3"])
        
            (InArcs, OutArcs) = Build_In_Out_Arcs(T)

            for i in 1:T.NNodes
                in_arcs_i = InArcs[i]  
                out_arcs_i = OutArcs[i]

                sum_in = sum(x_a[a] for a in in_arcs_i)
                sum_out = sum(x_a[a] for a in out_arcs_i)

                if abs(sum_in - 1) > epsilon
                    R3_1 = @build_constraint((sum(model[:x_a][a] for a in in_arcs_i)) == 1)
                    MOI.submit(model, MOI.LazyConstraint(cb_data), R3_1)
                    cuts_added += 1
                end

                if abs(sum_out - 1) > epsilon
                    R3_2 = @build_constraint((sum(model[:x_a][a] for a in out_arcs_i)) == 1)
                    MOI.submit(model, MOI.LazyConstraint(cb_data), R3_2)
                    cuts_added += 1
                end
            end
        end

        println("Added $cuts_added cuts.")
        return
    end  
    
    function Primal_Heuristic_Callback(cb_data)
        x_a_current = callback_value.(cb_data, model[:x_a])
        # println(x_a_current)
        
        cost, x_a_new, y_a, y_r, y_hat_r, u = Apply_2_OPT_Heuristic(T, x_a_current)

        if cost < Inf
            # println("Heuristic found solution of cost $cost")
            ValidateConstraints(T, x_a_new, y_a, y_r, y_hat_r, u)
            
            status = MOI.submit(model, MOI.HeuristicSolution(cb_data),
                [model[:x_a]; model[:y_a]; model[:y_r]; model[:y_hat_r]; model[:u]],
                [x_a_new; y_a; y_r; y_hat_r; u]
            )

            # println("Status of the solution: $status")
        end

        return
    end
        
    
    set_attribute(model, MOI.LazyConstraintCallback(), Cutting_Planes_Callback)    
    set_attribute(model, MOI.HeuristicCallback(), Primal_Heuristic_Callback)

    set_time_limit_sec(model, T.maxtime_ilp)
    
    optimize!(model)
    
    # println(primal_status(model))    
    
    if has_values(model)    
        x_a = model[:x_a]
        u = model[:u]
        
        T.lb_ilp = objective_bound(model) 
        T.time_ilp = solve_time(model)
        T.ub_ilp = objective_value(model)
        T.ub_ilp_arcs = [value(x_a[i]) for i in 1:T.NArcs]
        T.nn_ilp =  node_count(model)
        gap = relative_gap(model)
        u = [value(u[i]) for i in 1:T.NNodes]

        # is_feasible = VerifyIfSolutionIsFeasible(T, T.ub_ilp, T.ub_ilp_arcs)
        # println("Result of feasibility evaluation: $is_feasible")
    else
        println("No solution was found within time limit.")
        T.time_ilp = solve_time(model)
        gap = Inf
    end    
 
    
    println("Optimality Gap: ", gap)
    println("Time ILP: ", T.time_ilp)
    println("LB ILP: ", T.lb_ilp)
    println("UB ILP: ", T.ub_ilp)  
    println("UB ILP arcs: ", T.ub_ilp_arcs)
    println("Number of nodes: ", T.nn_ilp)

    println("Finished TriggerArcTSP_ilp ======== \n\n")
end
# --------------------------------------------------------------
