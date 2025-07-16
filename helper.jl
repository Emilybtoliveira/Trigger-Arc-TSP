function VerifyIfSolutionIsFeasible(T, objective, x)
    if count(s -> s == 1, x) != T.NNodes
        println("The number of arcs selected is different than the number of nodes.")
        return false
    end    
    
    route = GetRouteFromSolution(T, x)
    if isempty(route)
        return false
    end

    real_cost = GetRouteCost(T, route)	
    
    println("Solution is feasible.\nReal cost of the solution $real_cost")
    
    return real_cost == objective
end

function GetRouteFromSolution(T, x)
    visited = Set{Int}()
    route = Int[]
    nodes_arcs::Vector{IntType} = Vector{IntType}(undef,T.NNodes)

    for node in 1:T.NNodes
        nodes_arcs[node] = -1  
    end

    for arc in 1:T.NArcs
        if (x[arc] > 0.0001)
            if nodes_arcs[T.Arc[arc].u] != -1
                println("Solution sets two arcs leaving node $(T.Arc[arc].u)")
                return []
            end

            nodes_arcs[T.Arc[arc].u] = arc
        end
    end        

    current_node = 1
    while current_node ∉ visited
        push!(visited, current_node)
        arc = nodes_arcs[current_node]

        if arc == -1
            println("Node $current_node doesn't have outgoing arc.")
            return []
        end

        push!(route, arc)
        current_node = T.Arc[arc].v
    end

    if length(route) != T.NNodes || current_node != 1
        println("Not an Hamiltonian cycle.")
        return []
    end
    
    return route
end

function GetRouteCost(T, route)
    cost = 0.0
    for i in 1:length(route)
        arc = route[i]
        arc_cost = T.Arc[arc].cost
        triggers_of_arc_i = filter(t -> T.Trigger[t].target_arc_id == arc, 1:T.NTriggers)

        for j in 1:i-1  
            candidate_trigger = route[j]
            trigger_found = findfirst(t -> T.Trigger[t].trigger_arc_id == candidate_trigger, triggers_of_arc_i)
            if trigger_found !== nothing
                id_trigger = triggers_of_arc_i[trigger_found]
                arc_cost = T.Trigger[id_trigger].cost
                break
            end
        end
        cost += arc_cost
    end

    return cost
end


function BuildGreedyRouteFromFractionalSolution(T::TriggerArcTSP, x_frac::Vector{Float64})
    visited = Set{Int}()
    route = Int[]
    used_arcs = Set{Int}()
    
    current_node = 1
    push!(visited, current_node)

    while length(visited) < T.NNodes
        best_arc = -1
        best_score = -1.0

        for arc_id in 1:T.NArcs
            arc = T.Arc[arc_id]
            if arc.u == current_node && !(arc.v in visited) && !(arc_id in used_arcs)
                if x_frac[arc_id] > best_score 
                    best_score = x_frac[arc_id]
                    best_arc = arc_id
                end
            end
        end

        if best_arc == -1
            # println("No feasible arc leaving node $current_node")
            return []  
        end

        push!(route, best_arc)
        push!(visited, T.Arc[best_arc].v)
        push!(used_arcs, best_arc)

        current_node = T.Arc[best_arc].v
    end

    final_arc = findfirst(a -> T.Arc[a].u == current_node && T.Arc[a].v == 1, 1:T.NArcs)
    if final_arc == nothing
        # println("No arc to close cicle")
        return []
    end

    push!(route, final_arc)
    return route
end

function ValidateConstraints(T, x_a, y_a, y_r, y_hat_r, u)
    ε = 1e-4  # tolerância

    # 1. Verifica R1: número de arcos na solução
    if abs(sum(x_a) - T.NNodes) > ε
        println("R1 violada: sum(x_a) = $(sum(x_a)) ≠ $(T.NNodes)")
        
    end

    # 2. Verifica R2: ordem MTZ
    for arc in 1:T.NArcs
        u_i = T.Arc[arc].u
        v_i = T.Arc[arc].v
        if v_i == 1  # ignora MTZ com v=1
            continue
        end
        if x_a[arc] > 1 - ε && !(u[u_i] + 1 <= u[v_i] + ε)
            println("R2 violada para arco $arc: u[$u_i]+1 = $(u[u_i]+1) > u[$v_i] = $(u[v_i])")
            
        end
    end

    # 3. Verifica R3: conservação de fluxo
    in_deg = zeros(Float64, T.NNodes)
    out_deg = zeros(Float64, T.NNodes)
    for arc in 1:T.NArcs
        u_i = T.Arc[arc].u
        v_i = T.Arc[arc].v
        out_deg[u_i] += x_a[arc]
        in_deg[v_i]  += x_a[arc]
    end
    for i in 1:T.NNodes
        if abs(in_deg[i] - 1) > 0
            println("R3_in violada no nó $i: in_deg = $(in_deg[i])")
            
        end
        if abs(out_deg[i] - 1) > 0
            println("R3_out violada no nó $i: out_deg = $(out_deg[i])")
            
        end
    end

    # 4. Verifica R4: x_a = y_a + ∑y_r
    for arc in 1:T.NArcs
        related_triggers = filter(t -> T.Trigger[t].target_arc_id == arc, 1:T.NTriggers)
        
        y_sum = y_a[arc] + sum(y_r[t] for t in related_triggers; init = 0.0)
        if abs(x_a[arc] - y_sum) > 0
            println("R4 violada para arco $arc: x = $(x_a[arc]), y_sum = $y_sum")
            
        end
    end

    # 5. Verifica R5–R9: por trigger
    for t in 1:T.NTriggers
        trigger_id = T.Trigger[t].trigger_arc_id
        target_id  = T.Trigger[t].target_arc_id
        u_trigger = T.Arc[trigger_id].u
        u_target  = T.Arc[target_id].u
        
        # R5
        if y_r[t] > x_a[trigger_id]
            println("R5 violada: y_r[$t] > x_a[trigger] = $(x_a[trigger_id])")
            
        end

        # R6
        if !(u[u_trigger] + 1 <= u[u_target] + T.NNodes * (1 - y_r[t]))
            println("R6 violada: trigger $t, u_trigger = $(u[u_trigger]), u_target = $(u[u_target]), y_r = $(y_r[t])")
            
        end

        # R7
        if !(u[u_target] + 1 <= u[u_trigger] + T.NNodes * (1 - y_hat_r[t]))
            println("R7 violada: trigger $t, x_a $(x_a[target_id]), u_target = $(u[u_target]), u_trigger = $(u[u_trigger]), y_hat_r = $(y_hat_r[t])")
            
        end

        # R8
        lhs = x_a[trigger_id]
        rhs = (1 - x_a[target_id]) + (1 - y_a[target_id]) + y_hat_r[t]
        if lhs > rhs
            println("R8 violada: trigger $t, lhs = $lhs, rhs = $rhs")                    
        end

        # R9: múltiplos triggers para o mesmo target
        for t2 in (t+1):T.NTriggers
            if T.Trigger[t2].target_arc_id == target_id
                trigger2_id = T.Trigger[t2].trigger_arc_id
                u1 = u[T.Arc[trigger_id].u]
                u2 = u[T.Arc[trigger2_id].u]
                lhs = u2 - T.NNodes * y_hat_r[t2]
                rhs = u1 + T.NNodes * (2 - y_r[t] - x_a[trigger2_id]) - 1
                if lhs > rhs
                    println("R9 violada entre triggers $t e $t2")                            
                end
            end
        end
    end
end
