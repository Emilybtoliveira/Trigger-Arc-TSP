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