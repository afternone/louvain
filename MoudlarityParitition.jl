using Graphs

type ModularityGroup{V}
    nodes::Set{V}
    inner_weight::Float64
    in_weight::Float64
    out_weight::Float64
end

type ModularityGraph{V}
    graph::AbstractGraph{V}
    weights::Vector{Float64}
    correct_self_loops::Bool
end

type ModularityPartition{V}
    modularity_graph::ModularityGraph{V}
    membership::Vector{Int}
    community::Dict{Int,ModularityGroup{V}}
    total_inner_weight::Float64 # total weight in all community
    total_inner_possible_edges::Int # total possible edges in all community
end

# update partition when membership vector changed
function update_partition!{V}(partition::ModularityPartition{V})
    # reset partition
    empty!(partition.community)
    partition.total_inner_weight = 0.0
    partition.total_inner_possible_edges = 0
    g = partition.modularity_graph.graph
    weights = partition.modularity_graph.weights

    for u in vertices(g)
        u_idx = vertex_index(u, g)
        u_comm = partition.membership[u_idx]
        if haskey(partition.community, u_comm)
            # Add this node to the community sets
            push!(partition.community[u_comm].nodes, u)
            # loop over all incident edges
            for e in out_edges(u, g)
                e_idx = edge_index(e, g)
                v = target(e, g)
                v_idx = vertex_index(v, g)
                v_comm = partition.membership[v_idx]
                w = weights[e_idx]
                # Add weight to the outgoing weight of community of u
                partition.community[u_comm].out_weight += w
                # Add weight to the incoming weight of community of v
                if haskey(partition.community, v_comm)
                    partition.community[v_comm].in_weight += w
                else
                    out_weight = 0.0
                    for e in out_edges(v, g)
                        e_idx = edge_index(e, g)
                        out_weight += weights[e_idx]
                    end
                    in_weight = 0.0
                    for e in in_edges(v, g)
                        e_idx = edge_index(e, g)
                        in_weight += weights[e_idx]
                    end
                    partition.community[v_comm] = ModularityGroup(Set(v), 0.0, in_weight, out_weight)
                end
                # if it is an edge within a community
                if u_comm == v_comm
                    if !is_directed(g)
                        w /= 2.0
                    end
                    partition.community[u_comm].inner_weight += w
                    partition.total_inner_weight += w
                end
            end
        else
            out_weight = 0.0
            for e in out_edges(u, g)
                e_idx = edge_index(e, g)
                out_weight += weights[e_idx]
            end
            in_weight = 0.0
            for e in in_edges(u, g)
                e_idx = edge_index(e, g)
                in_weight += weights[e_idx]
            end
            partition.community[u_comm] = ModularityGroup(Set(u), 0.0, in_weight, out_weight)
        end
    end

    for grp in values(partition.community)
        n_c = length(grp.nodes)
        if partition.modularity_graph.correct_self_loops
            possible_edges = round(Int, n_c*n_c/(2.0-Float64(is_directed(g))))
        else
            possible_edges = round(Int, n_c*(n_c-1)/(2.0-Float64(is_directed(g))))
        end
        partition.total_inner_possible_edges += possible_edges
    end
end

function modularity_partition{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T})
    membership = collect(1:num_vertices(g))
    community = Dict{Int,ModularityGroup{V}}()
    for u in vertices(g)
        u_idx = vertex_index(u, g)
        out_weight = 0.0
        for e in out_edges(u, g)
            e_idx = edge_index(e, g)
            out_weight += weights[e_idx]
        end
        in_weight = 0.0
        for e in in_edges(u, g)
            e_idx = edge_index(e, g)
            in_weight += weights[e_idx]
        end
        community[u_idx] = ModularityGroup(Set(u), 0.0, in_weight, out_weight)
    end
    mg = ModularityGraph(g, weights, false)
    ModularityPartition(mg, membership, community, 0.0, 0)
end

# Move a node to a new community and update the partition
function move_node!{V}(partition::ModularityPartition{V}, u::V, new_comm::Int)
    g = partition.modularity_graph.graph
    u_idx = vertex_index(u, g)
    old_comm = partition.membership[u_idx]



using GraphPlot
g = graphfamous("karate")
mp = modularity_partition(g, ones(78))
mp.community
mp.total_inner_possible_edges
m  = [1,1,1,1,1,1,1,1,2,1,1,1,1,1,2,2,1,1,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2]
mp.membership = m
update_partition!(mp)
mp.total_inner_possible_edges
