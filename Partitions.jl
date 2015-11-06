type Group
    nodes::Set{Int}
    inner_prob::Float64
    exit_prob::Float64
end

type Partition{V}
    flowgraph::FlowGraph{V}
    membership::Vector{Int}
    community::Dict{Int,Group}
end

function init_partition{V}(fg::FlowGraph{V}, membership::Vector{Int})
    maximum(membership) ≤ num_vertices(fg.graph) || error("maximum(membership) must less than num_vertices(graph)")
    minimum(membership) > 0 || error("value of membership must be positive integer")
    community = Dict{Int,Group}()
    for u in vertices(fg.graph)
        u_idx = vertex_index(u, fg.graph)
        comm_idx = membership[u_idx]
        if haskey(community, comm_idx)
            push!(community[comm_idx].nodes, u_idx)
            community[comm_idx].inner_prob += fg.visit_prob[u_idx]
            for e in out_edges(u, fg.graph)
                e_idx = edge_index(e, fg.graph)
                v = target(e, fg.graph)
                v_idx = vertex_index(v, fg.graph)
                if !in(v_idx, community[comm_idx].nodes)
                    community[comm_idx].exit_prob += fg.trans_prob[e_idx]
                end
            end
            for e in in_edges(u, fg.graph)
                e_idx = edge_index(e, fg.graph)
                v = source(e, fg.graph)
                v_idx = vertex_index(v, fg.graph)
                if in(v_idx, community[comm_idx].nodes)
                    community[comm_idx].exit_prob -= fg.trans_prob[e_idx]
                end
            end
        else
            exit_prob = 0.0
            for e in out_edges(u, fg.graph)
                e_idx = edge_index(e, fg.graph)
                exit_prob += fg.trans_prob[e_idx]
            end
            community[comm_idx] = Group(Set(u_idx), fg.visit_prob[u_idx], exit_prob)
        end
    end
    Partition(fg, membership, community)
end

function init_partition{V}(fg::FlowGraph{V})
    n = num_vertices(fg.graph)
    membership = collect(1:n)
    groups = Array(Group, n)
    for u in vertices(fg.graph)
        u_idx = vertex_index(u, fg.graph)
        exit_prob = 0.0
        for e in out_edges(u, fg.graph)
            e_idx = edge_index(e, fg.graph)
            exit_prob += fg.trans_prob[e_idx]
        end
        groups[u_idx] = Group(Set(u_idx), fg.visit_prob[u_idx], exit_prob)
    end
    community = Dict{Int,Group}(zip(membership, groups))
    Partition(fg, membership, community)
end

function renumber_communities!{V}(partition::Partition{V})
    csizes = Int[length(partition.community[i].nodes) for i in keys(partition.community)]
    partition.community = Dict{Int, Group}(zip(sortperm(csizes), values(partition.community)))
    for (i, group) in partition.community
        for u_idx in group.nodes
            partition.membership[u_idx] = i
        end
    end
end

1/2
