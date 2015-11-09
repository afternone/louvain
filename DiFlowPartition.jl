using Graphs
import Graphs.graph

type DiFlowGroup{V}
    nodes::Set{V}
    inner_prob::Float64
    exit_prob::Float64
    iso_prob::Float64
end

type DiFlowPartition{V} <: AbstractPartition{V}
    flowgraph::DiFlowGraph{V}
    membership::Vector{Int}
    community::Dict{Int,DiFlowGroup{V}}
    total_exit_prob::Float64
end

g = simple_graph(3)
add_edge!(g,1,2)

fg = diflowgraph(g,τ=0.85)
fp = diflow_partition(g)
fp.community
m = [1,1,1]
fp1 = diflow_partition(fg, m)
fp1.community
fg.visit_prob

# construction
function diflow_partition{V}(fg::DiFlowGraph{V})
    n = num_vertices(fg.graph)
    membership = collect(1:n)
    groups = Array(DiFlowGroup{V}, n)

    for u in vertices(fg.graph)
        u_idx = vertex_index(u, fg.graph)
        if out_degree(u, fg.graph) > 0
            exit_prob = 0.0
            for e in out_edges(u, fg.graph)
                e_idx = edge_index(e, fg.graph)
                exit_prob += fg.tau*fg.visit_prob[u_idx]*(n-1)/n + (1.0-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
            end
            groups[u_idx] = DiFlowGroup(Set(u), fg.visit_prob[u_idx], exit_prob, 0.0)
        else
            groups[u_idx] = DiFlowGroup(Set(u), fg.visit_prob[u_idx], fg.visit_prob[u_idx]*(n-1)/n, fg.visit_prob[u_idx])
        end
    end

    community = Dict{Int,DiFlowGroup{V}}(zip(membership, groups))

    total_exit_prob = 0.0
    for group in values(community)
        total_exit_prob += group.exit_prob
    end

    DiFlowPartition(fg, membership, community, total_exit_prob)
end

function diflow_partition{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T}=ones(num_edges(g)); τ=0.15)
    fg = diflowgraph(g, weights, τ=τ)
    diflow_partition(fg)
end


function diflow_partition{V}(fg::DiFlowGraph{V}, membership::Vector{Int}; τ=0.15)
    maximum(membership) ≤ num_vertices(fg.graph) || error("maximum(membership) must less than num_vertices(g)")
    minimum(membership) > 0 || error("value of membership must be positive integer")

    g = fg.graph
    n = num_vertices(g)
    community = Dict{Int,DiFlowGroup{V}}()

    for u in vertices(g)
        u_idx = vertex_index(u, g)
        comm_idx = membership[u_idx]
        if haskey(community, comm_idx)
            ni = length(community[comm_idx].nodes)
            community[comm_idx].exit_prob += fg.tau*(n-ni-1)*fg.visit_prob[u_idx]/n - fg.tau*community[comm_idx].inner_prob/n -
                (1-fg.tau)*community[comm_idx].iso_prob/n
            push!(community[comm_idx].nodes, u)
            community[comm_idx].inner_prob += fg.visit_prob[u_idx]
            if out_degree(u, g) > 0
                for e in out_edges(u, g)
                    e_idx = edge_index(e, g)
                    v = target(e, g)
                    v_idx = vertex_index(v, g)
                    if !in(v, community[comm_idx].nodes)
                        community[comm_idx].exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                    end
                end
                for e in in_edges(u, g)
                    e_idx = edge_index(e, g)
                    v = source(e, g)
                    v_idx = vertex_index(v, g)
                    if in(v, community[comm_idx].nodes)
                        community[comm_idx].exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    end
                end
            else
                for e in in_edges(u, g)
                    e_idx = edge_index(e, g)
                    v = source(e, g)
                    v_idx = vertex_index(v, g)
                    if in(v, community[comm_idx].nodes)
                        community[comm_idx].exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    end
                end
                community[comm_idx].iso_prob += fg.visit_prob[u_idx]
                community[comm_idx].exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*(n-ni-1)/n
            end
        else
            if out_degree(u, fg.graph) > 0
                exit_prob = fg.tau*fg.visit_prob[u_idx]*(n-1)/n
                for e in out_edges(u, fg.graph)
                    e_idx = edge_index(e, fg.graph)
                    exit_prob += (1.0-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                end
                community[comm_idx] = DiFlowGroup(Set(u), fg.visit_prob[u_idx], exit_prob, 0.0)
            else
                community[comm_idx] = DiFlowGroup(Set(u), fg.visit_prob[u_idx], fg.visit_prob[u_idx]*(n-1)/n, fg.visit_prob[u_idx])
            end
        end
    end

    total_exit_prob = 0.0
    for group in values(community)
        total_exit_prob += group.exit_prob
    end

    DiFlowPartition(fg, membership, community, total_exit_prob)
end

function flow_partition{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T}, membership::Vector{Int})
    maximum(membership) ≤ num_vertices(g) || error("maximum(membership) must less than num_vertices(g)")
    minimum(membership) > 0 || error("value of membership must be positive integer")

    fg = diflowgraph(g, weights)
    diflow_partition(fg, membership)
end

function update_partition!{V}(partition::DiFlowPartition{V})
    maximum(partition.membership) ≤ num_vertices(partition.flowgraph.graph) || error("maximum(membership) must less than num_vertices(g)")
    minimum(partition.membership) > 0 || error("value of membership must be positive integer")

    fg = partition.flowgraph
    g = fg.graph

    n = num_vertices(g)
    empty!(partition.community)

    for u in vertices(g)
        u_idx = vertex_index(u, g)
        comm_idx = partition.membership[u_idx]
        if haskey(partition.community, comm_idx)
            ni = length(partition.community[comm_idx].nodes)
            partition.community[comm_idx].exit_prob += fg.tau*(n-ni-1)*fg.visit_prob[u_idx]/n -
                fg.tau*partition.community[comm_idx].inner_prob/n -
                (1-fg.tau)*partition.community[comm_idx].iso_prob/n
            push!(partition.community[comm_idx].nodes, u)
            partition.community[comm_idx].inner_prob += fg.visit_prob[u_idx]
            if out_degree(u, g) > 0
                for e in out_edges(u, g)
                    e_idx = edge_index(e, g)
                    v = target(e, g)
                    v_idx = vertex_index(v, g)
                    if !in(v, partition.community[comm_idx].nodes)
                        partition.community[comm_idx].exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                    end
                end
                for e in in_edges(u, g)
                    e_idx = edge_index(e, g)
                    v = source(e, g)
                    v_idx = vertex_index(v, g)
                    if in(v, partition.community[comm_idx].nodes)
                        partition.community[comm_idx].exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    end
                end
            else
                for e in in_edges(u, g)
                    e_idx = edge_index(e, g)
                    v = source(e, g)
                    v_idx = vertex_index(v, g)
                    if in(v, partition.community[comm_idx].nodes)
                        partition.community[comm_idx].exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    end
                end
                partition.community[comm_idx].iso_prob += fg.visit_prob[u_idx]
                partition.community[comm_idx].exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*(n-ni-1)/n
            end
        else
            if out_degree(u, fg.graph) > 0
                exit_prob = fg.tau*fg.visit_prob[u_idx]*(n-1)/n
                for e in out_edges(u, fg.graph)
                    e_idx = edge_index(e, fg.graph)
                    exit_prob += (1.0-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                end
                partition.community[comm_idx] = DiFlowGroup(Set(u), fg.visit_prob[u_idx], exit_prob, 0.0)
            else
                partition.community[comm_idx] = DiFlowGroup(Set(u), fg.visit_prob[u_idx], fg.visit_prob[u_idx]*(n-1)/n, fg.visit_prob[u_idx])
            end
        end
    end

    partition.total_exit_prob = 0.0
    for group in values(partition.community)
        partition.total_exit_prob += group.exit_prob
    end
end

# require interface
graph{V}(fp::DiFlowPartition{V}) = fp.flowgraph.graph
membership{V}(fp::DiFlowPartition{V}) = fp.membership
membership{V}(fp::DiFlowPartition{V}, u::V) = fp.membership[vertex_index(u, fp.flowgraph.graph)]
community{V}(fp::DiFlowPartition{V}) = keys(fp.community)
is_right_direction{V}(fp::DiFlowPartition{V}, quality, new_quality) = new_quality < quality
still_running{V}(fp::DiFlowPartition{V}, once_diff, ϵ) = once_diff < -ϵ

"Renumber the communities so that they are numbered 0,...,q-1 where q is the number of communities."
function renumber_communities!{V}(fp::DiFlowPartition{V})
    csizes = Int[length(fp.community[i].nodes) for i in keys(fp.community)]
    perm_idx = sortperm(csizes, rev=true)
    fp.community = Dict{Int, DiFlowGroup{V}}(zip(enumerate(collect(values(fp.community))[perm_idx])))
    for (i, group) in fp.community
        for u in group.nodes
            u_idx = vertex_index(u, fp.flowgraph.graph)
            fp.membership[u_idx] = i
        end
    end
end

"Move a node to a new community and update the partition, this also removes any empty communities."
function move_node!{V}(partition::DiFlowPartition{V}, u::V, new_comm::Int)
    haskey(partition.community, new_comm) || error("partition has no community $new_comm")

    fg = partition.flowgraph
    g = fg.graph
    n = num_vertices(g)
    u_idx = vertex_index(u, g)
    old_comm = partition.membership[u_idx]

    # just skip if move to a community the node already in
    if new_comm != old_comm
        n_old = length(partition.community[old_comm].nodes)
        n_new = length(partition.community[new_comm].nodes)


        partition.community[old_comm].exit_prob += -fg.tau*(n-n_old+1)*fg.visit_prob[u_idx]/n +
                fg.tau*partition.community[old_comm].inner_prob/n +
                (1-fg.tau)*partition.community[old_comm].iso_prob/n
        partition.total_exit_prob += -fg.tau*(n-n_old+1)*fg.visit_prob[u_idx]/n +
                fg.tau*partition.community[old_comm].inner_prob/n +
                (1-fg.tau)*partition.community[old_comm].iso_prob/n
        partition.community[new_comm].exit_prob += fg.tau*(n-n_new-1)*fg.visit_prob[u_idx]/n -
                fg.tau*partition.community[new_comm].inner_prob/n -
                (1-fg.tau)*partition.community[new_comm].iso_prob/n
        partition.total_exit_prob += fg.tau*(n-n_new-1)*fg.visit_prob[u_idx]/n -
                fg.tau*partition.community[new_comm].inner_prob/n -
                (1-fg.tau)*partition.community[new_comm].iso_prob/n

        # remove node from old community
        delete!(partition.community[old_comm].nodes, u)
        # change of visit probability of the old community
        partition.community[old_comm].inner_prob -= partition.flowgraph.visit_prob[u_idx]

        # add node to new community
        push!(partition.community[new_comm].nodes, u)
        # change of inner probability of the new community
        partition.community[new_comm].inner_prob += partition.flowgraph.visit_prob[u_idx]

        if out_degree(u, g) > 0
            for e in out_edges(u, g)
                e_idx = edge_index(e, g)
                v = target(e, g)
                v_idx = vertex_index(v, g)
                if !in(v, partition.community[new_comm].nodes)
                    partition.community[new_comm].exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                    partition.total_exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                else
                    partition.community[old_comm].exit_prob -= (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                    partition.total_exit_prob -= (1-fg.tau)*fg.visit_prob[u_idx]*fg.trans_prob[e_idx]
                end
            end
            for e in in_edges(u, g)
                e_idx = edge_index(e, g)
                v = source(e, g)
                v_idx = vertex_index(v, g)
                if in(v, partition.community[new_comm].nodes)
                    partition.community[new_comm].exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    partition.total_exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                else
                    partition.community[old_comm].exit_prob += (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    partition.total_exit_prob += (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                end
            end
        else
            for e in in_edges(u, g)
                e_idx = edge_index(e, g)
                v = source(e, g)
                v_idx = vertex_index(v, g)
                if in(v, partition.community[new_comm].nodes)
                    partition.community[new_comm].exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    partition.total_exit_prob -= (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                else
                    partition.community[old_comm].exit_prob += (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                    partition.taotal_exit_prob += (1-fg.tau)*fg.visit_prob[v_idx]*fg.trans_prob[e_idx]
                end
            end
            partition.community[new_comm].iso_prob += fg.visit_prob[u_idx]
            partition.community[old_comm].iso_prob -= fg.visit_prob[u_idx]
            partition.community[new_comm].exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*(n-n_new-1)/n
            partition.total_exit_prob += (1-fg.tau)*fg.visit_prob[u_idx]*(n-n_new-1)/n
            partition.community[old_comm].exit_prob += -(1-fg.tau)*fg.visit_prob[u_idx]*(n-n_old+1)/n
            partition.total_exit_prob += -(1-fg.tau)*fg.visit_prob[u_idx]*(n-n_old+1)/n
        end

        # if the old community is empty after remove node u, we remove it
        if isempty(partition.community[old_comm].nodes)
            delete!(partition.community, old_comm)
        end

        # update the membership vector
        partition.membership[u_idx] = new_comm
    end
end

fp.membership=[1,2,3]
update_partition!(fp)
fp.community
move_node!(fp, 3, 2)
fp.total_exit_prob
fp1 = diflow_partition(fg, [2,2,2])
fp1.community
