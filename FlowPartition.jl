using Graphs
import Graphs.graph

type FlowGroup{V}
    nodes::Set{V}
    inner_prob::Float64
    exit_prob::Float64
end

type FlowPartition{V} <: AbstractPartition{V}
    flowgraph::FlowGraph{V}
    membership::Vector{Int}
    community::Dict{Int,FlowGroup{V}}
    total_exit_prob::Float64
end

# construction
function flow_partition{V}(fg::FlowGraph{V})
    n = num_vertices(fg.graph)
    membership = collect(1:n)
    groups = Array(FlowGroup{V}, n)

    for u in vertices(fg.graph)
        u_idx = vertex_index(u, fg.graph)
        exit_prob = 0.0
        for e in out_edges(u, fg.graph)
            e_idx = edge_index(e, fg.graph)
            exit_prob += fg.trans_prob[e_idx]
        end
        groups[u_idx] = FlowGroup(Set(u), fg.visit_prob[u_idx], exit_prob)
    end

    community = Dict{Int,FlowGroup{V}}(zip(membership, groups))

    total_exit_prob = 0.0
    for group in values(community)
        total_exit_prob += group.exit_prob
    end

    FlowPartition(fg, membership, community, total_exit_prob)
end

function flow_partition{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T})
    fg = flow_graph(g, weights)
    flow_partition(fg)
end

function flow_partition{V}(fg::FlowGraph{V}, membership::Vector{Int})
    maximum(membership) ≤ num_vertices(fg.graph) || error("maximum(membership) must less than num_vertices(g)")
    minimum(membership) > 0 || error("value of membership must be positive integer")

    community = Dict{Int,FlowGroup{V}}()

    for u in vertices(g)
        u_idx = vertex_index(u, g)
        comm_idx = membership[u_idx]
        if haskey(community, comm_idx)
            push!(community[comm_idx].nodes, u)
            community[comm_idx].inner_prob += fg.visit_prob[u_idx]
            for e in out_edges(u, g)
                e_idx = edge_index(e, g)
                v = target(e, g)
                v_idx = vertex_index(v, g)
                if !in(v, community[comm_idx].nodes)
                    community[comm_idx].exit_prob += fg.trans_prob[e_idx]
                else
                    community[comm_idx].exit_prob -= fg.trans_prob[e_idx]
                end
            end
        else
            exit_prob = 0.0
            for e in out_edges(u, g)
                e_idx = edge_index(e, g)
                exit_prob += fg.trans_prob[e_idx]
            end
            community[comm_idx] = FlowGroup(Set(u), fg.visit_prob[u_idx], exit_prob)
        end
    end

    total_exit_prob = 0.0
    for group in values(community)
        total_exit_prob += group.exit_prob
    end

    FlowPartition(fg, membership, community, total_exit_prob)
end

function flow_partition{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T}, membership::Vector{Int})
    maximum(membership) ≤ num_vertices(g) || error("maximum(membership) must less than num_vertices(g)")
    minimum(membership) > 0 || error("value of membership must be positive integer")

    fg = flow_graph(g, weights)
    flow_partition(fg, membership)
end

# require interface
graph{V}(fp::FlowPartition{V}) = fp.flowgraph.graph
membership{V}(fp::FlowPartition{V}) = fp.membership
membership{V}(fp::FlowPartition{V}, u::V) = fp.membership[vertex_index(u, fp.flowgraph.graph)]

# mutation
"update partition based on the membership vector"
function update_partition!{V}(fp::FlowPartition{V}, membership::Vector{Int})
    maximum(membership) ≤ num_vertices(fp.flowgraph.graph) || error("maximum(membership) must less than num_vertices(graph)")
    minimum(membership) > 0 || error("value of membership must be positive integer")

    # update membership vector and communities
    fp.community = Dict{Int,FlowGroup{V}}()
    for u in vertices(fp.flowgraph.graph)
        u_idx = vertex_index(u, fp.flowgraph.graph)
        fp.membership[u_idx] = membership[u_idx]
        comm_idx = membership[u_idx]
        if haskey(fp.community, comm_idx)
            push!(fp.community[comm_idx].nodes, u)
            fp.community[comm_idx].inner_prob += fp.flowgraph.visit_prob[u_idx]
            for e in out_edges(u, fp.flowgraph.graph)
                e_idx = edge_index(e, fp.flowgraph.graph)
                v = target(e, fp.flowgraph.graph)
                v_idx = vertex_index(v, fp.flowgraph.graph)
                if !in(v, fp.community[comm_idx].nodes)
                    fp.community[comm_idx].exit_prob += fp.flowgraph.trans_prob[e_idx]
                else
                    fp.community[comm_idx].exit_prob -= fp.flowgraph.trans_prob[e_idx]
                end
            end
        else
            exit_prob = 0.0
            for e in out_edges(u, fp.flowgraph.graph)
                e_idx = edge_index(e, fp.flowgraph.graph)
                exit_prob += fp.flowgraph.trans_prob[e_idx]
            end
            fp.community[comm_idx] = FlowGroup(Set(u), fp.flowgraph.visit_prob[u_idx], exit_prob)
        end
    end

    fp.total_exit_prob = 0.0
    for group in values(fp.community)
        fp.total_exit_prob += group.exit_prob
    end
end

"Renumber the communities so that they are numbered 0,...,q-1 where q is the number of communities."
function renumber_communities!{V}(fp::FlowPartition{V})
    csizes = Int[length(fp.community[i].nodes) for i in keys(fp.community)]
    perm_idx = sortperm(csizes, rev=true)
    fp.community = Dict{Int, FlowGroup{V}}(zip(enumerate(collect(values(fp.community))[perm_idx])))
    for (i, group) in fp.community
        for u in group.nodes
            u_idx = vertex_index(u, fp.flowgraph.graph)
            fp.membership[u_idx] = i
        end
    end
end

"Move a node to a new community and update the partition, this also removes any empty communities."
function move_node!{V}(partition::FlowPartition{V}, u::V, new_comm::Int)
    haskey(partition.community, new_comm) || error("partition has no community $new_comm")

    u_idx = vertex_index(u, partition.flowgraph.graph)
    old_comm = partition.membership[u_idx]

    # just skip if move to a community the node already in
    if new_comm != old_comm
        # remove node from old community
        delete!(partition.community[old_comm].nodes, u)
        # change of visit probability of the old community
        partition.community[old_comm].inner_prob -= partition.flowgraph.visit_prob[u_idx]

        # add node to new community
        push!(partition.community[new_comm].nodes, u)
        # change of inner probability of the new community
        partition.community[new_comm].inner_prob += partition.flowgraph.visit_prob[u_idx]

        for e in out_edges(u, partition.flowgraph.graph)
            e_idx = edge_index(e, partition.flowgraph.graph)
            v = target(e, partition.flowgraph.graph)
            # change of exit probability of the old community
            if !in(v, partition.community[old_comm].nodes)
                partition.community[old_comm].exit_prob -= partition.flowgraph.trans_prob[e_idx]
                partition.total_exit_prob -= partition.flowgraph.trans_prob[e_idx]
            else
                partition.community[old_comm].exit_prob += partition.flowgraph.trans_prob[e_idx]
                partition.total_exit_prob += partition.flowgraph.trans_prob[e_idx]
            end
            # change of exit probability of the new community
            if !in(v, partition.community[new_comm].nodes)
                partition.community[new_comm].exit_prob += partition.flowgraph.trans_prob[e_idx]
                partition.total_exit_prob += partition.flowgraph.trans_prob[e_idx]
            else
                partition.community[new_comm].exit_prob -= partition.flowgraph.trans_prob[e_idx]
                partition.total_exit_prob -= partition.flowgraph.trans_prob[e_idx]
            end
        end

        # if the old community is empty after remove node u, we remove it
        if isempty(partition.community[old_comm].nodes)
            delete!(partition.community, old_comm)
        end

        # update the membership vector
        partition.membership[u_idx] = new_comm
    end
end

"Read new communities from coarser partition assuming that the community
 represents a node in the coarser partition (with the same index as the
 community number)."
function from_coarser_partition!{V}(partition::FlowPartition{V}, coarser_partition::FlowPartition{V})
    for u in vertices(partition.flowgraph.graph)
        u_idx = vertex_index(u, partition.flowgraph.graph)
        # what is the community of the node
        u_comm_level1 = partition.membership[u_idx]

        # In the coarser partition, the node should have the community id
        # so that the community of that node gives the coarser community.
        u_comm_level2 = coarser_partition.membership[u_comm_level1]
        partition.membership[u_idx] = u_comm_level2
    end
    update_partition!(partition, partition.membership)
end

"Read new partition from another partition."
function from_partition!{V}(partition::FlowPartition{V}, other_partition::FlowPartition{V})
    #Assign the membership of every node in the supplied partition to the one in this partition
    partition.membership[:] = other_partition.membership[:]
    update_partition!(partition, partition.membership)
end

"Calculate what is the total weight going from/to a node to a community."
function weight_to_comm{V}(partition::FlowPartition{V}, u::V, comm::Int)
    total_weight = 0.0
    u_idx = vertex_index(u, partition.flowgraph.graph)
    u_comm = partition.membership[u_idx]
    if u_comm != comm
        for e in out_edges(u, partition.flowgraph.graph)
            e_idx = edge_index(e, partition.flowgraph.graph)
            v = target(e, partition.flowgraph.graph)
            if in(v, partition.community[comm].nodes)
                total_weight += partition.flowgraph.trans_prob[e_idx]
            end
        end
    end
    total_weight
end
weight_from_comm{V}(partition::FlowPartition{V}, u::V, comm::Int) = weight_to_comm(partition, u, comm)

"get neighbor communities of node u"
function get_neigh_comms{V}(partition::FlowPartition{V}, u::V)
    neigh_comms = Set{Int}()
    for v in out_neighbors(u, partition.flowgraph.graph)
        v_idx = vertex_index(v, partition.flowgraph.graph)
        push!(neigh_comms, partition.membership[v_idx])
    end
    neigh_comms
end

# retrun p*log(p)
plogp(x) = x > 0.0 ? x*log(x) : 0.0
plogp(xs::Vector{Float64}) = Float64[plogp(x) for x in xs]

"Returns the difference in average decribe length if we move a node to a new community"
function diff_move{V}(partition::FlowPartition{V}, u::V, new_comm::Int)
    u_idx = vertex_index(u, partition.flowgraph.graph)
    old_comm = partition.membership[u_idx]

    δL = 0.0

    #println(partition.community[old_comm].exit_prob)

    if new_comm != old_comm
        δexit_prob_old_comm = 0.0
        δexit_prob_new_comm = 0.0
        for e in out_edges(u, partition.flowgraph.graph)
            e_idx = edge_index(e, partition.flowgraph.graph)
            v = target(e, partition.flowgraph.graph)
            if in(v, partition.community[old_comm].nodes)
                δexit_prob_old_comm += partition.flowgraph.trans_prob[e_idx]
            else
                δexit_prob_old_comm -= partition.flowgraph.trans_prob[e_idx]
            end
            if in(v, partition.community[new_comm].nodes)
                δexit_prob_new_comm -= partition.flowgraph.trans_prob[e_idx]
            else
                δexit_prob_new_comm += partition.flowgraph.trans_prob[e_idx]
            end
        end
        δtotal_exit_prob = δexit_prob_old_comm + δexit_prob_new_comm
        #println("*******")
        #println(δtotal_exit_prob, ", ", δexit_prob_old_comm, ", ", δexit_prob_new_comm)
        #println(partition.total_exit_prob)

        δL1 = plogp(partition.total_exit_prob + δtotal_exit_prob) - plogp(partition.total_exit_prob)
        δL2 = -2(plogp(partition.community[old_comm].exit_prob + δexit_prob_old_comm) - plogp(partition.community[old_comm].exit_prob))
        δL3 = -2(plogp(partition.community[new_comm].exit_prob + δexit_prob_new_comm) - plogp(partition.community[new_comm].exit_prob))
        δL4 = plogp(partition.community[old_comm].exit_prob + δexit_prob_old_comm + partition.community[old_comm].inner_prob - partition.flowgraph.visit_prob[u_idx]) -
            plogp(partition.community[old_comm].exit_prob + partition.community[old_comm].inner_prob)
        δL5 = plogp(partition.community[new_comm].exit_prob + δexit_prob_new_comm + partition.community[new_comm].inner_prob + partition.flowgraph.visit_prob[u_idx]) -
            plogp(partition.community[new_comm].exit_prob + partition.community[new_comm].inner_prob)

        #println("############")
        #println("1: ", partition.community[new_comm].exit_prob + δexit_prob_new_comm + partition.community[new_comm].inner_prob + partition.flowgraph.visit_prob[u_idx])
        println("3: ", plogp(partition.community[new_comm].exit_prob + δexit_prob_new_comm + partition.community[new_comm].inner_prob + partition.flowgraph.visit_prob[u_idx]))
        δL = δL1 + δL2 + δL3 + δL4 + δL5
    end

    δL
end

"Give the average decribe length of the partition."
function quality{V}(partition::FlowPartition{V})
    L1 = plogp(partition.total_exit_prob)
    L2 = -2sum(plogp([partition.community[i].exit_prob for i in keys(partition.community)]))
    L3 = -sum(plogp(partition.flowgraph.visit_prob))
    L4 = sum(plogp([partition.community[i].exit_prob+partition.community[i].inner_prob for i in keys(partition.community)]))

    L1 + L2 + L3 + L4
end
