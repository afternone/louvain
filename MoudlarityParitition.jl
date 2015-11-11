using Graphs

type ModularityGroup{V}
    nodes::Set{V}
    csize::Int
    inner_weight::Float64
    in_weight::Float64
    out_weight::Float64
end

type ModularityPartition{V} <: AbstractPartition{V}
    mod_graph::ModularityGraph{V}
    membership::Vector{Int}
    community::Dict{Int,ModularityGroup{V}}
    total_weight_in::Float64 # total weight in all community
    total_possible_edges_in::Int
end

function collapse_partition1{V}(partition::MPartition{V})
    num_comm = length(partition.community)
    collapsed_trans_prob = Dict{Int,Float64}[]
    for i=1:num_comm
        push!(collapsed_trans_prob, Dict{Int,Float64}())
    end

    for e in edges(partition.mgraph.graph)
        e_idx = edge_index(e, partition.mgraph.graph)
        u = source(e, partition.mgraph.graph)
        v = target(e, partition.mgraph.graph)
        u_idx = vertex_index(u, partition.mgraph.graph)
        v_idx = vertex_index(v, partition.mgraph.graph)
        u_comm = partition.membership[u_idx]
        v_comm = partition.membership[v_idx]

        w = partition.mgraph.edge_weights[e_idx]
        if haskey(collapsed_trans_prob[u_comm], v_comm)
            collapsed_trans_prob[u_comm][v_comm] += w
        else
            collapsed_trans_prob[u_comm][v_comm] = w
        end
    end

    graph = simple_graph(num_comm, is_directed=is_directed(partition.mgraph.graph))
    graph_edge_weights = Float64[]
    graph_node_sizes = Array(Int, num_comm)

    for u_comm=1:num_comm
        graph_node_sizes[u_comm] = partition.community[u_comm].csize
        for (v_comm, w) in collapsed_trans_prob[u_comm]
            add_edge!(graph, u_comm, v_comm)
            push!(graph_edge_weights, w)
        end
    end

    mpartition(graph, graph_edge_weights, graph_node_sizes, partition.mgraph.correct_self_loops)
end

function move_node1!{V}(mp::MPartition{V}, u::V, new_comm::Int)
    mg = mp.mgraph
    g = mg.graph
    u_idx = vertex_index(u,g)
    node_size = mg.node_sizes[u_idx]
    old_comm = mp.membership[u_idx]
    old_csize = mp.community[old_comm].csize
    new_csize = mp.community[new_comm].csize
    mp.total_possible_edges_in_all_comms += 2.0*node_size*(new_csize - old_csize + node_size)/(2.0 - Float64(is_directed(g)))

    # remove from old community
    delete!(mp.community[old_comm].nodes, u)
    mp.community[old_comm].csize -= node_size

    # add to new community
    push!(mp.community[new_comm].nodes, u)
    mp.community[new_comm].csize += node_size

    for e in out_edges(u,g)
        e_idx = edge_index(e,g)
        v = target(e,g)
        v_idx = vertex_index(v,g)
        v_comm = mp.membership[v_idx]
        w = mg.edge_weights[e_idx]
        int_weight = w/(is_directed(g) ? 1.0 : 2.0)
        if in(v, mp.community[old_comm].nodes)
            mp.community[old_comm].weight_in += w
            mp.community[old_comm].weight_inner -= int_weight
            mp.total_weight_in_all_comms -= int_weight
        else
            mp.community[old_comm].weight_out -= w
        end
        if in(v, mp.community[new_comm].nodes)
            mp.community[new_comm].weight_in -= w
            mp.community[new_comm].weight_inner += int_weight
            mp.total_weight_in_all_comms += int_weight
        else
            mp.community[new_comm].weight_out += w
        end
    end

    for e in in_edges(u,g)
        e_idx = edge_index(e,g)
        v = source(e,g)
        v_idx = vertex_index(v,g)
        v_comm = mp.membership[v_idx]
        w = mg.edge_weights[e_idx]
        int_weight = w/(is_directed(g) ? 1.0 : 2.0)
        if in(v, mp.community[old_comm].nodes)
            mp.community[old_comm].weight_out += w
            mp.community[old_comm].weight_inner -= int_weight
            mp.total_weight_in_all_comms -= int_weight
        else
            mp.community[old_comm].weight_in -= w
        end
        if in(v, mp.community[new_comm].nodes)
            mp.community[new_comm].weight_out -= w
            mp.community[new_comm].weight_inner += int_weight
            mp.total_weight_in_all_comms += int_weight
        else
            mp.community[new_comm].weight_in += w
        end
    end

    # if the old community is empty after remove node u, we remove it
    if isempty(mp.community[old_comm].nodes)
        delete!(mp.community, old_comm)
    end

    # update the membership vector
    mp.membership[u_idx] = new_comm
end





# construction
function modularity_partition{V}(mg::ModularityGraph{V})
    g = mg.graph
    weights = mg.edge_weights
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
        community[u_idx] = ModularityGroup(Set(u), mg.node_size[u_idx], mg.node_self_weights[u_idx], in_weight, out_weight)
    end
    total_weight_in = sum(mg.node_self_weights)
    total_weight = sum(weights)
    total_possible_edges_in = 0
    for n_c in mg.node_size
        possible_edges = 0
        if mg.correct_self_loops
            possible_edges = Int(n_c*n_c/(2 - Int(is_directed(g))))
        else
            possible_edges = Int(n_c*(n_c-1)/(2 - Int(is_directed(g))))
        end
        total_possible_edges_in += possible_edges
    end

    ModularityPartition(mg, membership, community, total_weight_in, total_possible_edges_in)
end

function update_partition!{V}(partition::ModularityPartition{V})
    mg = partition.mod_graph
    g = mg.graph
    membership = partition.membership

    empty!(partition.community)
    partition.total_weight_in = 0.0

    for u in vertices(g)
        u_idx = vertex_index(u,g)
        comm_idx = membership[u_idx]
        if haskey(partition.community, comm_idx)
            push!(partition.community[comm_idx].nodes, u)
            partition.community[comm_idx].csize += mg.node_size[u_idx]
            for e in out_edges(u,g)
                e_idx = edge_index(e,g)
                v = target(e,g)
                v_idx = vertex_index(v,g)
                v_comm = membership[v_idx]
                if !in(v, partition.community[comm_idx].nodes)
                    partition.community[comm_idx].out_weight += mg.edge_weights[e_idx]
                else
                    partition.community[comm_idx].inner_weight += is_directed(g) ? mg.edge_weights[e_idx] : mg.edge_weights[e_idx]/2
                    partition.total_weight_in += is_directed(g) ? mg.edge_weights[e_idx] : mg.edge_weights[e_idx]/2
                    partition.community[comm_idx].in_weight -= mg.edge_weights[e_idx]
                end
            end
            for e in in_edges(u,g)
                e_idx = edge_index(e,g)
                v = source(e,g)
                v_idx = vertex_index(v,g)
                v_comm = membership[v_idx]
                if in(v, partition.community[comm_idx].nodes)
                    partition.community[comm_idx].out_weight -= mg.edge_weights[e_idx]
                    partition.community[comm_idx].inner_weight += is_directed(g) ? mg.edge_weights[e_idx] : mg.edge_weights[e_idx]/2
                    partition.total_weight_in += is_directed(g) ? mg.edge_weights[e_idx] : mg.edge_weights[e_idx]/2
                else
                    partition.community[comm_idx].in_weight += mg.edge_weights[e_idx]
                end
            end
        else
            out_weight = 0.0
            for e in out_edges(u, g)
                e_idx = edge_index(e, g)
                out_weight += mg.edge_weights[e_idx]
            end
            in_weight = 0.0
            for e in in_edges(u, g)
                e_idx = edge_index(e, g)
                in_weight += mg.edge_weights[e_idx]
            end
            partition.community[comm_idx] = ModularityGroup(Set(u), mg.node_size[u_idx], 0.0, in_weight, out_weight)
        end
    end
end


# require interface
graph{V}(mp::ModularityPartition{V}) = mp.mod_graph.graph
membership{V}(mp::ModularityPartition{V}) = mp.membership
membership{V}(mp::ModularityPartition{V}, u::V) = mp.membership[vertex_index(u, mp.mod_graph.graph)]
community{V}(mp::ModularityPartition{V}) = keys(mp.community)
is_right_direction{V}(mp::ModularityPartition{V}, quality, new_quality) = new_quality > quality
still_running{V}(mp::ModularityPartition{V}, once_diff, ϵ) = once_diff > ϵ

"Renumber the communities so that they are numbered 0,...,q-1 where q is the number of communities."
function renumber_communities!{V}(mp::ModularityPartition{V})
    csizes = Int[length(mp.community[i].nodes) for i in keys(mp.community)]
    perm_idx = sortperm(csizes, rev=true)
    mp.community = Dict{Int, ModularityGroup{V}}(zip(enumerate(collect(values(mp.community))[perm_idx])))
    for (i, group) in mp.community
        for u in group.nodes
            u_idx = vertex_index(u, mp.mod_graph.graph)
            mp.membership[u_idx] = i
        end
    end
end




# TO DO
function move_node!{V}(mp::ModularityPartition{V}, u::V, new_comm::Int)
    haskey(mp.community, new_comm) || error("partition has no community $new_comm")
    mg = mp.mod_graph
    g = mg.graph
    u_idx = vertex_index(u, g)
    node_size = mg.node_size[u_idx]
    old_comm = mp.membership[u_idx]

    # just skip if move to a community the node already in
    if new_comm != old_comm
        mp.total_possible_edges_in += Int(2*node_size*(mp.community[new_comm].csize - mp.community[old_comm].csize + node_size)/(2 - Int(is_directed(g))))
        # remove node from old community
        delete!(mp.community[old_comm].nodes, u)
        mp.community[old_comm].csize -= node_size
        # add node to new community
        push!(mp.community[new_comm].nodes, u)
        mp.community[new_comm].csize += node_size

        for e in out_edges(u,g)
            e_idx = edge_index(e,g)
            v = target(e,g)
            v_idx = vertex_index(v,g)
            v_comm = mp.membership[v_idx]
            if !in(v, mp.community[old_comm].nodes)
                mp.community[old_comm].out_weight -= mg.edge_weights[e_idx]
            else
                w = is_directed(g) ? mg.edge_weights[e_idx] : mg.edge_weights[e_idx]/2
                mp.community[old_comm].inner_weight -= w
                mp.total_weight_in -= w
                mp.community[old_comm].in_weight += mg.edge_weights[e_idx]
            end
            if !in(v, mp.community[new_comm].nodes)
                mp.community[new_comm].out_weight += mg.edge_weights[e_idx]
            else
                w = is_directed(g) ? mg.edge_weights[e_idx] : mg.edge_weights[e_idx]/2
                mp.community[new_comm].inner_weight += w
                mp.total_weight_in += w
                mp.community[new_comm].in_weight -= mg.edge_weights[e_idx]
            end
        end
        for e in in_edges(u,g)
            e_idx = edge_index(e,g)
            v = source(e,g)
            v_idx = vertex_index(v,g)
            v_comm = mp.membership[v_idx]
            if in(v, mp.community[old_comm].nodes)
                mp.community[old_comm].out_weight += mg.edge_weights[e_idx]
                w = is_directed(g) ? mg.edge_weights[e_idx] : mg.edge_weights[e_idx]/2
                mp.community[old_comm].inner_weight -= w
                mp.total_weight_in -= w
            else
                mp.community[old_comm].in_weight -= mg.edge_weights[e_idx]
            end
            if in(v, mp.community[new_comm].nodes)
                mp.community[new_comm].out_weight -= mg.edge_weights[e_idx]
                w = is_directed(g) ? mg.edge_weights[e_idx] : mg.edge_weights[e_idx]/2
                mp.community[new_comm].inner_weight += w
                mp.total_weight_in += w
            else
                mp.community[new_comm].in_weight += mg.edge_weights[e_idx]
            end
        end

        # if the old community is empty after remove node u, we remove it
        if isempty(mp.community[old_comm].nodes)
            delete!(mp.community, old_comm)
        end

        # update the membership vector
        mp.membership[u_idx] = new_comm
    end
end

function from_coarser_partition!{V}(partition::ModularityPartition{V}, coarser_partition::ModularityPartition{V})
    for u in vertices(partition.mod_graph.graph)
        u_idx = vertex_index(u, partition.mod_graph.graph)
        # what is the community of the node
        u_comm_level1 = partition.membership[u_idx]

        # In the coarser partition, the node should have the community id
        # so that the community of that node gives the coarser community.
        u_comm_level2 = coarser_partition.membership[u_comm_level1]
        partition.membership[u_idx] = u_comm_level2
    end
    update_partition!(partition)
end

# weight from vertex to community
function weight_to_comm{V}(partition::ModularityPartition{V}, u::V, comm::Int)
    total_weight = 0.0
    u_idx = vertex_index(u, partition.mod_graph.graph)
    u_comm = partition.membership[u_idx]
    if u_comm != comm
        for e in out_edges(u, partition.mod_graph.graph)
            e_idx = edge_index(e, partition.mod_graph.graph)
            v = target(e, partition.mod_graph.graph)
            if in(v, partition.community[comm].nodes)
                total_weight += partition.mod_graph.edge_weights[e_idx]
            end
        end
    end
    total_weight
end

function weight_from_comm{V}(partition::ModularityPartition{V}, u::V, comm::Int)
    total_weight = 0.0
    u_idx = vertex_index(u, partition.mod_graph.graph)
    u_comm = partition.membership[u_idx]
    if u_comm != comm
        for e in in_edges(u, partition.mod_graph.graph)
            e_idx = edge_index(e, partition.mod_graph.graph)
            v = target(e, partition.mod_graph.graph)
            if in(v, partition.community[comm].nodes)
                total_weight += partition.mod_graph.edge_weights[e_idx]
            end
        end
    end
    total_weight
end

function in_strength{V}(mg::ModularityGraph{V}, u::V)
    instrength = 0.0
    for e in in_edges(u,g)
        e_idx = edge_index(e,g)
        instrength += mg.edge_weights[e_idx]
    end
    instrength
end

function out_strength{V}(mg::ModularityGraph{V}, u::V)
    instrength = 0.0
    for e in out_edges(u,g)
        e_idx = edge_index(e,g)
        instrength += mg.edge_weights[e_idx]
    end
    instrength
end

function diff_move{V}(mp::ModularityPartition{V}, u::V, new_comm::Int)
    mg = mp.mod_graph
    g = mg.graph
    u_idx = vertex_index(u, g)
    old_comm = mp.membership[u_idx]
    diff = 0.0
    if new_comm != old_comm
        w_to_old = weight_to_comm(mp, u, old_comm)
        w_from_old = weight_from_comm(mp, u, old_comm)
        w_to_new = weight_to_comm(mp, u, new_comm)
        w_from_new = weight_from_comm(mp, u, new_comm)
        k_out = out_strength(mg, u)
        k_in = in_strength(mg, u)
        self_weight = mg.node_self_weights[u_idx]
        K_out_old = mp.community[old_comm].out_weight
        K_in_old = mp.community[old_comm].in_weight
        K_out_new = mp.community[new_comm].out_weight
        K_in_new = mp.community[new_comm].in_weight
        total_weight = mg.total_weight*(2-Int(is_directed(g)))
        diff_old = (w_to_old - k_out*K_in_old/total_weight) +
            (w_from_old - k_in*K_out_old/total_weight)
        diff_new = (w_to_new + self_weight - k_out*K_in_new/total_weight) +
            (w_from_new + self_weight - k_in*K_out_new/total_weight)

        diff = diff_new - diff_old
    end
    diff
end

"get neighbor communities of node u"
function get_neigh_comms{V}(partition::ModularityPartition{V}, u::V)
    neigh_comms = Set{Int}()
    for v in out_neighbors(u, partition.mod_graph.graph)
        v_idx = vertex_index(v, partition.mod_graph.graph)
        push!(neigh_comms, partition.membership[v_idx])
    end
    neigh_comms
end

function quality{V}(mp::ModularityPartition{V})
    mg = mp.mod_graph
    g = mg.graph
    Q = 0.0
    for grp in values(mp.community)
        w = grp.inner_weight
        w_out = grp.out_weight
        w_in = grp.in_weight
        Q += w - w_out*w_in/((is_directed(g) ? 1.0 : 4.0)*mg.total_weight)
    end
    (2.0 - is_directed(g))*Q
end

function collapse_partition{V}(partition::ModularityPartition{V})
    mg = partition.mod_graph
    g = mg.graph
    num_comm = length(partition.community)
    collapsed_edge_weights = Dict{Int,Float64}[]
    for i=1:num_comm
        push!(collapsed_edge_weights, Dict{Int,Float64}())
    end

    for e in edges(g)
        e_idx = edge_index(e, g)
        u = source(e, g)
        v = target(e, g)
        u_idx = vertex_index(u, g)
        v_idx = vertex_index(v, g)
        u_comm = partition.membership[u_idx]
        v_comm = partition.membership[v_idx]

        # we just skip self loop
        if is_directed(g)
            if u_comm < v_comm
                if haskey(collapsed_edge_weights[u_comm], v_comm)
                    collapsed_edge_weights[u_comm][v_comm] += mg.edge_weights[e_idx]
                else
                    collapsed_edge_weights[u_comm][v_comm] = mg.edge_weights[e_idx]
                end
            end
        else
            if u_comm != v_comm
                if haskey(collapsed_edge_weights[u_comm], v_comm)
                    collapsed_edge_weights[u_comm][v_comm] += mg.edge_weights[e_idx]
                else
                    collapsed_edge_weights[u_comm][v_comm] = mg.edge_weights[e_idx]
                end
            end
        end
    end

    graph = simple_graph(num_comm, is_directed=is_directed(g))
    graph_edge_weights = Float64[]
    graph_node_size = Array(Int, num_comm)
    graph_self_weights = Array(Float64, num_comm)

    for u_comm=1:num_comm
        graph_node_size[u_comm] = partition.community[u_comm].csize
        graph_self_weights[u_comm] = partition.community[u_comm].inner_weight
        for (v_comm, edge_weight) in collapsed_edge_weights[u_comm]
            add_edge!(graph, u_comm, v_comm)
            push!(graph_edge_weights, edge_weight)
        end
    end
    mg = modularity_graph(graph, graph_edge_weights, graph_node_size, graph_self_weights)
    modularity_partition(mg)
end

mp.membership=collect(1:34)
update_partition!(mp)
mp.community
diff_move(mp, 1, 2)
quality(mp)
move_node!(mp, 2, 3)
optimize_partition(mp)
move_nodes(mp)
