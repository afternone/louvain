using Graphs

type MGroup{V}
    nodes::Set{V}
    csize::Int
    weight_inner::Float64
    weight_in::Float64
    weight_out::Float64
end

type MPartition{V} <: AbstractPartition{V}
    mgraph::MGraph{V}
    membership::Vector{Int}
    community::Dict{Int,MGroup{V}}
    total_weight_in_all_comms::Float64
    total_possible_edges_in_all_comms::Int
end

# require interface
graph{V}(mp::MPartition{V}) = mp.mgraph.graph
membership{V}(mp::MPartition{V}) = mp.membership
membership{V}(mp::MPartition{V}, u::V) = mp.membership[vertex_index(u, mp.mgraph.graph)]
community{V}(mp::MPartition{V}) = keys(mp.community)
is_right_direction{V}(mp::MPartition{V}, quality, new_quality) = new_quality > quality
still_running{V}(mp::MPartition{V}, once_diff, ϵ) = once_diff > ϵ

function mpartition{V,T<:Real}(g::AbstractGraph{V}, edge_weights::Vector{T}=ones(num_edges(g)),
                               node_sizes::Vector{Int}=ones(Int,num_vertices(g)), correct_self_loops::Bool=false)
    mg = mgraph(g, edge_weights, node_sizes, correct_self_loops)
    mpartition(mg)
end

function mpartition{V}(mg::MGraph{V})
    n = num_vertices(mg.graph)
    mp = MPartition{V}(mg, collect(1:n), Dict{Int,MGroup{V}}(), 0.0, 0)
    update_partition!(mp)
    mp
end

function update_partition1!{V}(partition::MPartition{V})
    mg = partition.mgraph
    g = mg.graph
    membership = partition.membership

    empty!(partition.community)
    partition.total_weight_in_all_comms = 0.0

    for u in vertices(g)
        u_idx = vertex_index(u,g)
        comm_idx = membership[u_idx]
        if haskey(partition.community, comm_idx)
            push!(partition.community[comm_idx].nodes, u)
            partition.community[comm_idx].csize += mg.node_sizes[u_idx]
            for e in out_edges(u,g)
                e_idx = edge_index(e,g)
                v = target(e,g)
                v_idx = vertex_index(v,g)
                v_comm = membership[v_idx]
                if !in(v, partition.community[comm_idx].nodes)
                    partition.community[comm_idx].weight_out += mg.edge_weights[e_idx]
                else
                    partition.community[comm_idx].weight_inner += is_directed(g) ? mg.edge_weights[e_idx] : mg.edge_weights[e_idx]/2
                    partition.total_weight_in_all_comms += is_directed(g) ? mg.edge_weights[e_idx] : mg.edge_weights[e_idx]/2
                    partition.community[comm_idx].weight_in -= mg.edge_weights[e_idx]
                end
            end
            for e in in_edges(u,g)
                e_idx = edge_index(e,g)
                v = source(e,g)
                v_idx = vertex_index(v,g)
                v_comm = membership[v_idx]
                if in(v, partition.community[comm_idx].nodes)
                    partition.community[comm_idx].weight_out -= mg.edge_weights[e_idx]
                    partition.community[comm_idx].weight_inner += is_directed(g) ? mg.edge_weights[e_idx] : mg.edge_weights[e_idx]/2
                    partition.total_weight_in_all_comms += is_directed(g) ? mg.edge_weights[e_idx] : mg.edge_weights[e_idx]/2
                else
                    partition.community[comm_idx].weight_in += mg.edge_weights[e_idx]
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
            partition.community[comm_idx] = MGroup(Set(u), mg.node_sizes[u_idx], mg.node_self_weights[u_idx], in_weight, out_weight)
        end
    end
end

function update_partition!{V}(mp::MPartition{V})
    mg = mp.mgraph
    g = mg.graph
    n = num_vertices(g)
    nb_comms = maximum(mp.membership)
    # reset partition
    node_sets = Array(Set{V}, nb_comms)
    for i=1:nb_comms
        node_sets[i] = Set{V}()
    end

    w_inner = zeros(nb_comms)
    w_in = zeros(nb_comms)
    w_out = zeros(nb_comms)
    csize = zeros(Int, nb_comms)

    mp.total_weight_in_all_comms = 0.0
    for u in vertices(g)
        u_idx = vertex_index(u,g)
        u_comm = mp.membership[u_idx]
        # add this node to the communtiy sets
        push!(node_sets[u_comm], u)
        # update the communtiy size
        csize[u_comm] += mg.node_sizes[u_idx]

        # loop over all out edges
        for e in out_edges(u,g)
            e_idx = edge_index(e,g)
            v = target(e,g)
            v_idx = vertex_index(v,g)
            v_comm = mp.membership[v_idx]
            # get the weight of the edge
            w = mg.edge_weights[e_idx]
            # Add weight to the outgoing weight of community of u
            w_out[u_comm] += w
            # Add weight to the incoming weight of community of v
            w_in[v_comm] += w
            # if it is an edge within a community
            if u_comm == v_comm
                if !is_directed(g)
                    w /= 2
                end
                w_inner[u_comm] += w
                mp.total_weight_in_all_comms += w
            end
        end

        # loop over all in edges
        for e in in_edges(u,g)
            e_idx = edge_index(e,g)
            v = source(e,g)
            v_idx = vertex_index(v,g)
            v_comm = mp.membership[v_idx]
            # get the weight of the edge
            w = mg.edge_weights[e_idx]
            # Add weight to the outgoing weight of community of u
            w_in[u_comm] += w
            # Add weight to the incoming weight of community of v
            w_out[v_comm] += w
            # if it is an edge within a community
            if u_comm == v_comm
                if !is_directed(g)
                    w /= 2
                end
                w_inner[u_comm] += w
                mp.total_weight_in_all_comms += w
            end
        end
    end
    mp.total_possible_edges_in_all_comms = 0
    for c=1:nb_comms
        n_c = csize[c]
        possible_edges = 0
        if mg.correct_self_loops
            possible_edges = round(Int, n_c*n_c/(2.0 - Int(is_directed(g))))
        else
            possible_edges = round(Int, n_c*(n_c-1)/(2.0 - Int(is_directed(g))))
        end
        mp.total_possible_edges_in_all_comms += possible_edges
    end

    empty!(mp.community)
    for i=1:nb_comms
        if !isempty(node_sets[i])
            mp.community[i] = MGroup{V}(node_sets[i], csize[i], w_inner[i], w_in[i], w_out[i])
        end
    end
end

"Renumber the communities so that they are numbered 0,...,q-1 where q is the number of communities."
function renumber_communities!{V}(mp::MPartition{V})
    csizes = Int[length(mp.community[i].nodes) for i in keys(mp.community)]
    perm_idx = sortperm(csizes, rev=true)
    mp.community = Dict{Int,MGroup{V}}(zip(enumerate(collect(values(mp.community))[perm_idx])))
    for (i, group) in mp.community
        for u in group.nodes
            u_idx = vertex_index(u, mp.mgraph.graph)
            mp.membership[u_idx] = i
        end
    end
end

"get neighbor communities of node u"
function get_neigh_comms{V}(mp::MPartition{V}, u::V)
    neigh_comms = Set{Int}()
    for v in out_neighbors(u, mp.mgraph.graph)
        v_idx = vertex_index(v, mp.mgraph.graph)
        push!(neigh_comms, mp.membership[v_idx])
    end
    neigh_comms
end

function from_coarser_partition!{V}(partition::MPartition{V}, coarser_partition::MPartition{V})
    for u in vertices(partition.mgraph.graph)
        u_idx = vertex_index(u, partition.mgraph.graph)
        # what is the community of the node
        u_comm_level1 = partition.membership[u_idx]

        # In the coarser partition, the node should have the community id
        # so that the community of that node gives the coarser community.
        u_comm_level2 = coarser_partition.membership[u_comm_level1]
        partition.membership[u_idx] = u_comm_level2
    end
    update_partition!(partition)
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

function move_node!{V}(mp::MPartition{V}, u::V, new_comm::Int)
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

    for e in out_edges(u, g)
        e_idx = edge_index(e,g)
        v = target(e,g)
        v_idx = vertex_index(v,g)
        v_comm = mp.membership[v_idx]
        w = mg.edge_weights[e_idx]
        mp.community[old_comm].weight_out -= w
        mp.community[new_comm].weight_out += w
        int_weight = w/(is_directed(g) ? 1.0 : 2.0)/( u == v ? 2.0 : 1.0)
        if old_comm == v_comm
            # remove the internal weight
            mp.community[old_comm].weight_inner -= int_weight
            mp.total_weight_in_all_comms -= int_weight
        end
        if new_comm == v_comm || u == v
            # add the internal weight
            mp.community[new_comm].weight_inner += int_weight
            mp.total_weight_in_all_comms += int_weight
        end
    end

    for e in in_edges(u, g)
        e_idx = edge_index(e,g)
        v = source(e,g)
        v_idx = vertex_index(v,g)
        v_comm = mp.membership[v_idx]
        w = mg.edge_weights[e_idx]
        mp.community[old_comm].weight_in -= w
        mp.community[new_comm].weight_in += w
        int_weight = w/(is_directed(g) ? 1.0 : 2.0)/( u == v ? 2.0 : 1.0)
        if old_comm == v_comm
            # remove the internal weight
            mp.community[old_comm].weight_inner -= int_weight
            mp.total_weight_in_all_comms -= int_weight
        end
        if new_comm == v_comm || u == v
            # add the internal weight
            mp.community[new_comm].weight_inner += int_weight
            mp.total_weight_in_all_comms += int_weight
        end
    end

    # if the old community is empty after remove node u, we remove it
    if isempty(mp.community[old_comm].nodes)
        delete!(mp.community, old_comm)
    end

    # update the membership vector
    mp.membership[u_idx] = new_comm
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

function collapse_partition{V}(mp::MPartition{V})
    mg = mp.mgraph
    g = mg.graph
    n = num_vertices(g)
    m = num_edges(g)
    nb_comm = length(mp.community)

    collapsed_edge_weights = Array(Dict{Int,Float64}, nb_comm)
    for i=1:nb_comm
        collapsed_edge_weights[i] = Dict{Int,Float64}()
    end

    for e in edges(g)
        e_idx = edge_index(e,g)
        w = mg.edge_weights[e_idx]
        u = source(e,g)
        v = target(e,g)
        u_idx = vertex_index(u,g)
        v_idx = vertex_index(v,g)
        u_comm = mp.membership[u_idx]
        v_comm = mp.membership[v_idx]
        if haskey(collapsed_edge_weights[u_comm], v_comm)
            collapsed_edge_weights[u_comm][v_comm] += w
        else
            collapsed_edge_weights[u_comm][v_comm] = w
        end
    end

    # Now create vector for edges, first determined the number of edges
    n_collapsed = nb_comm
    m_collapsed = 0

    for itr in collapsed_edge_weights
        m_collapsed += length(itr)
    end

    collapsed_graph = simple_graph(n_collapsed, is_directed=is_directed(g))

    collapsed_weights = Float64[]
    total_collapsed_weight = 0.0

    for u=1:n_collapsed
        for (v,w) in collapsed_edge_weights[u]
            add_edge!(collapsed_graph, u, v)
            push!(collapsed_weights, w)
            total_collapsed_weight += w
        end
    end

    if abs(total_collapsed_weight-mg.total_weight) > 1.0e-6
        error("Total collapsed weight is not equal to original weight.")
    end

    if num_vertices(collapsed_graph) != nb_comm
        error("Something went wrong with collapsing the graph.")
    end

    # calculate new node sizes
    csizes = zeros(Int, n_collapsed)
    for (i,c) in mp.community
        csizes[i] = c.csize
    end

    collapsed_mg = mgraph(collapsed_graph, collapsed_weights, csizes, mg.correct_self_loops)
    mpartition(collapsed_mg)
end

function diff_move{V}(mp::MPartition{V}, u::V, new_comm::Int)
    mg = mp.mgraph
    g = mg.graph
    u_idx = vertex_index(u,g)
    old_comm = mp.membership[u_idx]
    diff = 0.0
    if new_comm != old_comm
        w_to_old = weight_to_community(mg, u, old_comm, mp.membership)
        w_from_old = weight_from_community(mg, u, old_comm, mp.membership)
        w_to_new = weight_to_community(mg, u, new_comm, mp.membership)
        w_from_new = weight_from_community(mg, u, new_comm, mp.membership)
        k_out = mg.strength_out[u_idx]
        k_in = mg.strength_in[u_idx]
        self_weight = mg.node_self_weights[u_idx]
        K_out_old = mp.community[old_comm].weight_out
        K_in_old = mp.community[old_comm].weight_in
        K_out_new = mp.community[new_comm].weight_out + k_out
        K_in_new = mp.community[new_comm].weight_in + k_in
        total_weight = mg.total_weight*(2.0 - Float64(is_directed(g)))
        diff_old = (w_to_old - k_out*K_in_old/total_weight) + (w_from_old - k_in*K_out_old/total_weight)
        diff_new = (w_to_new + self_weight - k_out*K_in_new/total_weight) +
            (w_from_new + self_weight - k_in*K_out_new/total_weight)
        diff = diff_new - diff_old
    end
    diff
end

function quality{V}(mp::MPartition{V})
    mg = mp.mgraph
    g = mg.graph
    Q = 0.0
    for c in values(mp.community)
        w = c.weight_inner
        w_out = c.weight_out
        w_in = c.weight_in
        Q += w - w_out*w_in/((is_directed(g) ? 1.0 : 4.0)*mg.total_weight)
    end
    (2.0 - Float64(is_directed(g)))*Q
end

using GraphPlot
out_degree(11,g)
g = graphfamous("karate")
mp = mpartition(g)
mp.community
mp.membership = fill(1,34)
update_partition1!(mp)
mp.community
diff_move(mp, 1, 2)
move_node!(mp, 1, 2)
quality(mp)
optimize_partition(mp)
mp.community
move_nodes(mp)
from_coarser_partition!(mp, mp)
renumber_communities!(mp)
mp1 = collapse_partition(mp)
move_nodes(mp)
edges(mp1.mgraph.graph)
mp1.mgraph.edge_weights
mp2 = collapse_partition(mp1)
move_nodes(mp2)
from_coarser_partition!(mp1, mp2)
mp.community
mp1.community
quality(mp)
optimize_partition(mp)
mp.community
g = graphfamous("karate")
mp = mpartition(g)
82511/2/num_edges(g)
quality(mp2)
move_nodes(mp)
mp2 = collapse_partition(mp1)
move_nodes(mp2)
mp2.community
mp1.community
mp1.mgraph.edge_weights


