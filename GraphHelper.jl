"""
Creates a graph with communities as node and links as weights between communities.

The weight of the edges in the new graph is simply the sum of the weight
of the edges between the communities. The size of a node in the new graph
is simply the size of the community in the old graph.
"""
function collapse_graph{V}(partition::FlowPartition{V})
    num_comm = length(partition.community)
    collapsed_trans_prob = fill(Dict{Int,Float64}(), num_comm)

    for e in edges(partition.flowgraph.graph)
        e_idx = edge_index(e, partition.flowgraph.graph)
        u = source(e, partition.flowgraph.graph)
        v = target(e, partition.flowgraph.graph)
        u_idx = vertex_index(u, partition.flowgraph.graph)
        v_idx = vertex_index(v, partition.flowgraph.graph)
        u_comm = partition.membership[u_idx]
        v_comm = partition.membership[v_idx]

        # we just skip self loop
        if haskey(collapsed_trans_prob[u_comm], v_comm)
            collapsed_trans_prob[u_comm][v_comm] += partition.flowgraph.trans_prob[e_idx]
        else
            collapsed_trans_prob[u_comm][v_comm] = partition.flowgraph.trans_prob[e_idx]
        end
    end

    graph = simple_graph(num_comm, is_directed=false)
    graph_trans_prob = Float64[]
    graph_visit_prob = Array(Float64, num_comm)

    for u_comm=1:num_comm
        graph_visit_prob[u_comm] = partition.community[u_comm].inner_prob
        for (v_comm, trans_prob) in collapsed_trans_prob[u_comm]
            if u_comm < v_comm
                add_edge!(graph, u_comm, v_comm)
                push!(graph_trans_prob, trans_prob)
            end
        end
    end

    fg = FlowGraph(graph, graph_visit_prob, graph_trans_prob)
    flow_partition(fg)
end
