using Graphs

# type of directed flow graph
type DiFlowGraph{V}
    graph::AbstractGraph{V}
    tau::Float64 # damping parameter
    visit_prob::Vector{Float64} # nodes visit probability
    trans_prob::Vector{Float64} # edges transform probability
end

"construct directed DiFlowGraph{V} from directed AbstractGraph{V}"
function diflowgraph{V,T<:Real}(graph::AbstractGraph{V}, weights::Vector{T}=ones(num_edges(graph)); τ=0.15)
    relative_trans_prob = relative_trans_prob_directed(graph, weights)
    visit_prob = visit_prob_directed(graph, relative_trans_prob, τ=τ, ϵ=sqrt(eps()), N=1000)
    trans_prob = trans_prob_directed(graph, relative_trans_prob, visit_prob)
    DiFlowGraph(graph, τ, visit_prob, trans_prob)
end

"calculate relative transform probability of edges in a directed graphs"
function relative_trans_prob_directed{V,T<:Real}(graph::AbstractGraph{V}, weights::Vector{T}=ones(num_edges(graph)))
    @graph_requires graph vertex_list vertex_map adjacency_list edge_map

    is_directed(graph) || error("graph must be directed.")
    length(weights)==num_edges(graph) || error("length(weights) must equal num_edges(graph)")

    relative_trans_prob = zeros(num_edges(graph))

    for u in vertices(graph)
        sumw = 0.0
        for e in out_edges(u, graph)
            sumw += weights[edge_index(e, graph)]
        end
        if sumw > 0.0
            for e in out_edges(u, graph)
                e_idx = edge_index(e, graph)
                relative_trans_prob[e_idx] = weights[e_idx]/sumw
            end
        end
    end

    relative_trans_prob
end

"""
calculate visit probability of each node in directed networks
`tau` is damping parameter,  `epsilon` is convergence precision, `N` is maximum iterations
"""
function visit_prob_directed{V,T<:Real}(graph::AbstractGraph{V}, relative_trans_prob::Vector{T}; τ=0.15, ϵ=10e-8, N=1000)
    @graph_requires graph vertex_list vertex_map adjacency_list edge_map

    is_directed(graph) || error("graph must be directed.")
    length(relative_trans_prob)==num_edges(graph) || error("length(relative_trans_prob) must equal num_edges(graph)")

    n = num_vertices(graph)
    p = fill(1/n, n) # initial visit probability
    p1 = zeros(n)
    i = 0 # initial iterations
    δp = 1.0 # δp > ϵ for into while loop at first time

    while δp > ϵ && i < N
        i += 1
        dp = 0.0
        for u in vertices(graph)
            if out_degree(u, graph) == 0
                dp += (1.0 - τ) * p[vertex_index(u, graph)] / n
            end
        end
        for u in vertices(graph)
            u_idx = vertex_index(u, graph)
            p1[u_idx] = dp + τ/n
            for e in in_edges(u, graph)
                v = source(e, graph)
                v_idx = vertex_index(v, graph)
                e_idx = edge_index(e, graph)
                p1[u_idx] += (1.0-τ) * relative_trans_prob[e_idx] * p[v_idx]
            end
        end
        δp = maxabs(p1-p)
        p[:] = p1[:]
    end

    p
end

"calculate transform probability of edges in directed graphs"
function trans_prob_directed{V,T<:Real}(graph::AbstractGraph{V}, relative_trans_prob::Vector{T}, visit_prob::Vector{T})
    @graph_requires graph vertex_list vertex_map adjacency_list edge_map

    is_directed(graph) || error("graph must be directed.")
    length(visit_prob )== num_vertices(graph) || error("length(visit_prob) must equal num_vertices(graph)")
    length(relative_trans_prob )== num_edges(graph) || error("length(relative_trans_prob) must equal num_edges(graph)")

    trans_prob = zeros(T, num_edges(graph))

    for e in edges(graph)
        e_idx = edge_index(e, graph)
        u = source(e, graph)
        u_idx = vertex_index(u, graph)
        trans_prob[e_idx] = visit_prob[u_idx] * relative_trans_prob[e_idx]
    end

    trans_prob
end
