using Graphs
import Graphs.graph

type FlowGroup{V}
    nodes::Set{V}
    inner_prob::Float64
    exit_prob::Float64
end

type DiFlowPartition{V} <: AbstractPartition{V}
    flowgraph::DiFlowGraph{V}
    membership::Vector{Int}
    community::Dict{Int,FlowGroup{V}}
    total_exit_prob::Float64
end

# construction
function diflow_partition{V}(fg::DiFlowGraph{V})
    n = num_vertices(fg.graph)
    membership = collect(1:n)
    groups = Array(FlowGroup{V}, n)

    for u in vertices(fg.graph)
        u_idx = vertex_index(u, fg.graph)
        exit_prob = 0.0
        if out_degree(u, fg.graph) > 0
            for e in out_edges(u, fg.graph)
                e_idx = edge_index(e, fg.graph)
                exit_prob += fg.tau*fg.visit_prob[u_idx]*(n-1)/n + (1.0-fg.tau)*fg.trans_prob[e_idx]
            end
        else
            exit_prob += fg.visit_prob[u_idx]*(n-1)/n
        end
        groups[u_idx] = FlowGroup(Set(u), fg.visit_prob[u_idx], exit_prob)
    end

    community = Dict{Int,FlowGroup{V}}(zip(membership, groups))

    total_exit_prob = 0.0
    for group in values(community)
        total_exit_prob += group.exit_prob
    end

    DiFlowPartition(fg, membership, community, total_exit_prob)
end

function diflow_partition{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T}; τ=0.15)
    fg = diflowgraph(g, weights, τ=τ)
    diflow_partition(fg)
end

function diflow_partition{V}(fg::DiFlowGraph{V}, membership::Vector{Int}; τ=0.15)
    maximum(membership) ≤ num_vertices(fg.graph) || error("maximum(membership) must less than num_vertices(g)")
    minimum(membership) > 0 || error("value of membership must be positive integer")

    g = fg.graph
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
                    for w in community[comm_idx].nodes
                        w_idx = vertex_index(w, g)
                        if out_degree(w, g) > 0
                            community[comm_idx] -= fg.tau*fg.visit_prob[w_idx]/n
                        else
                            community[comm_idx] -= fg.visit_prob[w_idx]/n
                        end
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
