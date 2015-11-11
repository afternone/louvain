type ModularityGraph{V}
    graph::AbstractGraph{V}
    edge_weights::Vector{Float64}
    node_size::Vector{Int}
    node_self_weights::Vector{Float64}
    total_weight::Float64
    total_size::Int
    correct_self_loops::Bool
    density::Float64
end

using Graphs
function modularity_graph{V,T<:Real}(g::AbstractGraph{V},
                                     edge_weights::Vector{T}=ones(num_edges(g)),
                                     node_size::Vector{Int}=ones(Int,num_vertices(g)),
                                     node_self_weights::Vector{T}=zeros(num_vertices(g)),
                                     correct_self_loops::Bool=false)
    w = sum(edge_weights)
    n_size = sum(node_size)
    normalise = correct_self_loops ? n_size*n_size : n_size*(n_size-1)
    density = is_directed(g) ? w/normalise : 2w/normalise
    ModularityGraph(g, edge_weights, node_size, node_self_weights, w, n_size, correct_self_loops, density)
end
