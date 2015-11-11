type MGraph{V}
    graph::AbstractGraph{V}
    node_sizes::Vector{Int}
    edge_weights::Vector{Float64}
    strength_in::Vector{Float64}
    strength_out::Vector{Float64}
    node_self_weights::Vector{Float64}
    total_weight::Float64
    total_size::Int
    correct_self_loops::Bool
    density::Float64
end

type MGroup{V}
    nodes::Set{V}
    weight_inner::Float64
    weight_in::Float64
    weight_out::Float64
end

type MPartition{V} <: AbstractPartition{V}
    mgraph::MGraph{V}
    membership::Vector{Int}
    community::Dict{Int,MGroup{V}}
    total_weight_in_all_comms::Float64
    total_possible_edges_in_all_comms::Float64
end
