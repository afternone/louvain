type ModularityGroup{V}
    nodes::Set{V}
    inner_weight::Float64
    in_weight::Float64
    out_weight::Float64
end

type ModularityPartition{V}
    graph::AbstractGraph{V}
    membership::Vector{Int}
    community::Dict{Int,ModularityGroup{V}}
    total_weight::Float64 # total weight in all community
    total_possible_edges::Int # total possible edges in all community
end

function modularity_partition{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T})
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
        community[u_idx] = ModularityGroup(Set(u), 0.0, in_weight, out_weight)
    end
    ModularityPartition(g, membership, community, 0.0, 0)
end

function modularity_partition{V,T<:Real}(g::AbstractGraph{V}, weights::Vector{T}, membership::Vector{Int})
    community = Dict{Int,ModularityGroup{V}}()
    total_weight = 0.0
    total_possible_edges = 0

    for u in vertices(g)
        u_idx = vertex_index(u, g)
        comm_idx = membership[u_idx]
        if haskey(community, comm_idx)
            push!(community[comm_idx].nodes, u)
            for e in out_edges(u, g)
                e_idx = edge_index(e, g)
                v = target(e, g)
                v_idx = vertex_index(v, g)
                if !in(v, community[comm_idx].nodes)
                    community[comm_idx].out_weight += weights[e_idx]
                else
                    community[comm_idx].inner_weight += weights[e_idx]
                    total_weight += weights[e_idx]
                    community[comm_idx].in_weight -= weights[e_idx]
                end
            end
            for e in in_edges(u, g)
                e_idx = edge_index(e, g)
                v = source(e, g)
                v_idx = vertex_index(v, g)
                if !in(v, community[comm_idx].nodes)
                    community[comm_idx].in_weight += weights[e_idx]
                else
                    community[comm_idx].inner_weight += weights[e_idx]
                    total_weight += weights[e_idx]
                    community[comm_idx].out_weight -= weights[e_idx]
                end
            end
        else
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
            community[u_idx] = ModularityGroup(Set(u), 0.0, in_weight, out_weight)
        end
    end

    ModularityPartition(g, membership, community, 0.0, 0)


modularity_partition(g, ones(78))
