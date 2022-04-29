# This is Asher's version of the footrule

using LinearAlgebra

"""
    RankingComparison(data_ranking, model_ranking)
    RankingComparison(treatment_model, ranking_data, treatments)
A comparison of two rankings.  Rankings are stored as a mapping from treatments to ranks.
"""
struct RankingComparison{T}
    data_ranking::Dict{Set{T}, Int}
    model_ranking::Dict{Set{T}, Int}
    RankingComparison(data_ranking, model_ranking) = begin
        T = promote_type(
            eltype(eltype(data_ranking).parameters[1]),
            eltype(eltype(model_ranking).parameters[1]),
        )
        new{T}(data_ranking, model_ranking)
    end
end
Base.show(io::IO, comp::RankingComparison) = print(io, "RankingComparison(n_data = $(length(comp.data_ranking)), n_model = $(length(comp.model_ranking)))")

function spearman_footrule_normalization(n)
    return sum([abs(2 * i - n - 1) for i in 1:n])
end

function spearman_footrule(σ::AbstractVector{<:Integer})
    n = length(σ)
    F = sum(i -> abs(i - σ[i]), 1:n)
    return F / spearman_footrule_normalization(n)
end

function spearman_footrule(σ1::AbstractVector{<:Integer}, σ2::AbstractVector{<:Integer})
    length(σ1) != length(σ2) && throw("ranking arrays must have the same length")
    return spearman_footrule(σ2[sortperm(σ1)])
end

function spearman_footrule(comp::RankingComparison)
    data_dict = comp.data_ranking
    model_dict = comp.model_ranking
    data_max_rank = maximum(values(data_dict))
    model_max_rank = maximum(values(model_dict))
    combined_txs = collect(union(collect(keys(data_dict)), collect(keys(model_dict))))
    σ1 = [get(data_dict, tx, data_max_rank + 1) for tx in combined_txs]
    σ2 = [get(model_dict, tx, model_max_rank + 1) for tx in combined_txs]
    return spearman_footrule(σ1, σ2)
end

function weighted_spearman_footrule_normalization(n, w)
    return sum([w[i] * abs(2 * i - n - 1) for i in 1:n])
end

function weighted_spearman_footrule(
    σ::AbstractVector{<:Integer},
    position_weights::AbstractVector{<:Real},
    element_weights::AbstractVector{<:Real},
    similarity_weights::AbstractMatrix{<:Real}
    )
    n = length(σ)
    w = Float64.(element_weights)
    D = similarity_weights
    if !all(==(1), position_weights)
        p_i = prepend!(1 .+ cumsum(position_weights), [one(eltype(position_weights))])
        for i in 1:n
            s = σ[i]
            if i != s
                w[i] *= (p_i[i] - p_i[s]) / (i - s)
            end
        end
    end
    F = 0
    for i in 1:n
        u = sum([w[j] * D[i, j] for j in 1:i])
        v = sum([w[j] * D[i, j] for j in filter(j -> σ[j] <= σ[i], 1:n)])
        F += w[i] * abs(u - v)
    end
    return F / weighted_spearman_footrule_normalization(n, w)
end

function symmetrized_weighted_spearman_footrule(
    σ::AbstractVector{<:Integer},
    position_weights::AbstractVector{<:Real},
    element_weights::AbstractVector{<:Real},
    similarity_weights::AbstractMatrix{<:Real}
    )
    σ_inv = sortperm(σ)
    F = weighted_spearman_footrule(σ, position_weights, element_weights, similarity_weights)
    F_inv = weighted_spearman_footrule(σ_inv, position_weights, element_weights, similarity_weights)
    return 0.5 * (F + F_inv)
end

function generalized_spearman_footrule(
    σ::AbstractVector{<:Integer},
    position_weights::AbstractVector{<:Real},
    element_weights::AbstractVector{<:Real},
    similarity_weights::AbstractMatrix{<:Real},
    )
    # size checks
    length(σ) != (length(position_weights) + 1) && throw("position weights must have length one less than that of the ranking arrays")
    length(σ) != length(element_weights) && throw("element weights must have the same length as that of the ranking arrays")
    size(similarity_weights, 1) != size(similarity_weights, 2) && throw("similarity weights must be square matrix")
    length(σ) != size(similarity_weights, 1) && throw("similarity weights matrix must have the same length as that of the ranking arrays")
    # value checks
    !all(0 .<= position_weights .<= 1) && throw("must have all position weight values between zero and one")
    !all(0 .<= element_weights .<= 1) && throw("must have all element weight values between zero and one")
    !all(0 .<= similarity_weights .<= 1) && throw("must have all similarity weight values between zero and one")
    !issymmetric(similarity_weights) && throw("similarity matrix must be symmetric")
    return symmetrized_weighted_spearman_footrule(σ, position_weights, element_weights,similarity_weights)
end

"""
    generalized_spearman_footrule(σ1, σ2; position_weights, element_weights, similarity_weights)
    generalized_spearman_footrule(σ; position_weights, element_weights, similarity_weights)
Calculate the generalized Spearman's footrule distance between rankings σ1, σ2.
If only one ranking array is provided, it is assumed that the other ranking array is the identity permutation.
Position weight refer to the cost of swapping the item at index i with the item at index i - 1, and should be a vector of length n - 1 (starting at i = 2).
Element weight refers to the cost of moving item i, and should be a vector of length n.
Similarity weight refers to the cost of swapping item i with item j, and should be a matrix of size (n, n).
Default weights are chosen to recover the original Spearman's Footrule.
Ref. [Kumar & Vassilvitskii 2010](https://dl.acm.org/doi/abs/10.1145/1772690.1772749).
"""
function generalized_spearman_footrule(
    σ::AbstractVector{<:Integer};
    position_weights::Union{AbstractVector{<:Real}, Missing} = missing,
    element_weights::Union{AbstractVector{<:Real}, Missing} = missing,
    similarity_weights::Union{AbstractMatrix{<:Real}, Missing} = missing,
    )
    n = length(σ)
    position_weights = if ismissing(position_weights)
        [1 for i in 2:n]
    else
        position_weights
    end
    element_weights = if ismissing(element_weights)
        [1 for i in 1:n]
    else
        element_weights
    end
    similarity_weights = if ismissing(similarity_weights)
        [1 for i in 1:n, j in 1:n]
    else
        similarity_weights
    end
    return generalized_spearman_footrule(σ, position_weights, element_weights, similarity_weights)
end

function generalized_spearman_footrule(
    σ1::AbstractVector{<:Integer},
    σ2::AbstractVector{<:Integer};
    position_weights::Union{AbstractVector{<:Real}, Missing} = missing,
    element_weights::Union{AbstractVector{<:Real}, Missing} = missing,
    similarity_weights::Union{AbstractMatrix{<:Real}, Missing} = missing,
    )
    length(σ1) != length(σ2) && throw("length of ranking arrays must match")
    n = length(σ1)
    p = sortperm(σ1)
    position_weights = if ismissing(position_weights)
        [1 for i in 2:n]
    else
        position_weights
    end
    element_weights = if ismissing(element_weights)
        [1 for i in 1:n]
    else
        element_weights[p]
    end
    similarity_weights = if ismissing(similarity_weights)
        [1 for i in 1:n, j in 1:n]
    else
        similarity_weights[p, p]
    end
    return generalized_spearman_footrule(σ2[p], position_weights, element_weights, similarity_weights)
end

function generalized_spearman_footrule(
    comp::RankingComparison;
    position_weights::Union{AbstractVector{<:Real}, Missing} = missing,
    element_weights::Union{AbstractVector{<:Real}, Missing} = missing,
    similarity_weights::Union{AbstractMatrix{<:Real}, Missing} = missing,
    )
    data_dict = comp.data_ranking
    model_dict = comp.model_ranking
    data_max_rank = maximum(values(data_dict))
    model_max_rank = maximum(values(model_dict))
    combined_txs = collect(union(collect(keys(data_dict)), collect(keys(model_dict))))
    σ1 = [get(data_dict, tx, data_max_rank + 1) for tx in combined_txs]
    σ2 = [get(model_dict, tx, model_max_rank + 1) for tx in combined_txs]
    return generalized_spearman_footrule(σ1, σ2; position_weights, element_weights, similarity_weights)
end
