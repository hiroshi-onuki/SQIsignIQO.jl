
const chi_eval = Dict{Tuple{Int, Int}, Int}(
    (0, 0) => 1,
    (0, 1) => 1,
    (0, 2) => 1,
    (0, 3) => 1,
    (1, 0) => 1,
    (1, 1) => -1,
    (1, 2) => 1,
    (1, 3) => -1,
    (2, 0) => 1,
    (2, 1) => 1,
    (2, 2) => -1,
    (2, 3) => -1,
    (3, 0) => 1,
    (3, 1) => -1,
    (3, 2) => -1,
    (3, 3) => 1
)

const even_indices = [
    [0, 0],
    [0, 1],
    [0, 2],
    [0, 3],
    [1, 0],
    [1, 2],
    [2, 0],
    [2, 1],
    [3, 0],
    [3, 3]
]

const splitting_map = Dict{Tuple{Int, Int}, Vector{Int}}(
    (0, 2) => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0],
    (3, 3) => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
    (0, 3) => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1],
    (2, 1) => [1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, 1, -1, -1],
    (0, 1) => [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0],
    (1, 2) => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0],
    (2, 0) => [1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, -1, -1, 1, 1],
    (3, 0) => [1, 1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1],
    (1, 0) => [1, 1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, -1, 1, -1, 1]
)

function level_22_constants_sqr(tnull::ThetaNullLv2, chi::Integer, i::Integer)
    Ucontant = 0
    for t in 0:3
        Ucontant += chi_eval[(chi, t)] * tnull[t + 1] * tnull[(i ⊻ t) + 1]
    end
    return Ucontant
end

function identify_even_index(null_point::ThetaNullLv2)
    for (chi, i) in even_indices
        U_sqr = level_22_constants_sqr(null_point, chi, i)
        if U_sqr == 0
            return chi, i
        end
    end
    error("Not a product of elliptic_curves")
end

function compute_splitting_matrix(tnull::ThetaNullLv2)
    chi, i = identify_even_index(tnull)
    if chi == 0 && i == 0
        zeta = gen(parent(tnull[1]))    # zeta = sqrt(-1)
        return [1, zeta, 1, zeta, 1, -zeta, -1, zeta, 1, zeta, -1, -zeta, -1, zeta, -1, zeta]
    else
        return splitting_map[(chi, i)]
    end
end

function splitting_isomorphism(tnull::ThetaNullLv2{T}, image_points::Vector{ThetaPtLv2{T}}) where T <: RingElem
    M = compute_splitting_matrix(tnull)
    tnull = apply_base_chagne(tnull, M)
    
    ret = Vector{ThetaPtLv2{T}}(undef, length(image_points))
    for i in 1:length(image_points)
        ret[i] = apply_base_chagne(image_points[i], M)
    end
    return ThetaNullLv2(tnull), ret
end