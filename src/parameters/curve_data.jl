export DlogData, CurveData

struct DlogData
    e::Int
    window_size::Int
    T1::Vector{Vector{FqFieldElem}}
    T2::Vector{Vector{FqFieldElem}}
    strategy::Vector{Int}
end

# structure for precomputed values
struct CurveData{T <: RingElem}
    A0::T
    A0d::T
    A0dd::T
    a24_0::Proj1{T}
    j0::T
    P2e::Point{T}
    Q2e::Point{T}
    xP2e::Proj1{T}
    xQ2e::Proj1{T}
    xPQ2e::Proj1{T}
    xP2e_short::Proj1{T}
    xQ2e_short::Proj1{T}
    xPQ2e_short::Proj1{T}
    DegreesOddTorsionBases::Vector{Int}
    ExponentsOddTorsionBases::Vector{Int}
    OddTorsionBases::Vector{Vector{Proj1{T}}}
    Matrices_2e::Vector{Matrix{BigInt}}
    Matrix_2ed_inv::Matrix{BigInt}
    Matrices_odd::Vector{Vector{Matrix{Int}}}
    isomorphism_to_A0::Function
    dlog_data_full::DlogData
    dlog_data_chall::DlogData
end
