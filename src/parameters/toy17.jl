const p = BigInt(2)^17 - 1
const ExponentFull = 17
const ExponentForIsogeny = 6
const ExponentForTorsion = 11
const Cofactor = 1
const ExtraDegree = 3 * 5 * 17

const KLPT_repres_num_gamma_trial = 16384
const KLPT_equiv_bound_coeff = 6
const KLPT_equiv_num_iter = 28561
const SmallPrimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

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
    wp_P2e_Q2e::T
    DegreesOddTorsionBases::Vector{Int}
    DegreesOddTorsionBasesTwist::Vector{Int}
    OddTorsionBases::Vector{Vector{Proj1{T}}}
    OddTorsionBasesTwist::Vector{Vector{Proj1{T}}}
    Matrices_2e::Vector{Matrix{BigInt}}
    Matrices_odd::Vector{Vector{Matrix{Int}}}
    Matrices_odd_twist::Vector{Vector{Matrix{Int}}}
    isomorphism_to_A0::Function
end
