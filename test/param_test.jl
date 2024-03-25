using Nemo
using KaniSQIsign

function check_torsion_orders(e::Int, E0)
    A0 = E0.A0
    a24 = E0.a24_0
    P2e = E0.P2e
    Q2e = E0.Q2e

    @test is_infinity(mult(BigInt(2)^e, P2e, Proj1(A0)))
    @test !is_infinity(mult(BigInt(2)^(e-1), P2e, Proj1(A0)))
    @test is_infinity(mult(BigInt(2)^e, Q2e, Proj1(A0)))
    @test !is_infinity(mult(BigInt(2)^(e-1), Q2e, Proj1(A0)))
    for i in 1:length(E0.OddTorsionBases)
        l = E0.DegreesOddTorsionBases[i]
        el = E0.ExponentsOddTorsionBases[i]
        xP, xQ, xPQ = E0.OddTorsionBases[i]
        @test is_infinity(ladder(l^el, xP, a24))
        @test is_infinity(ladder(l^el, xQ, a24))
        @test is_infinity(ladder(l^el, xPQ, a24))
    end
end

# [a]P + [b]Q
function linear_comb(a::Integer, b::Integer, xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, a24::Proj1{T}, order::Integer) where T <: RingElem
    xaP = ladder(a, xP, a24) # x(aP)
    xQaP = ladder3pt(order - a, xQ, xP, xPQ, a24) # x(Q - aP)
    xaPbQ = ladder3pt(b, xaP, xQ, xQaP, a24)
    return xaPbQ
end

function matrix_action(basis::Vector{Proj1{T}}, M::Matrix{S}, a24::Proj1{T}, order::Integer) where T <: RingElem where S <: Integer
    xP, xQ, xPQ = basis
    MPQ1 = linear_comb(M[1, 1], M[2, 1], xP, xQ, xPQ, a24, order)
    MPQ2 = linear_comb(M[1, 2], M[2, 2], xP, xQ, xPQ, a24, order)
    return MPQ1, MPQ2
end

# check x([i]P) = -x(P)
function check_i_action(basis::Vector{Proj1{T}}, Mi::Matrix{S}, a24::Proj1{T}, order::S) where T <: RingElem where S <: Integer
    MPQ1, MPQ2 = matrix_action(basis, Mi, a24, order)
    return MPQ1 == -basis[1] && MPQ2 == -basis[2]
end

# check x([i + j]P) = x([i]P) + x([j]P)
function check_ij_action(p::BigInt, basis::Vector{Proj1{T}}, Mij::Matrix{S}, a24::Proj1{T}, order::S) where T <: RingElem where S <: Integer
    MPQ1, MPQ2 = matrix_action(basis, Mij, a24, order)
    MPQ1 = xDBL(MPQ1, a24) # x([i + j]P)
    MPQ2 = xDBL(MPQ2, a24) # x([i + j]Q)

    FrobP = Proj1(basis[1].X^p, basis[1].Z^p)
    xijP = x_add_sub(-basis[1], FrobP, a24) # x([i + j]P) or x([i - j]P)
    if MPQ1 != xijP
        MPQ1 != xADD(-basis[1], FrobP, xijP) && return false
    end

    FrobQ = Proj1(basis[2].X^p, basis[2].Z^p)
    xijQ = x_add_sub(-basis[2], FrobQ, a24) # x([i + j]Q) or x([i - j]Q)
    if MPQ2 != xijQ
        MPQ2 != xADD(-basis[2], FrobQ, xijQ) && return false
    end
    return true
end

# check x([1 + k]P) = x(P) + x([i][j]P)
function check_1k_action(p::BigInt, basis::Vector{Proj1{T}}, M1k::Matrix{S}, a24::Proj1{T}, order::S) where T <: RingElem where S <: Integer
    MPQ1, MPQ2 = matrix_action(basis, M1k, a24, order)
    MPQ1 = xDBL(MPQ1, a24) # x([1 + k]P)
    MPQ2 = xDBL(MPQ2, a24) # x([1 + k]Q)

    xijP = -Proj1(basis[1].X^p, basis[1].Z^p)
    x1ijP = x_add_sub(basis[1], xijP, a24) # x(P) + x([i][j]P) or x(P) - x([i][j]P)
    if MPQ1 != x1ijP
        MPQ1 != xADD(basis[1], xijP, x1ijP) && return false
    end

    xijQ = -Proj1(basis[2].X^p, basis[2].Z^p)
    x1ijQ = x_add_sub(basis[2], xijQ, a24) # x(Q) + x([i][j]Q) or x(Q) - x([i][j]Q)
    if MPQ2 != x1ijQ
        MPQ2 != xADD(basis[2], xijQ, x1ijQ) && return false
    end
    return true
end

function check_matrices_actions(p::BigInt, e::Int, ed::Int, E0)
    a24 = E0.a24_0

    # check actions on 2^e-torsion
    xP2e, xQ2e, xPQ2e = E0.xP2e, E0.xQ2e, E0.xPQ2e
    @test check_i_action([xP2e, xQ2e, xPQ2e], E0.Matrices_2e[1], a24, BigInt(2)^e)
    @test check_ij_action(p, [xP2e, xQ2e, xPQ2e], E0.Matrices_2e[2], a24, BigInt(2)^e)
    @test check_1k_action(p, [xP2e, xQ2e, xPQ2e], E0.Matrices_2e[3], a24, BigInt(2)^e)

    m11 = rand(1:BigInt(2)^ed)
    m21 = rand(1:BigInt(2)^ed)
    m12 = rand(1:BigInt(2)^ed)
    m22 = rand(1:BigInt(2)^ed)
    a, b, c, d = E0.Matrix_2ed_inv * [m11, m21, m12, m22]
    @test ([a 0; 0 a] + b * E0.Matrices_2e[1] + c * E0.Matrices_2e[2] + d * E0.Matrices_2e[3]) .% BigInt(2)^ed == [m11 m12; m21 m22]

    # check actions on odd-torsion
    for i in 1:length(E0.OddTorsionBases)
        l = E0.DegreesOddTorsionBases[i]
        el = E0.ExponentsOddTorsionBases[i]
        xP, xQ, xPQ = E0.OddTorsionBases[i]
        @test check_i_action([xP, xQ, xPQ], E0.Matrices_odd[i][1], a24, l^el)
        @test check_ij_action(p, [xP, xQ, xPQ], E0.Matrices_odd[i][2], a24, l^el)
        @test check_1k_action(p, [xP, xQ, xPQ], E0.Matrices_odd[i][3], a24, l^el)
    end
end

function param_check(param::Module)
    _, _, global_data = param.make_precomputed_values()
    E0 = global_data.E0

    check_torsion_orders(param.ExponentFull, E0)
    check_matrices_actions(param.p, param.ExponentFull, param.SQISIGN_challenge_length, E0)
end

param_check(KaniSQIsign.Level1)
#param_check(KaniSQIsign.Level3)
#param_check(KaniSQIsign.Level5)