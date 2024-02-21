
# return (a, b) s.t. M(a, b)^t = 0 mod l^e
function kernel_coefficients(M::Matrix{T}, l::Int, e::Int) where T <: Integer
    N = T(l)^e
    if M[1, 1] % l != 0
        return (M[1, 2] *  invmod(-M[1, 1], N)) % N, 1
    elseif M[2, 1] % l != 0
        return (M[2, 2] *  invmod(-M[2, 1], N)) % N, 1
    elseif M[1, 2] % l != 0
        return 1, (M[1, 1] *  invmod(-M[1, 2], N)) % N
    else
        return 1, (M[2, 1] *  invmod(-M[2, 2], N)) % N
    end
end

# return the kernel of the l^e-isogeny determined by the kernel matrix M
function kernel_gen_power_of_prime(xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, a24::Proj1{T},
    M::Matrix{S}, l::Int, e::Int) where T <: RingElem where S <: Integer
    a, b = kernel_coefficients(M, l, e)
    if a == 1
        b < 0 && (b += S(l)^e)
        return ladder3pt(b, xP, xQ, xPQ, a24)
    else
        a < 0 && (a += S(l)^e)
        return ladder3pt(a, xQ, xP, xPQ, a24)
    end
end

# I : left ideal of O0 s.t. I = I_1 I_2, n(I_1) is odd and n(I_2) = 2^ExponentForIsogeny
# a24 : the coefficient of E := E_0 / E_0[I_1]
# xP, xQ, xPQ : basis of E[2^ExponentFull] s.t. (P, Q)^t = M(P0, Q0)^t
# return the coefficient a24d of E' := E_0 / E_0[I_1 I_2],
# the fixed basis (P', Q') of E'[2^ExponentFull], and M' s.t. (P', Q')^t = M'(P0, Q0)^t
# If is_special is true, then n(I_2) = 2^ExponentForIsogeny * (small odd number)
function short_ideal_to_isogeny(I::LeftIdeal, a24::Proj1{T}, xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T},
    M::Matrix{BigInt}, tdata::TorsionData, is_special::Bool=false) where T <: RingElem
    e = ExponentForIsogeny
    xPd = xDBLe(xP, a24, ExponentFull - e)
    xQd = xDBLe(xQ, a24, ExponentFull - e)
    xPQd = xDBLe(xPQ, a24, ExponentFull - e)

    # 2^e-isogeny corresponding to I_2
    alpha = primitive_element(I)
    M0 = alpha[1]*[1 0; 0 1] + alpha[2]*tdata.Matrices_2e[1] + alpha[3]*tdata.Matrices_2e[2] + alpha[4]*tdata.Matrices_2e[3]
    ker = kernel_gen_power_of_prime(xPd, xQd, xPQd, a24, M*M0, 2, e)
    eval_points = [xP, xQ, xPQ]
    if is_special
        push!(eval_points, ker)
        degs = Int[]
        for i in 1:length(tdata.DegreesOddTorsionBases)
            l = tdata.DegreesOddTorsionBases[i]
            push!(degs, l)
            xPl, xQl, xPQl = tdata.OddTorsionBases[i]
            Ml = alpha[1] * [1 0; 0 1] + alpha[2] * tdata.Matrices_odd[i][1] + alpha[3] * tdata.Matrices_odd[i][2] + alpha[4] * tdata.Matrices_odd[i][3]
            ker_l = kernel_gen_power_of_prime(xPl, xQl, xPQl, a24, Ml, l, 1)
            push!(eval_points, ker_l)
        end
        for i in 1:length(tdata.DegreesOddTorsionBasesTwist)
            l = tdata.DegreesOddTorsionBasesTwist[i]
            push!(degs, l)
            xPl, xQl, xPQl = tdata.OddTorsionBasesTwist[i]
            Ml = alpha[1] * [1 0; 0 1] + alpha[2] * tdata.Matrices_odd_twist[i][1] + alpha[3] * tdata.Matrices_odd_twist[i][2] + alpha[4] * tdata.Matrices_odd_twist[i][3]
            ker_l = kernel_gen_power_of_prime(xPl, xQl, xPQl, a24, Ml, l, 1)
            push!(eval_points, ker_l)
        end
        while length(degs) > 0
            l = pop!(degs)
            ker = pop!(eval_points)
            a24, eval_points = odd_isogeny(a24, ker, l, eval_points)
        end
        ker = pop!(eval_points)
    end
    a24d, images = two_e_iso(a24, ker, e, eval_points)

    # compute beta in I s.t. J := I*\bar{beta}/n(I) has norm 2^ExpTor - a^2 - b^2
    beta, a, b, found = two_e_good_element(I, ExponentForTorsion)
    !found && throw(ArgumentError("No good element found"))

    # compute the images of the basis of E_0[2^ExponentFull] under the isogeny corresponding to J
    xP2t, xQ2t, xPQ2t = images
    beta_bar = involution(beta) # beta_bar = J \bar{I}
    M_beta_bar = beta_bar[1]*[1 0; 0 1] + beta_bar[2]*tdata.Matrices_2e[1] + beta_bar[3]*tdata.Matrices_2e[2] + beta_bar[4]*tdata.Matrices_2e[3]
    c11, c21, c12, c22 = M_beta_bar * M
    xP2 = linear_comb_2_e(c11, c21, xP2t, xQ2t, xPQ2t, a24d, ExponentFull)
    xQ2 = linear_comb_2_e(c12, c22, xP2t, xQ2t, xPQ2t, a24d, ExponentFull)
    xPQ2 = linear_comb_2_e(c11-c12, c21-c22, xP2t, xQ2t, xPQ2t, a24d, ExponentFull)
    @assert is_infinity(xDBLe(xP2, a24d, ExponentFull - e))
    @assert is_infinity(xDBLe(xQ2, a24d, ExponentFull - e))
    @assert is_infinity(xDBLe(xPQ2, a24d, ExponentFull - e))

    # compute the images of the basis of E_0[2^ExponentFull] under norm(I)(a + bi)
    c11, c21, c12, c22 = [a 0; 0 a] + b * tdata.Matrices_2e[1]
    a24_0 = tdata.a24_0
    xP0 = tdata.xP2e
    xQ0 = tdata.xQ2e
    xPQ0 = tdata.xPQ2e
    xP0 = xDBLe(xP0, a24_0, e)
    xQ0 = xDBLe(xQ0, a24_0, e)
    xPQ0 = xDBLe(xPQ0, a24_0, e)
    N = (norm(I) >> e)
    xP0 = ladder(N, xP0, a24_0)
    xQ0 = ladder(N, xQ0, a24_0)
    xPQ0 = ladder(N, xPQ0, a24_0)
    xP1 = linear_comb_2_e(c11, c21, xP0, xQ0, xPQ0, a24_0, ExponentForTorsion)
    xQ1 = linear_comb_2_e(c12, c22, xP0, xQ0, xPQ0, a24_0, ExponentForTorsion)
    xPQ1 = linear_comb_2_e(c11-c12, c21-c22, xP0, xQ0, xPQ0, a24_0, ExponentForTorsion)

    # pairing check
    A1 = tdata.A0
    P1 = Point(A1, xP1)
    Q1 = Point(A1, xQ1)
    PQ1 = add(P1, -Q1, Proj1(A1))
    if xPQ1 != Proj1(PQ1.X, PQ1.Z)
        Q1 = -Q1
    end
    PQ1 = add(P1, -Q1, Proj1(A1))
    @assert xPQ1 == Proj1(PQ1.X, PQ1.Z)

    A2 = Montgomery_coeff(a24d)
    P2 = Point(A2, xP2)
    Q2 = Point(A2, xQ2)
    PQ2 = add(P2, -Q2, Proj1(A2))
    if xPQ2 != Proj1(PQ2.X, PQ2.Z)
        Q2 = -Q2
    end
    PQ2 = add(P2, -Q2, Proj1(A2))
    @assert xPQ2 == Proj1(PQ2.X, PQ2.Z)

    @assert Weil_pairing_2power(A1, P1, Q1, ExponentForTorsion) * Weil_pairing_2power(A2, P2, Q2, ExponentForTorsion) == 1

    # compute (2,2)-isogenies
    P1P2 = CouplePoint(xP1, xP2)
    Q1Q2 = CouplePoint(xQ1, xQ2)
    PQ1PQ2 = CouplePoint(xPQ1, xPQ2)
    Es, images = product_isogeny_sqrt_no_strategy(a24_0, a24d, P1P2, Q1Q2, PQ1PQ2, CouplePoint{FqFieldElem}[], ExponentForTorsion)

    @assert jInvariant_A(Es[1]) == jInvariant_a24(a24_0) || jInvariant_A(Es[2]) == jInvariant_a24(a24_0)

    return Es
end 