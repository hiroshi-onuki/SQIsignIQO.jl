
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

# return the kernel of the l^e-isogeny determined by the kernel matrix Mker
function kernel_gen_power_of_prime(xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, a24::Proj1{T},
    Mker::Matrix{S}, Mimage::Matrix{S}, l::Int, e::Int) where T <: RingElem where S <: Integer
    a, b = kernel_coefficients(Mker, l, e)
    a, b = Mimage * [a, b]
    if a % l != 0
        b = (b * invmod(a, S(l)^e)) % S(l)^e
        b < 0 && (b += S(l)^e)
        return ladder3pt(b, xP, xQ, xPQ, a24)
    else
        a = (a * invmod(b, S(l)^e)) % S(l)^e
        a < 0 && (a += S(l)^e)
        return ladder3pt(a, xQ, xP, xPQ, a24)
    end
end

# I : left ideal of O0 s.t. I = I_1 I_2, n(I_1) is odd and n(I_2) = 2^e for e <= ExponentForIsogeny
# a24 : the coefficient of E := E_0 / E_0[I_1]
# xP, xQ, xPQ : the fixed basis of E[2^ExponentFull] s.t. phi_I_2(P_0, Q_0)^t = M(P, Q)^t
# return the coefficient a24d of E' := E_0 / E_0[I_1 I_2],
# the fixed basis (P', Q') of E'[2^ExponentFull], and M' s.t. phi_J(P_0, Q_0)^t = M'(P', Q')^t
# If is_special is true, then n(I_2) = 2^ExponentForIsogeny * ExtraDegree
function short_ideal_to_isogeny(I::LeftIdeal, a24::Proj1{T}, xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T},
    M::Matrix{BigInt}, e::Int, cdata::CurveData, is_special::Bool=false) where T <: RingElem
    xPd = xDBLe(xP, a24, ExponentFull - e)
    xQd = xDBLe(xQ, a24, ExponentFull - e)
    xPQd = xDBLe(xPQ, a24, ExponentFull - e)

    # 2^e-isogeny corresponding to I_2
    alpha = primitive_element(I)
    M0 = alpha[1]*[1 0; 0 1] + alpha[2]*cdata.Matrices_2e[1] + alpha[3]*cdata.Matrices_2e[2] + alpha[4]*cdata.Matrices_2e[3]
    ker = kernel_gen_power_of_prime(xPd, xQd, xPQd, a24, M0, M, 2, e)
    eval_points = [xP, xQ, xPQ]
    if is_special
        push!(eval_points, ker)
        degs = Int[]
        for i in 1:length(cdata.DegreesOddTorsionBases)
            l = cdata.DegreesOddTorsionBases[i]
            push!(degs, l)
            xPl, xQl, xPQl = cdata.OddTorsionBases[i]
            Ml = alpha[1] * [1 0; 0 1] + alpha[2] * cdata.Matrices_odd[i][1] + alpha[3] * cdata.Matrices_odd[i][2] + alpha[4] * cdata.Matrices_odd[i][3]
            ker_l = kernel_gen_power_of_prime(xPl, xQl, xPQl, a24, Ml, M, l, 1)
            push!(eval_points, ker_l)
        end
        for i in 1:length(cdata.DegreesOddTorsionBasesTwist)
            l = cdata.DegreesOddTorsionBasesTwist[i]
            push!(degs, l)
            xPl, xQl, xPQl = cdata.OddTorsionBasesTwist[i]
            Ml = alpha[1] * [1 0; 0 1] + alpha[2] * cdata.Matrices_odd_twist[i][1] + alpha[3] * cdata.Matrices_odd_twist[i][2] + alpha[4] * cdata.Matrices_odd_twist[i][3]
            ker_l = kernel_gen_power_of_prime(xPl, xQl, xPQl, a24, Ml, M, l, 1)
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
    M_beta_bar = beta_bar[1]*[1 0; 0 1] + beta_bar[2]*cdata.Matrices_2e[1] + beta_bar[3]*cdata.Matrices_2e[2] + beta_bar[4]*cdata.Matrices_2e[3]
    c11, c21, c12, c22 = M * M_beta_bar
    xP2 = linear_comb_2_e(c11, c21, xP2t, xQ2t, xPQ2t, a24d, ExponentFull)
    xQ2 = linear_comb_2_e(c12, c22, xP2t, xQ2t, xPQ2t, a24d, ExponentFull)
    xPQ2 = linear_comb_2_e(c11-c12, c21-c22, xP2t, xQ2t, xPQ2t, a24d, ExponentFull)
    @assert is_infinity(xDBLe(xP2, a24d, ExponentFull - e))
    @assert is_infinity(xDBLe(xQ2, a24d, ExponentFull - e))
    @assert is_infinity(xDBLe(xPQ2, a24d, ExponentFull - e))

    # compute the images of the basis of E_0[2^ExponentFull] under norm(I)(a + bi)
    c11, c21, c12, c22 = [a 0; 0 a] + b * cdata.Matrices_2e[1]
    a24_0 = cdata.a24_0
    xP0 = cdata.xP2e
    xQ0 = cdata.xQ2e
    xPQ0 = cdata.xPQ2e
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
    A1 = cdata.A0
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

    @assert Weil_pairing_2power(A1, P1, Q1, ExponentForTorsion)^BigInt(2)^(ExponentForTorsion - 1) != 1
    @assert Weil_pairing_2power(A2, P2, Q2, ExponentForTorsion)^BigInt(2)^(ExponentForTorsion - 1) != 1
    @assert Weil_pairing_2power(A1, P1, Q1, ExponentForTorsion) * Weil_pairing_2power(A2, P2, Q2, ExponentForTorsion) == 1
    # end of pairing check

    # fixed basis of E'[2^ExponentFull]
    xPd, xQd, xPQd = torsion_basis(a24d, ExponentFull)

    # compute (2,2)-isogenies
    P1P2 = CouplePoint(xP1, xP2)
    Q1Q2 = CouplePoint(xQ1, xQ2)
    PQ1PQ2 = CouplePoint(xPQ1, xPQ2)
    O1 = infinity_point(parent(a24_0.X))
    O1Pd = CouplePoint(O1, xPd)
    O1Qd = CouplePoint(O1, xQd)
    O1PQd = CouplePoint(O1, xPQd)
    Es, images = product_isogeny_sqrt_no_strategy(a24_0, a24d, P1P2, Q1Q2, PQ1PQ2, [O1Pd, O1Qd, O1PQd], ExponentForTorsion)

    # isomorphism to A0
    if Es[1] == Proj1(cdata.A0) || Es[1] == Proj1(cdata.A0d) || Es[1] == Proj1(cdata.A0dd)
        idx = 1
    else
        idx = 2
    end
    xPdd = images[1][idx]
    xQdd = images[2][idx]
    xPQdd = images[3][idx]
    xPdd = cdata.isomorphism_to_A0(Es[idx], xPdd)
    xQdd = cdata.isomorphism_to_A0(Es[idx], xQdd)
    xPQdd = cdata.isomorphism_to_A0(Es[idx], xPQdd)
    @assert is_infinity(xDBLe(xPdd, a24_0, ExponentFull))
    @assert is_infinity(xDBLe(xQdd, a24_0, ExponentFull))
    @assert is_infinity(xDBLe(xPQdd, a24_0, ExponentFull))

    # compute the matrix M' s.t. phi_J(P0, Q0) = (Pd, Qd)M'
    c11, c21, c12, c22 = ec_dlog_power_of_2(xPdd, xQdd, xPQdd, cdata.P2e, cdata.Q2e, cdata.A0, ExponentFull)
    D = BigInt(2)^ExponentForTorsion - a^2 - b^2
    Md = [c22 -c12; -c21 c11] * invmod((c11 * c22 - c12 * c21), BigInt(2)^ExponentFull) * D

    return a24d, xPd, xQd, xPQd, Md, beta, D
end

# isogeny E0 to E0/E0[I], where n(I) = ExtraDegree*2^e
function ideal_to_isogeny_from_O0(I::LeftIdeal, e::Int, cdata::CurveData)
    a24 = cdata.a24_0
    xP = cdata.xP2e
    xQ = cdata.xQ2e
    xPQ = cdata.xPQ2e
    M = BigInt[1 0; 0 1]

    # the first isogeny is the special case
    I_d = larger_ideal(I, ExtraDegree * BigInt(2)^ExponentForIsogeny)
    a24, xP, xQ, xPQ, M, beta, D = short_ideal_to_isogeny(I_d, a24, xP, xQ, xPQ, M, ExponentForIsogeny, cdata, true)
    I = ideal_transform(I, beta, ExtraDegree * BigInt(2)^ExponentForIsogeny)
    e -= ExponentForIsogeny

    while e > 0
        e_d = min(e, ExponentForIsogeny)
        I_d = larger_ideal(I, D*BigInt(2)^e_d)
        a24, xP, xQ, xPQ, M, beta, D_new = short_ideal_to_isogeny(I_d, a24, xP, xQ, xPQ, M, e_d, cdata)
        I = ideal_transform(I, beta, D*BigInt(2)^e_d)
        e -= e_d
        D = D_new
    end

    return a24, xP, xQ, xPQ, M
end