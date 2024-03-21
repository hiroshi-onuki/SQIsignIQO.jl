

# return alpha = (ai + bj + cij)/N s.t. alpha^2 = -d
function sqrt_in_quaternion(d::Integer)
    while true
        N = rand(p+1:10*p)
        !is_prime(N) && continue
        a = (BigInt(sqrtmod(ZZ(d), ZZ(p))) * N) % p
        M = div(d * N^2 - a^2, p)
        b, c, found = sum_of_two_squares(M)
        if found
            return a, b, c, N
        end
    end
end

function compute_order2(E0::E0Data)
    order2e = BigInt(2)^Level1.ExponentFull
    d = 2

    found = false
    a24 = E0.a24_0
    xP, xQ, xPQ = E0.xP2e, E0.xQ2e, E0.xPQ2e
    alpha = QOrderElem(0)
    I = LeftIdeal(Quaternion_0, Quaternion_0, Quaternion_0, Quaternion_0)
    M = BigInt[1 0; 0 1]
    Ma = BigInt[1 0; 0 1]
    N = 1
    D = 1
    while !found
        a, b, c, N = sqrt_in_quaternion(d)
        alpha = Level1.order_elem_from_standard_basis(0, a, b, c)
        Ma = alpha[1] * [1 0; 0 1] + alpha[2] * E0.Matrices_2e[1] + alpha[3] * E0.Matrices_2e[2] + alpha[4] * E0.Matrices_2e[3]
        @assert alpha * alpha == QOrderElem(-d * N^2)

        I = Level1.LeftIdeal(alpha, N)
        alpha, found = KLPT(I, N)
        alpha = div(alpha, gcd(alpha))
        @assert found
        J = Level1.ideal_transform(I, alpha, N)

        e = Int(log(2, div(Level1.norm(J), ExtraDegree)))
        alpha = involution(alpha)
        a24 = E0.a24_0
        xP, xQ, xPQ = E0.xP2e, E0.xQ2e, E0.xPQ2e
        M = BigInt[1 0; 0 1]
        is_first = true
        extdeg = ExtraDegree
        D = 1
        while e > 0
            ed = min(e, ExponentForIsogeny)
            n_I_d = D * extdeg * BigInt(2)^ed
            I_d = larger_ideal(J, n_I_d)
            a24, xP, xQ, xPQ, M, beta, D, found = short_ideal_to_isogeny(I_d, a24, xP, xQ, xPQ, M, D, ed, E0, is_first, Quaternion_0, 0, 0)
            !found && break
            J = ideal_transform(J, beta, n_I_d)
            alpha = div(alpha * involution(beta), n_I_d)
            e -= ed
            is_first = false
            extdeg = 1
        end
    end

    alpha = involution(alpha)
    Malpha = alpha[1] * [1 0; 0 1] + alpha[2] * E0.Matrices_2e[1] + alpha[3] * E0.Matrices_2e[2] + alpha[4] * E0.Matrices_2e[3]
    M = (M * Malpha * invmod(D, order2e)) .% order2e
    Minv = invmod_2x2(M, order2e)
    Mdual = (Minv * N) .% order2e
    Msqrt2 = (M * Ma * Minv * invmod(N, order2e)) .% order2e

    xPs = xDBLe(xP, a24, ExponentFull - ExponentForTorsion)
    xQs = xDBLe(xQ, a24, ExponentFull - ExponentForTorsion)
    xPQs = xDBLe(xPQ, a24, ExponentFull - ExponentForTorsion)

    A = Montgomery_coeff(a24)
    P = Point(A, xP)
    Q = Point(A, xQ)
    PQ = add(P, -Q, Proj1(A))
    if !(xPQ == Proj1(PQ.X, PQ.Z))
        Q = -Q
    end
    tp_table_P = make_pairing_table(A, P, ExponentFull)
    tp_table_Q = make_pairing_table(A, Q, ExponentFull)
    tp_PQ = Tate_pairing_P0(Q, tp_table_P, Cofactor)
    base_dlog = fq_dlog_power_of_2_opt(tp_PQ, E0.dlog_data_full)
    base_dlog = invmod(base_dlog, order2e)

    return OrderData(d, A, jInvariant_a24(a24), a24, xP, xQ, xPQ, xPs, xQs, xPQs, I, Mdual, N, Msqrt2, tp_table_P, tp_table_Q, base_dlog)
end