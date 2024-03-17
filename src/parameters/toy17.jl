include("toy17/prime.jl")
include("toy17/klpt_constants.jl")

include("../quaternion/order.jl")
include("../quaternion/cornacchia.jl")
include("../quaternion/ideal.jl")
include("../quaternion/klpt.jl")

include("../elliptic_curves/dlog.jl")

include("../ideal_to_isogeny/ideal_to_isogeny.jl")

include("../sqisign/sqisign.jl")

# Fp2 and values in Fp2
function make_field_curve_torsions()
    _, T = polynomial_ring(GF(p), "T")
    Fp2, Fp2_i = finite_field(T^2 + 1, "i")
    
    A0 = Fp2(0)

    # constatns from precompute/toy17.sage
    P2e = Point(49608*Fp2_i + 124440, 67996*Fp2_i + 35927)
    Q2e = Point(81463*Fp2_i + 6631, 35927*Fp2_i + 63075)
    M_i_2e = BigInt[0 131071; 1 0]
    M_ij_2e = BigInt[5990 24780; 24781 125082]
    M_1k_2e = BigInt[106292 5990; 5990 24781]
    M44inv = BigInt[205 154 154 308; 154 308 307 358; 204 409 409 308; 103 204 204 409]
    xP3 = Proj1(107897*Fp2_i + 119558)
    xQ3 = Proj1(23174*Fp2_i + 11513)
    xPQ3 = Proj1(23174*Fp2_i + 119558)
    M_i_3 = [0 1; 2 0]
    M_ij_3 = [2 0; 2 1]
    M_1k_3 = [0 1; 1 1]
    xP5 = Proj1(97807*Fp2_i + 60197)
    xQ5 = Proj1(52971*Fp2_i + 103834)
    xPQ5 = Proj1(33264*Fp2_i + 70874)
    M_i_5 = [1 2; 4 4]
    M_ij_5 = [2 1; 3 3]
    M_1k_5 = [4 2; 0 2]
    xP17 = Proj1(9763*Fp2_i + 109052)
    xQ17 = Proj1(61020*Fp2_i + 119882)
    xPQ17 = Proj1(104068*Fp2_i + 2991)
    M_i_17 = [5 8; 1 12]
    M_ij_17 = [1 1; 7 16]
    M_1k_17 = [11 14; 0 7]

    a24_0 = A_to_a24(A0)
    xP2e = Proj1(P2e.X, P2e.Z)
    xQ2e = Proj1(Q2e.X, Q2e.Z)
    PQ2e = add(P2e, -Q2e, Proj1(A0))
    xPQ2e = Proj1(PQ2e.X, PQ2e.Z)
    xP2e_short = xDBLe(xP2e, a24_0, ExponentForIsogeny)
    xQ2e_short = xDBLe(xQ2e, a24_0, ExponentForIsogeny)
    xPQ2e_short = xDBLe(xPQ2e, a24_0, ExponentForIsogeny)

    # precomputed values for discrete logarithm
    wp_P2e_Q2e = Weil_pairing_2power(A0, P2e, Q2e, ExponentFull)
    window_size = 3
    fq_dlog_table1, fq_dlog_table2 = make_dlog_table(wp_P2e_Q2e, ExponentForIsogeny, window_size)
    strategy_dlog = compute_strategy(div(ExponentFull, window_size) - 1, window_size, 1)
    dlog_data_full = DlogData(ExponentFull, window_size, fq_dlog_table1, fq_dlog_table2, strategy_dlog)
    base = wp_P2e_Q2e^(BigInt(2)^(ExponentFull - SQISIGN_challenge_length))
    fq_dlog_table1_c, fq_dlog_table2_c = make_dlog_table(base, SQISIGN_challenge_length, window_size)
    strategy_dlog_c = compute_strategy(div(SQISIGN_challenge_length, window_size) - 1, window_size, 1)
    dlog_data_chall = DlogData(SQISIGN_challenge_length, window_size, fq_dlog_table1_c, fq_dlog_table2_c, strategy_dlog_c)

    DegreesOddTorsionBases = Int[3, 5, 17]
    ExponentsOddTorsionBases = Int[1, 1, 1]
    OddTorsionBases = [[xP3, xQ3, xPQ3], [xP5, xQ5, xPQ5], [xP17, xQ17, xPQ17]]

    Matrices_2e = [M_i_2e, M_ij_2e, M_1k_2e]
    Matrices_odd = [[M_i_3, M_ij_3, M_1k_3], [M_i_5, M_ij_5, M_1k_5], [M_i_17, M_ij_17, M_1k_17]]

    # make constants for isomorphism to the curve E_A0
    _, T = polynomial_ring(Fp2, "T")
    As = roots((256 * (T^2 - 3)^3 - 1728 * (T^2 - 4))/T^2)
    A0d = As[1]
    beta = -A0d/3
    gamma = square_root(1 / (1 - 3*beta^2))
    gamma = gamma[1]/gamma[2]
    A0dd = As[2]
    beta_d = -A0dd/3
    gamma_d = square_root(1 / (1 - 3*beta_d^2))
    gamma_d = gamma_d[1]/gamma_d[2]
    function isomorphism_to_A0(A::Proj1{FqFieldElem}, P::Proj1{FqFieldElem})
        if A == Proj1(A0)
            return P
        elseif A == Proj1(A0d)
            return Proj1(gamma*(P.X - beta*P.Z), P.Z)
        elseif A == Proj1(A0dd)
            return Proj1(gamma_d*(P.X - beta_d*P.Z), P.Z)
        else
            throw(ArgumentError("A is not A0d or A0dd"))
        end
    end

    return Fp2, Fp2_i, CurveData(A0, A0d, A0dd, a24_0, jInvariant_A(A0), P2e, Q2e, xP2e, xQ2e, xPQ2e, xP2e_short, xQ2e_short, xPQ2e_short, DegreesOddTorsionBases, ExponentsOddTorsionBases, OddTorsionBases, Matrices_2e, M44inv, Matrices_odd, isomorphism_to_A0, dlog_data_full, dlog_data_chall)
end