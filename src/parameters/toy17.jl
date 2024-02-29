include("toy17/prime.jl")
include("toy17/klpt_constants.jl")

include("../quoternion/order.jl")
include("../quoternion/cornacchia.jl")
include("../quoternion/ideal.jl")
include("../quoternion/klpt.jl")

include("../ideal_to_isogeny/ideal_to_isogeny.jl")

include("../sqisign/sqisign.jl")

# Fp2 and values in Fp2
function make_field_curve_torsions()
    _, T = polynomial_ring(GF(p), "T")
    Fp2, Fp2_i = finite_field(T^2 + 1, "i")
    
    A0 = Fp2(0)

    # constatns from precompute/toy17.sage
    P2e = Point(126351*Fp2_i + 57101, 75454*Fp2_i + 34373)
    Q2e = Point(125461*Fp2_i + 44066, 7433*Fp2_i + 62895)
    M_i_2e = BigInt[74133 82715; 109618 56939]
    M_ij_2e = BigInt[119801 127723; 3501 11271]
    M_1k_2e = BigInt[78125 75652; 94193 52948]
    xP3 = Proj1(107897*Fp2_i + 11513)
    xQ3 = Proj1(23174*Fp2_i + 11513)
    xPQ3 = Proj1(23174*Fp2_i + 119558)
    M_i_3 = [1 2; 2 2]
    M_ij_3 = [2 0; 2 1]
    M_1k_3 = [1 2; 2 0]
    xP5 = Proj1(65913*Fp2_i + 71629)
    xQ5 = Proj1(52971*Fp2_i + 103834)
    xPQ5 = Proj1(52971*Fp2_i + 27237)
    M_i_5 = [2 2; 0 3]
    M_ij_5 = [0 1; 2 0]
    M_1k_5 = [0 2; 1 1]
    xP17 = Proj1(82856*Fp2_i + 49506)
    xQ17 = Proj1(68076*Fp2_i + 66534)
    xPQ17 = Proj1(62990*Fp2_i + 68206)
    M_i_17 = [11 3; 16 6]
    M_ij_17 = [16 3; 8 1]
    M_1k_17 = [14 2; 15 4]

    a24_0 = A_to_a24(A0)
    xP2e = Proj1(P2e.X, P2e.Z)
    xQ2e = Proj1(Q2e.X, Q2e.Z)
    PQ2e = add(P2e, -Q2e, Proj1(A0))
    xPQ2e = Proj1(PQ2e.X, PQ2e.Z)
    xP2e_short = xDBLe(xP2e, a24_0, ExponentForIsogeny)
    xQ2e_short = xDBLe(xQ2e, a24_0, ExponentForIsogeny)
    xPQ2e_short = xDBLe(xPQ2e, a24_0, ExponentForIsogeny)
    wp_P2e_Q2e = Weil_pairing_2power(A0, P2e, Q2e, ExponentFull)

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

    return Fp2, Fp2_i, CurveData(A0, A0d, A0dd, a24_0, jInvariant_A(A0), P2e, Q2e, xP2e, xQ2e, xPQ2e, xP2e_short, xQ2e_short, xPQ2e_short, wp_P2e_Q2e, DegreesOddTorsionBases, ExponentsOddTorsionBases, OddTorsionBases, Matrices_2e, Matrices_odd, isomorphism_to_A0)
end