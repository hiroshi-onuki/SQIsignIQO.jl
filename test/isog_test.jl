using Nemo
import KaniSQIsign: Proj1, random_point, xDBL, is_infinity,
    Montgomery_to_theta, Montgomery_point_to_theta,
    product_theta_null, product_theta_pt, ladder, dbl,
    two_two_isogeny_8torsion, ThetaLv1, theta_to_Montgomery, jInvariant, Hadamard,
    level_22_constants_sqr, CouplePoint, InfPoint, gluing_isogeny

function order8point(A::T, p::ZZRingElem) where T <: RingElem
    @assert p % 8 == 7
    F = parent(A)
    a24 = Proj1(A + 2, F(4))
    while true
        P = random_point(A)
        n = div(p + 1, 8)
        P = ladder(n, P, a24)
        if !is_infinity(ladder(ZZ(4), P, a24))
            return P
        end
    end
end

p = ZZ(5119)
R, T = polynomial_ring(GF(p), "T")
Fp2, i = FiniteField(T^2 + 1, "i")

A1 = Fp2(0)
A2 = Fp2(6)
a24_1 = Proj1(A1 + 2, Fp2(4))
a24_2 = Proj1(A2 + 2, Fp2(4))
T1 = order8point(A1, p)
T2 = order8point(A2, p)

P1P2 = CouplePoint(T1, T2)
Q1Q2 = CouplePoint(T1, T2)

codomain, _ = gluing_isogeny(a24_1, a24_2, P1P2, Q1Q2, [P1P2])

even_indices = [
            [0, 0],
            [0, 1],
            [0, 2],
            [0, 3],
            [1, 0],
            [1, 2],
            [2, 0],
            [2, 1],
            [3, 0],
            [3, 3],
        ]
for (chi, i) in even_indices
    println(level_22_constants_sqr(codomain, chi, i))
end

A3 = theta_to_Montgomery(ThetaLv1(a, b))
A4 = theta_to_Montgomery(ThetaLv1(a, c))
println(jInvariant(Proj1((A3 + 2)/4)))
println(256*(A3^2 - 3)^3/(A3^2 - 4))
println(jInvariant(Proj1((A4 + 2)/4)))
