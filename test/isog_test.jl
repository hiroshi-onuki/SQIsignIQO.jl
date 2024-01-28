using Nemo
import KaniSQIsign: Weil_pairing_2power, random_point, random_point_order_2power,
    Proj1, Point, odd_isogeny, is_infinity, ladder, Montgomery_coeff

function basis_2power_torsion(A::T, e::Integer) where T <: RingElem
    p = characteristic(parent(A))
    @assert (p + 1) % ZZ(2)^e == 0
    xP = random_point_order_2power(A, p + 1, e)
    P = Point(A, xP)
    while true
        xQ = random_point_order_2power(A, p + 1, e)
        Q = Point(A, xQ)
        weil = Weil_pairing_2power(A, P, Q, e)
        if weil^(ZZ(2)^(e-1)) != 1
            @assert weil^(ZZ(2)^e) == 1
            return P, Q
        end
    end
end

p = ZZ(6143)
f = 3
e = 11
@assert p + 1 == ZZ(2)^e * f
R, T = polynomial_ring(GF(p), "T")
Fp2, i = FiniteField(T^2 + 1, "i")

A1 = Fp2(0)
a24 = Proj1(A1 + 2, Fp2(4))
K = random_point(A1)
K = ladder(div(p+1, 3), K, a24)
while is_infinity(K)
    global K = random_point(A1)
    global K = ladder(div(p+1, 3), K, a24)
end
P1, Q1 = basis_2power_torsion(A1, 2)
xP1 = Proj1(P1.X, P1.Z)
xQ1 = Proj1(Q1.X, Q1.Z)
a24_2, (xP2, xQ2) = odd_isogeny(a24, K, 3, [xP1, xQ1])
A2 = Montgomery_coeff(a24_2)
P2 = Point(A2, xP2)
Q2 = Point(A2, xQ2)
@assert Weil_pairing_2power(A1, P1, Q1, 2) * Weil_pairing_2power(A2, P2, Q2, 2) == 1 || Weil_pairing_2power(A1, P1, Q1, 2) == Weil_pairing_2power(A2, P2, Q2, 2)
