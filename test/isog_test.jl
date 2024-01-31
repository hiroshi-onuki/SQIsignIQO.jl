using Nemo
import KaniSQIsign: Weil_pairing_2power, random_point, random_point_order_2power,
    Proj1, Point, odd_isogeny, is_infinity, ladder, Montgomery_coeff,
    xDBLe, CouplePoint, gluing_isogeny, product_isogeny_no_strategy

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

p = ZZ(21503)
e = 3
N = 7
M = invmod(-N, 2^(e+2))
R, T = polynomial_ring(GF(p), "T")
Fp2, i = FiniteField(T^2 + 1, "i")

A1 = Fp2(0)
a24_1 = Proj1(A1 + 2, Fp2(4))
K = random_point(A1)
K = ladder(div(p+1, N), K, a24_1)
while is_infinity(K)
    global K = random_point(A1)
    global K = ladder(div(p+1, N), K, a24_1)
end
P1, Q1 = basis_2power_torsion(A1, e + 2)
xP1 = Proj1(P1.X, P1.Z)
xQ1 = Proj1(Q1.X, Q1.Z)
a24_2, (xP2, xQ2) = odd_isogeny(a24_1, K, N, [xP1, xQ1])
A2 = Montgomery_coeff(a24_2)
P2 = Point(A2, xP2)
xQ2 = ladder(ZZ(M), xQ2, a24_2)
Q2 = Point(A2, xQ2)
@assert Weil_pairing_2power(A1, P1, Q1, e + 2) == Weil_pairing_2power(A2, P2, Q2, e + 2) || Weil_pairing_2power(A1, P1, Q1, e + 2)*Weil_pairing_2power(A2, P2, Q2, e + 2) == 1

xP1 = Proj1(9869*i)
xP2 = Proj1(13690*i + 4536)
xQ1 = Proj1(15061*i + 18352)
xQ2 = Proj1(3576*i + 5041)
a24_2 = Proj1(17409*i+4362 + 2, Fp2(4))
P1P2 = CouplePoint(xP1, xP2)
Q1Q2 = CouplePoint(xQ1, xQ2)
codomain, images = product_isogeny_no_strategy(a24_1, a24_2, P1P2, Q1Q2, CouplePoint{FqFieldElem}[], e)

