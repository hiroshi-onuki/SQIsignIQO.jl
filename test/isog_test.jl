using Nemo
import KaniSQIsign: Weil_pairing_2power, random_point_order_2power, Point

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
A2 = Fp2(6)

P1, Q1 = basis_2power_torsion(A1, e)
