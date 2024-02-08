using Nemo
import KaniSQIsign: Proj1, random_point, ladder, x_add_sub

p = ZZ(21503)
R, T = polynomial_ring(GF(p), "T")
Fp2, i = finite_field(T^2 + 1, "i")

A = Fp2(6)
a24 = Proj1(A + 2, Fp2(4))
P = random_point(A)

Q1 = ladder(ZZ(3), P, a24)
Q2 = ladder(ZZ(5), P, a24)

R = x_add_sub(Q1, Q2, a24)
@assert R == ladder(ZZ(8), P, a24) || R == ladder(ZZ(2), P, a24)