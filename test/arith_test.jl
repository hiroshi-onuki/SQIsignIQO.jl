using Nemo
import KaniSQIsign: Proj1, random_point, xDBL,
    Montgomery_to_theta, Montgomery_point_to_theta,
    product_theta_null, product_theta_pt, ladder, dbl

p = ZZ(103)
R, T = polynomial_ring(GF(p), "T")
Fp2, i = FiniteField(T^2 + 1, "i")
A1 = Fp2(0)
A2 = Fp2(6)
t1 = Montgomery_to_theta(A1)
t2 = Montgomery_to_theta(A2)

P1 = random_point(A1)
Q1 = xDBL(P1, Proj1(A1 + 2, Fp2(4)))
P2 = random_point(A2)
Q2 = xDBL(P2, Proj1(A2 + 2, Fp2(4)))
tP1 = Montgomery_point_to_theta(t1, P1)
tQ1 = Montgomery_point_to_theta(t1, Q1)
tP2 = Montgomery_point_to_theta(t2, P2)
tQ2 = Montgomery_point_to_theta(t2, Q2)

tnull = product_theta_null(t1, t2)
P = product_theta_pt(tP1, tP2)
Q = product_theta_pt(tQ1, tQ2)

@test ladder(tnull, p + 1, P) == tnull
@test Q == dbl(tnull, P)