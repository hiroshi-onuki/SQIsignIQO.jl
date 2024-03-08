using Nemo
import KaniSQIsign: Proj1, ec_dlog_power_of_2, random_point_order_2power, add, mult, make_dlog_table

p = BigInt(2)^247*79 - 1
R, T = polynomial_ring(GF(p), "T")
Fp2, i = finite_field(T^2 + 1, "i")
e = 247

A = Fp2(0)
P, Q = basis_2power_torsion(A, e)
R, S = basis_2power_torsion(A, e)
PQ = add(P, Q, Proj1(A))

xP = Proj1(P.X, P.Z)
xQ = Proj1(Q.X, Q.Z)
xPQ = Proj1(PQ.X, PQ.Z)

n1, n2, n3, n4 = @time ec_dlog_power_of_2(xP, xQ, xPQ, R, S, A, e)
n1R = mult(n1, R, Proj1(A))
n2S = mult(n2, S, Proj1(A))
n3R = mult(n3, R, Proj1(A))
n4S = mult(n4, S, Proj1(A))
@test P == add(n1R, n2S, Proj1(A)) || P == add(-n1R, -n2S, Proj1(A))
@test Q == add(n3R, n4S, Proj1(A)) || Q == add(-n3R, -n4S, Proj1(A))

T1, T2 = make_dlog_table(i, e, 3)
