using Nemo
using KaniSQIsign
import KaniSQIsign.Level1: ec_bi_dlog_E0

p = KaniSQIsign.Level1.p
_, _, cdata = KaniSQIsign.Level1.make_field_curve_torsions()

A = cdata.A0
P0, Q0 = cdata.P2e, cdata.Q2e
e = KaniSQIsign.Level1.ExponentFull

P, Q = basis_2power_torsion(A, e)
PQ = add(P, Q, Proj1(A))
xP = Proj1(P.X, P.Z)
xQ = Proj1(Q.X, Q.Z)
xPQ = Proj1(PQ.X, PQ.Z)

n1, n2, n3, n4 = ec_bi_dlog_E0(xP, xQ, xPQ, cdata)
n1, n2, n3, n4 = @time ec_bi_dlog_E0(xP, xQ, xPQ, cdata)
n1P0 = mult(n1, P0, Proj1(A))
n2Q0 = mult(n2, Q0, Proj1(A))
n3P0 = mult(n3, P0, Proj1(A))
n4Q0 = mult(n4, Q0, Proj1(A))
@test P == add(n1P0, n2Q0, Proj1(A)) || P == add(-n1P0, -n2Q0, Proj1(A))
@test Q == add(n3P0, n4Q0, Proj1(A)) || Q == add(-n3P0, -n4Q0, Proj1(A))

