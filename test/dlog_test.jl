using Nemo
using KaniSQIsign
import KaniSQIsign.Level1: ec_bi_dlog_E0

p = KaniSQIsign.Level1.p
global_data = KaniSQIsign.Level1.make_precomputed_values()
E0 = global_data.E0

A = E0.A0
P0, Q0 = E0.P2e, E0.Q2e
e = KaniSQIsign.Level1.ExponentFull

P, Q = basis_2power_torsion(A, e)
PQ = add(P, Q, Proj1(A))
xP = Proj1(P.X, P.Z)
xQ = Proj1(Q.X, Q.Z)
xPQ = Proj1(PQ.X, PQ.Z)

n1, n2, n3, n4 = ec_bi_dlog_E0(xP, xQ, xPQ, E0)
n1, n2, n3, n4 = @time ec_bi_dlog_E0(xP, xQ, xPQ, E0)
n1P0 = mult(n1, P0, Proj1(A))
n2Q0 = mult(n2, Q0, Proj1(A))
n3P0 = mult(n3, P0, Proj1(A))
n4Q0 = mult(n4, Q0, Proj1(A))
@test P == add(n1P0, n2Q0, Proj1(A)) || P == add(-n1P0, -n2Q0, Proj1(A))
@test Q == add(n3P0, n4Q0, Proj1(A)) || Q == add(-n3P0, -n4Q0, Proj1(A))

E0d = global_data.orders_data[1]
A = E0d.A
e = KaniSQIsign.Level1.ExponentFull
xP, xQ, xPQ = E0d.xP2e, E0d.xQ2e, E0d.xPQ2e
P = Point(A, xP)
Q = Point(A, xQ)
tp_table = make_pairing_table(A, P, e)

t1 = Tate_pairing_P0(Q, tp_table, KaniSQIsign.Level1.Cofactor)
Q = mult(5, Q, Proj1(A))
t2 = Tate_pairing_P0(Q, tp_table, KaniSQIsign.Level1.Cofactor)
@assert t1^5 == t2

