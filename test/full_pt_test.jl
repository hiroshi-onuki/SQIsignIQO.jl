using Nemo
import KaniSQIsign: Proj1, Point, add, double


p = ZZ(21503)
R, T = polynomial_ring(GF(p), "T")
Fp2, i = FiniteField(T^2 + 1, "i")

A = Proj1(Fp2(12*i), Fp2(2*i))

P = Point(Fp2(4), 9287*i)
Q = Point(Fp2(3), Fp2(1344))
PQ = add(P, Q, A)
println(PQ.X/PQ.Z, " ", PQ.Y/PQ.Z)
P2 = double(P, A)
println(P2.X/P2.Z, " ", P2.Y/P2.Z)
Q2 = double(Q, A)
println(Q2.X/Q2.Z, " ", Q2.Y/Q2.Z)
P4 = double(P2, A)
println(P4.X/P4.Z, " ", P4.Y/P4.Z)
P2Q2 = add(P2, Q2, A)
println(P2Q2.X/P2Q2.Z, " ", P2Q2.Y/P2Q2.Z)