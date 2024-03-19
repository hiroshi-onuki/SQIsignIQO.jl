using Nemo
import KaniSQIsign: Proj1, random_point, ladder, x_add_sub, torsion_basis, xDBLe, is_infinity,
    two_e_iso, compute_strategy

p = BigInt(21503)
R, T = polynomial_ring(GF(p), "T")
Fp2, i = finite_field(T^2 + 1, "i")

A = Fp2(6)
a24 = Proj1(A + 2, Fp2(4))
P = random_point(A)

Q1 = ladder(3, P, a24)
Q2 = ladder(5, P, a24)

R = x_add_sub(Q1, Q2, a24)
@assert R == ladder(8, P, a24) || R == ladder(2, P, a24)

P, Q = torsion_basis(a24, 10)
@assert is_infinity(xDBLe(P, a24, 10))
@assert is_infinity(xDBLe(Q, a24, 10))
@assert !is_infinity(xDBLe(P, a24, 9))
@assert !is_infinity(xDBLe(Q, a24, 9))
@assert xDBLe(P, a24, 9) != xDBLe(Q, a24, 9)

strategy = compute_strategy(4, 1, 1)
ker, _, _ = torsion_basis(a24, 10)
a24d, im = two_e_iso(a24, ker, 10, [Q1])
a24dd, imd = two_e_iso(a24, ker, 10, [Q1], strategy)
@assert a24d == a24dd
@assert im[1] == imd[1]
