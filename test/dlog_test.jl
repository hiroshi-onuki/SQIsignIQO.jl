using Nemo
import KaniSQIsign: fq_dlog_power_of_2


p = BigInt(2^8 * 3 * 5 * 7 - 1)
R, T = polynomial_ring(GF(p), "T")
Fp2, i = finite_field(T^2 + 1, "i")

x = rand(Fp2)
base = rand(Fp2)
e = 8
n = fq_dlog_power_of_2(x, base, e)
@assert x == base^n
