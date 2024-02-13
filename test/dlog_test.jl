using Nemo
import KaniSQIsign: fq_dlog_power_of_2


p = BigInt(2^8 * 3 * 5 * 7 - 1)
R, T = polynomial_ring(GF(p), "T")
Fp2, i = finite_field(T^2 + 1, "i")
e = 9

base = rand(Fp2)^div(p^2 - 1, 2^e)
while base^(2^(e-1)) == 1
    global base = rand(Fp2)^div(p^2 - 1, 2^e)
end

x = rand(Fp2)^div(p^2 - 1, 2^e)
@assert base^(2^e) == 1
@assert x^(2^e) == 1

n = fq_dlog_power_of_2(x, base, e)
@assert x == base^n
