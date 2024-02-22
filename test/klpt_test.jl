import KaniSQIsign
import KaniSQIsign.Level1: LeftIdeal, FullRepresentInteger, RandomEquivalentPrimeIdeal

p = KaniSQIsign.Level1.p

N = BigInt(2)^1000
a, found = FullRepresentInteger(N)
a = div(a, gcd(a))
@assert norm(a) % BigInt(2)^256 == 0
I = LeftIdeal(a, BigInt(2)^256)
println(factor(ZZ(norm(I))))
J, found = RandomEquivalentPrimeIdeal(I)
println(factor(ZZ(norm(J))))
