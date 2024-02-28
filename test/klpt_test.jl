import KaniSQIsign
import KaniSQIsign.Level1: LeftIdeal, FullRepresentInteger, RandomEquivalentPrimeIdeal, EichlerModConstraint, QOrderElem,
    isin, Quoternion_j, Quoternion_ij

p = KaniSQIsign.Level1.p

N = BigInt(2)^1000
a, found = FullRepresentInteger(N)
a = div(a, gcd(a))
@assert norm(a) % BigInt(2)^256 == 0
I = LeftIdeal(a, BigInt(2)^256)

J, N, found = RandomEquivalentPrimeIdeal(I)
gamma, found = FullRepresentInteger(N * BigInt(2)^256)
C, D = EichlerModConstraint(J, N, gamma, QOrderElem(1), true)
@assert isin(gamma * (C * Quoternion_j - D * Quoternion_ij), J)
