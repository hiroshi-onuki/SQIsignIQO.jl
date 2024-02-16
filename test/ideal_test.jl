using Nemo
import KaniSQIsign
import KaniSQIsign.Level1: QOrderElem, norm, LeftIdeal, FullRepresentInteger, two_e_good_element

e1 = KaniSQIsign.Level1.ExponentForIsogeny
e2 = KaniSQIsign.Level1.ExponentForTorsion
f = KaniSQIsign.Level1.Cofactor

N = 101 * f
e = Int(ceil(log(2, p))) + 10
a, found = FullRepresentInteger(N*BigInt(2)^e)
if found
    println(a)
    println(div(norm(a), BigInt(2)^e))
end
while gcd(a) != 1
    global a = div(a, gcd(a))
end

I = LeftIdeal(a, BigInt(2)^e1 * f)
println(factor(ZZ(norm(I))))
@assert norm(I) == BigInt(2)^e1 * f

x, a, b, found = two_e_good_element(I, e2)
println(x)
println(factor(ZZ(a^2 + b^2)))
@assert norm(x) % norm(I) == 0
@assert a^2 + b^2 == BigInt(2)^e2 - div(norm(x), norm(I))
println(log(2, div(norm(x), norm(I))))
