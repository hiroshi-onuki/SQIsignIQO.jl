using Nemo
import KaniSQIsign: QOrderElem, norm, LeftIdeal, FullRepresentInteger, two_e_good_element

p = BigInt(2)^247 * 79 - 1
p = BigInt(2)^377 * 3 * 7 * 11 - 1
p = BigInt(2)^496 * 5^3 - 1
e1 = 238
e2 = 258

N = 79 * 3 * 7 * 11 * 5^3
e = Int(ceil(log(2, p))) + 10
a, found = FullRepresentInteger(N*BigInt(2)^e, p)
if found
    println(a)
    println(div(norm(a), BigInt(2)^e))
end
while gcd(a) != 1
    global a = div(a, gcd(a))
end

I = LeftIdeal(a, BigInt(2)^e1 * 5^3)
println(factor(ZZ(norm(I))))
@assert norm(I) == BigInt(2)^e1 * 5^3

x, a, b, found = two_e_good_element(I, e2)
println(x)
println(factor(ZZ(a^2 + b^2)))
@assert norm(x) % norm(I) == 0
@assert a^2 + b^2 == BigInt(2)^e2 - div(norm(x), norm(I))
println(log(2, div(norm(x), norm(I))))
