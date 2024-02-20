using Nemo
import KaniSQIsign
import KaniSQIsign.Level1: QOrderElem, norm, LeftIdeal, FullRepresentInteger, two_e_good_element

e1 = KaniSQIsign.Level1.ExponentForIsogeny
e2 = KaniSQIsign.Level1.ExponentForTorsion
f = KaniSQIsign.Level1.Cofactor

function check_short_element()
    N = 101 * f * 3 * 5
    e = Int(ceil(log(2, p))) + e1

    a = QOrderElem(1)
    while norm(a) % (BigInt(2)^e1) != 0
        a, found = FullRepresentInteger(N*BigInt(2)^e)
        if found
            while gcd(a) != 1
                a = div(a, gcd(a))
            end
        end
    end

    I = LeftIdeal(a, BigInt(2)^e1 * 3 * 5 * f)
    @assert norm(I) == BigInt(2)^e1 * 3 * 5 * f

    x, a, b, found = two_e_good_element(I, e2)
    if found
        @assert norm(x) % norm(I) == 0
        @assert a^2 + b^2 == BigInt(2)^e2 - div(norm(x), norm(I))
    else
        println("Not found")
    end
end

for _ in 1:100
    check_short_element()
end
println("All tests passed!")
