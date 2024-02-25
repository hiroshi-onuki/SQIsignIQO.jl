using Nemo
import KaniSQIsign

e1 = KaniSQIsign.Level1.ExponentForIsogeny
e2 = KaniSQIsign.Level1.ExponentForTorsion
f = KaniSQIsign.Level1.Cofactor

function check_short_element(param::Module)
    N = 101 * f * 3 * 5
    e = Int(ceil(log(2, p))) + e1

    a = param.QOrderElem(1)
    while param.norm(a) % (BigInt(2)^e1) != 0
        a, found = param.FullRepresentInteger(N*BigInt(2)^e)
        if found
            while param.gcd(a) != 1
                a = div(a, param.gcd(a))
            end
        end
    end

    I = param.LeftIdeal(a, BigInt(2)^e1 * 3 * 5 * f)
    @assert norm(I) == BigInt(2)^e1 * 3 * 5 * f

    x, a, b, found = param.two_e_good_element(I, e2)
    if found
        @assert param.norm(x) % norm(I) == 0
        @assert a^2 + b^2 == BigInt(2)^e2 - div(param.norm(x), param.norm(I))
    else
        println("Not found")
    end
end

for _ in 1:1
    check_short_element(KaniSQIsign.Level1)
end
println("Short ideal tests passed!")
