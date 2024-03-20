using Nemo
import KaniSQIsign

function check_short_element(param::Module)

    e1 = param.ExponentForIsogeny
    e2 = param.ExponentForTorsion
    ext_factor = param.ExtraDegree

    N = 101 * ext_factor
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

    I = param.LeftIdeal(a, BigInt(2)^e1 * ext_factor)
    @assert param.norm(I) == BigInt(2)^e1 * ext_factor

    cor_func(argN) = param.sum_of_two_squares(BigInt(2)^e2 - argN)
    x, a, b, found = param.two_e_good_element(I, norm(I), cor_func, norm(I) << e2)
    if found
        @assert param.norm(x) % param.norm(I) == 0
        @assert a^2 + b^2 == BigInt(2)^e2 - div(param.norm(x), param.norm(I))
    else
        println("Not found")
    end
end

for _ in 1:10
    check_short_element(KaniSQIsign.Level1)
end
println("Short ideal tests passed!")
