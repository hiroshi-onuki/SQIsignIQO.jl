using Nemo
using KaniSQIsign

function test_id2iso(param::Module, n::Int)
    p = param.p
    e1 = param.ExponentForIsogenyDim1
    global_data = param.make_precomputed_values()

    N = 337416411430778090000393699717 # large prime
    e = n*Int(ceil(log(2, p)))

    a, found = param.FullRepresentInteger(N*BigInt(2)^e * param.ExtraDegree)
    if !found
        println("Not found")
        return
    end
    a = div(a, gcd(a))
    println(factor(ZZ(param.norm(a))))
    I = param.LeftIdeal(a, BigInt(2)^(n*e1) * param.ExtraDegree)
    I = div(I, param.gcd(I))
    println(factor(ZZ(param.norm(I))))
    param.ideal_to_isogeny_from_O0(I, n*e1 - 10, global_data)
end

strategy = KaniSQIsign.Level1.StrategyDim2
test_id2iso(KaniSQIsign.Level1, 4)
