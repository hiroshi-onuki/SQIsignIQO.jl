using KaniSQIsign

function test_KeyGenKLPT(param::Module)
    N = BigInt(35497637494684292293)
    e = param.Log2p

    a, found = param.FullRepresentInteger(N*BigInt(2)^e)
    I = param.LeftIdeal(a, N)
    beta, found = param.KeyGenKLPT(I, N)
    println(found)
    println(param.isin(beta, I))
    N = div(param.norm(beta), param.gcd(beta))
    println(factor(ZZ(N)))
    return N == BigInt(2)^param.KLPT_keygen_length
end

@test test_KeyGenKLPT(KaniSQIsign.Level1)
