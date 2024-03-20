using KaniSQIsign

function signing_test(param::Module, num::Int)
    _, _, global_data = param.make_precomputed_values()
    E0 = global_data.E0
    println("Signing test for $(param)")
    for _ in 1:num
        println("Generate keys")
        pk, sk, found = @time param.key_gen(E0)

        println("Sign message")
        m = "message to sign"
        sign = @time param.signing(pk, sk, m, E0)
        println("Signature: $(length(sign)) bytes")

        println("Verify signature")
        verif = @time param.verify(pk, m, sign)
        @assert verif
        println("Verified")
    end
end

signing_test(KaniSQIsign.Level1, 10)
