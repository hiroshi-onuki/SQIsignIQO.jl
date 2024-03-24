using KaniSQIsign

function signing_test(param::Module, num::Int)
    _, _, global_data = param.make_precomputed_values()
    println("Signing test for $(param)")
    for _ in 1:num
        println("Generate keys")
        pk, sk, found = @time param.key_gen(global_data)

        println("Sign message")
        m = "message to sign"
        sign = @time param.signing(pk, sk, m, global_data)
        println("Signature: $(length(sign)) bytes")

        println("Verify signature")
        verif = @time param.verify(pk, m, sign)
        @assert verif
        println("Verified")
    end
end

signing_test(KaniSQIsign.Level1, 10)
signing_test(KaniSQIsign.Level3, 10)
signing_test(KaniSQIsign.Level5, 10)