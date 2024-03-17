using KaniSQIsign

function signing_test(param::Module, num::Int)
    _, _, cdata = param.make_field_curve_torsions()
    println("Signing test for $(param)")
    for _ in 1:num
        println("Generate keys")
        pk, sk, found = @time param.key_gen(cdata)

        println("Sign message")
        m = "message to sign"
        sign = @time param.signing(pk, sk, m, cdata)
        println("Signature: $(length(sign)) bytes")

        println("Verify signature")
        verif = @time param.verify(pk, m, sign)
        @assert verif
        println("Verified")
    end
end

signing_test(KaniSQIsign.Level1, 10)
