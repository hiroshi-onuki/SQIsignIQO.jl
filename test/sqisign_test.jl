using KaniSQIsign

function signing_test(param::Module, num::Int)
    _, _, cdata = param.make_field_curve_torsions()
    println("Signing test for $(param)")
    for _ in 1:num
        println("Generate keys")
        pk, sk, found = @time param.key_gen(cdata)

        println("Sign message")
        m = "message to sign"
        sign, s1, s2, r = @time param.signing(pk, sk, m, cdata)

        println("Verify signature")
        verif = @time param.verify(pk, m, sign, s1, s2, r)
        println("Verified: $(verif)")
    end
end

signing_test(KaniSQIsign.Level1, 100)
