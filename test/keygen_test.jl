using KaniSQIsign

function keygen_test(param::Module)
    Fp2, Fp2_i, cdata = param.make_field_curve_torsions()
    param.key_gen(cdata)
end

function gen_ideals(param::Module, e::Int)
    S = param.LeftIdeal[]
    for _ in 1:(1000 * 2^e)
        I = param.sample_random_ideal_2e(e)
        @assert gcd(I) == 1
        @assert param.norm(I) == BigInt(2)^e
        if !(I in S)
            push!(S, I)
        end
        println(length(S))
    end
    println("Generated $(length(S)) ideals of norm 2^$e")
end

#keygen_test(KaniSQIsign.Level1)
gen_ideals(KaniSQIsign.Level1, 5)
