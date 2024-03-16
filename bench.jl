using Nemo
import Pkg
Pkg.activate(@__DIR__)
using KaniSQIsign
using KaniSQIsign.Level1

function benchmark_test(param::Module, num::Int)
    _, _, cdata = param.make_field_curve_torsions()

    # for compilation
    pk, sk, found = param.key_gen(cdata)
    m = "message to sign"
    sign, s1, s2, r = param.signing(pk, sk, m, cdata)
    verif = param.verify(pk, m, sign, s1, s2, r)

    t_gen = 0
    t_sign = 0
    t_verif = 0
    println("Benchmark test for $(param) start")
    for i in 1:num
        (pk, sk, found), t, _, _, _ = @timed param.key_gen(cdata)
        t_gen += t

        m = "message to sign"
        (sign, s1, s2, r), t, _, _, _ = @timed param.signing(pk, sk, m, cdata)
        t_sign += t

        verif, t, _, _, _ = @timed param.verify(pk, m, sign, s1, s2, r)
        t_verif += t

        @assert verif
        print("\r($i/$num) done.")
    end
    print("\n")
    println("Average time for key generation: ", t_gen / num)
    println("Average time for signing: ", t_sign / num)
    println("Average time for verification: ", t_verif / num)
end

benchmark_test(KaniSQIsign.Level1, 10)
