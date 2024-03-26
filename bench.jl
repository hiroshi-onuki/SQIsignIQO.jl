using Nemo
import Pkg
Pkg.activate(@__DIR__)
using KaniSQIsign

function benchmark_test(param::Module, num::Int)
    global_data = param.make_precomputed_values()

    # for compilation
    pk, sk, found = param.key_gen(global_data)
    m = "message to sign"
    sign = param.signing(pk, sk, m, global_data)
    verif = param.verify(pk, m, sign)

    t_gen = 0
    t_sign = 0
    t_verif = 0
    println("Benchmark test for $(param) start")
    for i in 1:num
        (pk, sk, found), t, _, _, _ = @timed param.key_gen(global_data)
        t_gen += t

        m = "message to sign"
        sign, t, _, _, _ = @timed param.signing(pk, sk, m, global_data)
        t_sign += t

        verif, t, _, _, _ = @timed param.verify(pk, m, sign)
        t_verif += t

        @assert verif
        print("\r($i/$num) done.")
    end
    print("\n")
    println("Average time for key generation: ", t_gen / num)
    println("Average time for signing: ", t_sign / num)
    println("Average time for verification: ", t_verif / num)
end

benchmark_test(KaniSQIsign.Level1, 100) 
benchmark_test(KaniSQIsign.Level3, 100) 
benchmark_test(KaniSQIsign.Level5, 100) 