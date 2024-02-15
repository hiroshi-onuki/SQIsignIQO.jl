using Nemo
import KaniSQIsign: QOrderElem, norm, LeftIdeal, FullRepresentInteger

p = BigInt(2)^247*79 - 1

N = 103
e = Int(ceil(log(2, p/N))) + 5
a, found = FullRepresentInteger(N*BigInt(2)^e, p)
if found
    println(a)
    println(div(norm(a), BigInt(2)^e))
end

I = LeftIdeal(a, N)
println(norm(I))
