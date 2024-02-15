using Nemo
import KaniSQIsign: QOrderElem, norm, LeftIdeal, FullRepresentInteger, small_element

p = BigInt(2)^247*79 - 1
p = 103

N = 13
e = Int(ceil(log(2, p/N))) + 5
a, found = FullRepresentInteger(N*BigInt(2)^e, p)
if found
    println(a)
    println(div(norm(a), BigInt(2)^e))
end

I = LeftIdeal(a, N)
println(norm(I))
x = small_element(I)
println(x)
@assert norm(x) % norm(I) == 0
println(log(2, div(norm(x), norm(I))))
