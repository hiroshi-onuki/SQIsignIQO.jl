using Nemo
import KaniSQIsign: QOrderElem, norm, involution

p = ZZ(21503)

x = QOrderElem(ZZ(1), ZZ(2), ZZ(3), ZZ(4), p)
y = QOrderElem(ZZ(1), ZZ(3), ZZ(12), ZZ(1), p)
println(x*involution(x))
@assert x*involution(x) == QOrderElem(norm(x), ZZ(0), ZZ(0), ZZ(0), p)
@assert norm(x) * norm(y) == norm(x*y)