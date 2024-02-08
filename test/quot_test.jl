using Nemo
import KaniSQIsign: QOrderElem, norm

p = 21503

x = QOrderElem(ZZ(1), ZZ(2), ZZ(3), ZZ(4), ZZ(p))
y = QOrderElem(ZZ(5), ZZ(6), ZZ(7), ZZ(8), ZZ(p))
@assert norm(x) * norm(y) == norm(x*y)