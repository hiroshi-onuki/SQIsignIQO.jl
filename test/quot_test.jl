using Nemo
import KaniSQIsign.Level1: QOrderElem, norm, involution

x = QOrderElem(1, 2, 3, 4)
y = QOrderElem(1, 3, 12, 1)
@assert x*involution(x) == QOrderElem(norm(x), 0, 0, 0)
@assert norm(x) * norm(y) == norm(x*y)