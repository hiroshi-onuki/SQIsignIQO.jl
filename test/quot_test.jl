using Nemo
import KaniSQIsign: QOrderElem, norm, involution

p = 21503

x = QOrderElem(1, 2, 3, 4, p)
y = QOrderElem(1, 3, 12, 1, p)
@assert x*involution(x) == QOrderElem(norm(x), 0, 0, 0, p)
@assert norm(x) * norm(y) == norm(x*y)