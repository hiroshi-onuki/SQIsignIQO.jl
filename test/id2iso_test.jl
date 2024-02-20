using Nemo
using KaniSQIsign
using KaniSQIsign.Level1
import KaniSQIsign.Level1: short_ideal_to_isogeny

p = KaniSQIsign.Level1.p
e1 = KaniSQIsign.Level1.ExponentForIsogeny
e2 = KaniSQIsign.Level1.ExponentForTorsion
Fp2, i, tdata = KaniSQIsign.Level1.make_field_curve_torsions()
a24_0 = tdata.a24_0
xP0 = tdata.xPe2
xQ0 = tdata.xQe2
xPQ0 = tdata.xPQe2

N = 3 * 5 * 79
e = Int(ceil(log(2, p))) + e1

a = QOrderElem(1)
a, found = FullRepresentInteger(N*BigInt(2)^e)
a = div(a, gcd(a))

I = LeftIdeal(a, BigInt(2)^e1 * 3 * 5 * 79)
while gcd(I) != 1
    global I = div(I, gcd(I))
end
println(factor(ZZ(norm(I))))

M0 = BigInt[1 0; 0 1]
a24 = short_ideal_to_isogeny(I, a24_0, xP0, xQ0, xPQ0, M0, tdata)