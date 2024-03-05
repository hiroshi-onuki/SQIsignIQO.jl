using Nemo
import Pkg
Pkg.activate(@__DIR__)
using KaniSQIsign
using KaniSQIsign.Level1

Fp2, Fp2_i, cdata = Level1.make_field_curve_torsions()
a = -1350115946779945062706551282873903371109
b = 14196007264819976798155349539053604713891080213218451346004006012776352581
c = 1686799866919031686392406573617475223258
d = 2700231893559890125413102565747806742218
N = 150458946602379417989372322849839053131658421176126008782090514836129534625889
alpha = Level1.QOrderElem(a, b, c, d)
println(alpha * alpha)
@assert alpha * alpha == Level1.QOrderElem(-2 * N^2)
I = Level1.LeftIdeal(alpha, N)
alpha, found = @time Level1.KLPT(I, N)
#alpha = div(alpha, Level1.gcd(alpha))
@assert found
J = Level1.ideal_transform(I, alpha, N)
println(factor(ZZ(Level1.norm(J))))
e = Int(log(2, div(Level1.norm(J), Level1.ExtraDegree)))
a24, xP, xQ, xPQ, M = Level1.ideal_to_isogeny_from_O0(J, e, cdata)
println(jInvariant_a24(a24))