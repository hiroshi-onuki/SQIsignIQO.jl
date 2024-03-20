using Nemo
import Pkg
Pkg.activate(@__DIR__)
using KaniSQIsign
using KaniSQIsign.Level1

struct OtherOrderData
    A::FqFieldElem
    a24_0::Proj1{FqFieldElem}
    j0::FqFieldElem
    P2e::Point{FqFieldElem}
    Q2e::Point{FqFieldElem}
    xP2e::Proj1{FqFieldElem}
    xQ2e::Proj1{FqFieldElem}
    xPQ2e::Proj1{FqFieldElem}
    xP2e_short::Proj1{FqFieldElem}
    xQ2e_short::Proj1{FqFieldElem}
    xPQ2e_short::Proj1{FqFieldElem}
    M::Matrix{BigInt}
    connecting_deg::BigInt
end

# return alpha = (ai + bj + cij)/N s.t. alpha^2 = -d
function sqrt_in_quaternion(d::Integer, param::Module)
    p = param.p
    while true
        N = rand(p+1:10*p)
        !is_prime(N) && continue
        a = (BigInt(sqrtmod(ZZ(d), ZZ(p))) * N) % p
        M = div(d * N^2 - a^2, p)
        b, c, found = param.sum_of_two_squares(M)
        if found
            return a, b, c, N
        end
    end
end

Fp2, Fp2_i, cdata = Level1.make_field_curve_torsions()
d = 2
a, b, c, N = sqrt_in_quaternion(d, Level1)
alpha = Level1.order_elem_from_standard_basis(0, a, b, c)
println(alpha * alpha)
@assert alpha * alpha == Level1.QOrderElem(-d * N^2)
I = Level1.LeftIdeal(alpha, N)
alpha, found = @time Level1.KLPT(I, N)
alpha = div(alpha, Level1.gcd(alpha))
@assert found
J = Level1.ideal_transform(I, alpha, N)
println(factor(ZZ(Level1.norm(J))))
e = Int(log(2, div(Level1.norm(J), Level1.ExtraDegree)))
a24, xP, xQ, xPQ, M = Level1.ideal_to_isogeny_from_O0(J, e, cdata)
println(jInvariant_a24(a24))

for _ in 1:100
    n = rand(BigInt(2)^128:BigInt(2)^130)
    x, y, found_xy = Level1.sum_of_two_squares_2(n)
    if found_xy
        @assert x^2 + 2*y^2 == n
        println("n = $n, x = $x, y = $y")
    end
end