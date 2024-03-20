using Nemo
import Pkg
Pkg.activate(@__DIR__)
using KaniSQIsign
using KaniSQIsign.Level1

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

Fp2, Fp2_i, E0 = Level1.make_field_curve_torsions()
order2e = BigInt(2)^Level1.ExponentFull
d = 2
a, b, c, N = sqrt_in_quaternion(d, Level1)
alpha = Level1.order_elem_from_standard_basis(0, a, b, c)
Ma = alpha[1] * [1 0; 0 1] + alpha[2] * E0.Matrices_2e[1] + alpha[3] * E0.Matrices_2e[2] + alpha[4] * E0.Matrices_2e[3]
println(alpha * alpha)
@assert alpha * alpha == Level1.QOrderElem(-d * N^2)
I = Level1.LeftIdeal(alpha, N)
alpha, found = @time Level1.KLPT(I, N)
alpha = div(alpha, Level1.gcd(alpha))
@assert found
J = Level1.ideal_transform(I, alpha, N)
println(factor(ZZ(Level1.norm(J))))
e = Int(log(2, div(Level1.norm(J), Level1.ExtraDegree)))
a24, xP, xQ, xPQ, M = Level1.ideal_to_isogeny_from_O0(J, e, E0)
println(jInvariant_a24(a24))
Minv = [M[2, 2] -M[1, 2]; -M[2, 1] M[1, 1]] * invmod(M[1, 1]*M[2, 2] - M[1, 2]*M[2, 1], order2e)
Msqrt2 = (M * Ma * Minv * invmod(N, order2e)) .% order2e
println((Msqrt2 * Msqrt2 - [-2 0; 0 -2]) .% order2e)

for _ in 1:100
    n = rand(BigInt(2)^128:BigInt(2)^130)
    x, y, found_xy = Level1.sum_of_two_squares_2(n)
    if found_xy
        @assert x^2 + 2*y^2 == n
        println("n = $n, x = $x, y = $y")
    end
end