using Nemo
import KaniSQIsign: Proj1, Point, add, double, mult, infinity_full_point, random_point, affine


p = BigInt(21503)
R, T = polynomial_ring(GF(p), "T")
Fp2, i = finite_field(T^2 + 1, "i")

A = Proj1(Fp2(12*i), Fp2(2*i))
xP = random_point(affine(A))
P = Point(affine(A), xP)
@test mult(p + 1, P, A) == infinity_full_point(Fp2)