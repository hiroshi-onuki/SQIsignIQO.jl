import endomorphism as end
import torsion_points as tp

def matrix_to_str(M):
    return "[%s %s; %s %s]" % (M[0,0], M[0,1], M[1,0], M[1,1])

set_random_seed(0)

p = 2^247 * 79 - 1
e = 247
Fp4, Fp2, zeta4 = tp.calcFields(p)
E0 = EllipticCurve(Fp4, [1, 0])
Fp2d.<Fp2_i> = GF(p^2, modulus=x^2+1)

P, Q = tp.basis(E0, Fp2, False, 2, e)
Ms = end.action_matrices([P, Q], 2^e, zeta4, Fp4)
Px, Py = [tp.Fp2ToFp2d(v, zeta4, Fp2_i) for v in P.xy()]
Qx, Qy = [tp.Fp2ToFp2d(v, zeta4, Fp2_i) for v in Q.xy()]

with open("level1torsion.txt", "w") as f:
    f.write("const P2e = Point(%s, %s)\n" % (Px, Py))
    f.write("const Q2e = Point(%s, %s)\n" % (Qx, Qy))
    f.write("const M_i_2e = %s\n" % matrix_to_str(Ms[0]))
    f.write("const M_ij_2e = %s\n" % matrix_to_str(Ms[1]))
    f.write("const M_1k_2e = %s\n" % matrix_to_str(Ms[2]))

# odd torsion in E(Fp2)
P79, Q79 = tp.basis(E0, Fp2, False, 79, 1)
Ms = end.action_matrices([P79, Q79], 79, zeta4, Fp4)


P3, Q3 = tp.basis(E0, Fp2, True, 3, 1)
Ms = end.action_matrices([P3, Q3], 3, zeta4, Fp4)
print([det(M) % 3 for M in Ms])

P5, Q5 = tp.basis(E0, Fp2, True, 5, 1)
Ms = end.action_matrices([P5, Q5], 5, zeta4, Fp4)
print([det(M) % 5 for M in Ms])

