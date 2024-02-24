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

out_file = open("level1torsion.txt", "w")

# 2^e-torsion in E(Fp2)
P, Q = tp.basis(E0, Fp2, False, 2, e)
while not (2^(e-2)*P).xy()[0]^2 == 1:
    P, Q = tp.basis(E0, Fp2, False, 2, e)
Ms = end.action_matrices([P, Q], 2^e, zeta4, Fp4)
Px, Py = [tp.Fp2ToFp2d(v, zeta4, Fp2_i) for v in P.xy()]
Qx, Qy = [tp.Fp2ToFp2d(v, zeta4, Fp2_i) for v in Q.xy()]
out_file.write("P2e = Point(%s, %s)\n" % (Px, Py))
out_file.write("Q2e = Point(%s, %s)\n" % (Qx, Qy))
out_file.write("M_i_2e = %s\n" % matrix_to_str(Ms[0]))
out_file.write("M_ij_2e = %s\n" % matrix_to_str(Ms[1]))
out_file.write("M_1k_2e = %s\n" % matrix_to_str(Ms[2]))

# odd torsion in E(Fp2)
P79, Q79 = tp.basis(E0, Fp2, False, 79, 1)
Ms = end.action_matrices([P79, Q79], 79, zeta4, Fp4)
xP79 = tp.Fp2ToFp2d(P79.xy()[0], zeta4, Fp2_i)
xQ79 = tp.Fp2ToFp2d(Q79.xy()[0], zeta4, Fp2_i)
xPQ79 = tp.Fp2ToFp2d((P79 - Q79).xy()[0], zeta4, Fp2_i)
out_file.write("xP79 = Proj1(%s)\n" % xP79)
out_file.write("xQ79 = Proj1(%s)\n" % xQ79)
out_file.write("xPQ79 = Proj1(%s)\n" % xPQ79)
out_file.write("M_i_79 = %s\n" % matrix_to_str(Ms[0]))
out_file.write("M_ij_79 = %s\n" % matrix_to_str(Ms[1]))
out_file.write("M_1k_79 = %s\n" % matrix_to_str(Ms[2]))

# odd torsion in E^t(Fp2)
P3, Q3 = tp.basis(E0, Fp2, True, 3, 1)
Ms = end.action_matrices([P3, Q3], 3, zeta4, Fp4)
xP3 = tp.Fp2ToFp2d(P3.xy()[0], zeta4, Fp2_i)
xQ3 = tp.Fp2ToFp2d(Q3.xy()[0], zeta4, Fp2_i)
xPQ3 = tp.Fp2ToFp2d((P3 - Q3).xy()[0], zeta4, Fp2_i)
out_file.write("xP3 = Proj1(%s)\n" % xP3)
out_file.write("xQ3 = Proj1(%s)\n" % xQ3)
out_file.write("xPQ3 = Proj1(%s)\n" % xPQ3)
out_file.write("M_i_3 = %s\n" % matrix_to_str(Ms[0]))
out_file.write("M_ij_3 = %s\n" % matrix_to_str(Ms[1]))
out_file.write("M_1k_3 = %s\n" % matrix_to_str(Ms[2]))

P5, Q5 = tp.basis(E0, Fp2, True, 5, 1)
Ms = end.action_matrices([P5, Q5], 5, zeta4, Fp4)
xP5 = tp.Fp2ToFp2d(P5.xy()[0], zeta4, Fp2_i)
xQ5 = tp.Fp2ToFp2d(Q5.xy()[0], zeta4, Fp2_i)
xPQ5 = tp.Fp2ToFp2d((P5 - Q5).xy()[0], zeta4, Fp2_i)
out_file.write("xP5 = Proj1(%s)\n" % xP5)
out_file.write("xQ5 = Proj1(%s)\n" % xQ5)
out_file.write("xPQ5 = Proj1(%s)\n" % xPQ5)
out_file.write("M_i_5 = %s\n" % matrix_to_str(Ms[0]))
out_file.write("M_ij_5 = %s\n" % matrix_to_str(Ms[1]))
out_file.write("M_1k_5 = %s\n" % matrix_to_str(Ms[2]))

out_file.close()