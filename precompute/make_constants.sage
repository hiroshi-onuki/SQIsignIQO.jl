import endomorphism as end
import torsion_points as tp

def matrix_to_str(M):
    return "[%s %s; %s %s]" % (M[0,0], M[0,1], M[1,0], M[1,1])

def matrix44_to_str(M):
    return "[%s %s %s %s; %s %s %s %s; %s %s %s %s; %s %s %s %s]" % (
        M[0,0], M[0,1], M[0,2], M[0,3],
        M[1,0], M[1,1], M[1,2], M[1,3],
        M[2,0], M[2,1], M[2,2], M[2,3],
        M[3,0], M[3,1], M[3,2], M[3,3]
    )

def make_constants(p, e, ed, degs, degs_d, file_name):
    Fp4, Fp2, zeta4 = tp.calcFields(p)
    E0 = EllipticCurve(Fp4, [1, 0])
    Fp2d.<Fp2_i> = GF(p^2, modulus=x^2+1)

    out_file = open(file_name, "w")

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
    Msd = [[[1, 0], [0, 1]]] + Ms
    Z2e = quotient(ZZ, 2^ed)
    M44 = matrix(Z2e, [[Msd[i][j][k] for i in range(4)] for j, k in [(0, 0), (1, 0), (0, 1), (1, 1)]])
    out_file.write("M44inv = %s\n" % matrix44_to_str(M44^(-1)))

    # odd torsion in E(Fp2)
    for l, e in factor(degs):
        P, Q = tp.basis(E0, Fp2, False, l, e)
        Ms = end.action_matrices([P, Q], l, zeta4, Fp4)
        xP = tp.Fp2ToFp2d(P.xy()[0], zeta4, Fp2_i)
        xQ = tp.Fp2ToFp2d(Q.xy()[0], zeta4, Fp2_i)
        xPQ = tp.Fp2ToFp2d((P - Q).xy()[0], zeta4, Fp2_i)
        out_file.write("xP%d = Proj1(%s)\n" % (l, xP))
        out_file.write("xQ%d = Proj1(%s)\n" % (l, xQ))
        out_file.write("xPQ%d = Proj1(%s)\n" % (l, xPQ))
        out_file.write("M_i_%d = %s\n" % (l, matrix_to_str(Ms[0])))
        out_file.write("M_ij_%d = %s\n" % (l, matrix_to_str(Ms[1])))
        out_file.write("M_1k_%d = %s\n" % (l, matrix_to_str(Ms[2])))

    # odd torsion in E^t(Fp2)
    for l, e in factor(degs_d):
        P, Q = tp.basis(E0, Fp2, true, l, e)
        Ms = end.action_matrices([P, Q], l, zeta4, Fp4)
        xP = tp.Fp2ToFp2d(P.xy()[0], zeta4, Fp2_i)
        xQ = tp.Fp2ToFp2d(Q.xy()[0], zeta4, Fp2_i)
        xPQ = tp.Fp2ToFp2d((P - Q).xy()[0], zeta4, Fp2_i)
        out_file.write("xP%d = Proj1(%s)\n" % (l, xP))
        out_file.write("xQ%d = Proj1(%s)\n" % (l, xQ))
        out_file.write("xPQ%d = Proj1(%s)\n" % (l, xPQ))
        out_file.write("M_i_%d = %s\n" % (l, matrix_to_str(Ms[0])))
        out_file.write("M_ij_%d = %s\n" % (l, matrix_to_str(Ms[1])))
        out_file.write("M_1k_%d = %s\n" % (l, matrix_to_str(Ms[2])))

    out_file.close()

# toy17
set_random_seed(0)
p = 2^17 - 1
e = 17
ed = 9
degs = 1
degs_d = 3*5*17
make_constants(p, e, ed, degs, degs_d, "toy17.txt")

# level1
set_random_seed(0)
p = 2^247 * 79 - 1
e = 247
ed = 128
degs = 79
degs_d = 3*5
make_constants(p, e, ed, degs, degs_d, "level1torsion.txt")