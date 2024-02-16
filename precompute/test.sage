import endomorphism as end

# return GF(p^4), GF(p^2), square root of -1 in GF(p^2)
def calcFields(p):
    assert p % 4 == 3
    R = PolynomialRing(GF(p), name="x")
    x = R.gens()[0]
    Fp2 = GF(p**2, modulus=x**2+1, name="i")
    z = Fp2.random_element()
    while z.is_square():
        z = Fp2.random_element()
    t = ZZ(z + z**p)
    n = ZZ(z**(p+1))
    R = PolynomialRing(ZZ, name="x")
    x = R.gens()[0]
    Fp4 = GF(p**4, modulus=x**4 - t*x**2 + n, name="z")
    Fp2 = Fp4.subfield(2)
    i = Fp2(-1).sqrt()
    return Fp4, Fp2, i

# return a point on E of order l^e, where l is prime.
# Assume E is a supersingular elliptic curve s.t. E(Fp^2) = (Z/(p+1)Z)^2
def point_ord(E, l, e):
    assert is_prime(l)
    p = E.base_ring().characteristic()
    assert (p + 1) % l**e == 0
    P = (p + 1)//(l**e)*E.random_point()
    while (l**(e-1)*P).is_zero():
        P = (p + 1)//(l**e)*E.random_point()
    assert (l**e*P).is_zero()
    if P[1][0] >= (p + 1)//2 or (P[1][0] == 0 and P[1][1] >= (p + 1)//2):   # even if the seed of random_point is fixed, the sign of point is random.
        P = -P
    return P

# return a basis of E[l^e], where l is prime.
def basis(E, l, e):
    P = point_ord(E, l, e)
    Q = point_ord(E, l, e)
    while (l**(e-1)*P).weil_pairing(l**(e-1)*Q, l) == 1:
        Q = point_ord(E, l, e)
    return P, Q

p = 2^247 * 79 - 1
e = 247
Fp4, Fp2, zeta2 = calcFields(p)

E0 = EllipticCurve(Fp2, [1, 0])
P, Q = basis(E0, 2, e)
Ms = end.action_matrices([P, Q], 2^e, zeta2, Fp4)
print([det(M) % 2^e for M in Ms])