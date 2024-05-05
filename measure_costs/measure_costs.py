from math import log, prod

# constants from the SQIsign
TORSION_ODD_PRIMES = [0x3, 0x17, 0x3b, 0x65, 0x6d, 0xc5, 0x1eb, 0x2e7, 0x779, 0x7, 0xb, 0xd, 0x25, 0x59, 0x61, 0x6b, 0x83, 0x89, 0xdf, 0xef, 0x17f, 0x185, 0x1f3, 0x25f, 0x409, 0x419, 0x4a9, 0x7b5]
TORSION_ODD_POWERS = [0x24, 0x2, 0x2, 0x2, 0x2, 0x2, 0x2, 0x2, 0x2, 0x4, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1]
TORSION_PLUS_ODD_PRIMES = [0x3, 0x17, 0x3b, 0x65, 0x6d, 0xc5, 0x1eb, 0x2e7, 0x779]
TORSION_PLUS_ODD_POWERS = [0x24, 0x2, 0x2, 0x2, 0x2, 0x2, 0x2, 0x2, 0x2]
TORSION_MINUS_ODD_PRIMES = [0x7, 0xb, 0xd, 0x25, 0x59, 0x61, 0x6b, 0x83, 0x89, 0xdf, 0xef, 0x17f, 0x185, 0x1f3, 0x25f, 0x409, 0x419, 0x4a9, 0x7b5]
TORSION_MINUS_ODD_POWERS = [0x4, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1, 0x1]

M = 3
S = 2
p = 2**247*79 - 1
I = 2 + M + 1.5*log(p, 2)   # '(x + iy)^(-1) = (x^2 + y^2)^(p - 2) * (x - iy)
Fpsqrt = 1.5*log(p, 2)      # sqrt(x) = x^((p + 1)/4) for x in Fp

# compute the cost of a strategy with mul. cost cl and eval. cost cr
def strategy_cost(n, cl, cr):
    C = [0 for i in range(n+1)]
    for k in range(2,n+1):
        j = 1
        z = k - 1
        while j < z:
            m = j + (z-j)//2
            w = m + 1
            t1 = C[m] + C[k-m] + (k-m)*cl + m*cr
            t2 = C[w] + C[k-w] + (k-w)*cl + w*cr
            if t1 <= t2:
                z = m
            else:
                j = w
        C[k] = C[j] + C[k-j] + (k-j)*cl + j*cr
    return C[n]

# the cost of xDBLADD, Algorithm 5 in the SIKE document
def DblAdd():
    return 7*M + 4*S

# the cost of the multiplication of a point by n
def Ladder(n):
    return (len(bin(n)) - 3) * DblAdd()

# the cost of computing a square root in Fp^2
# Algorithm 9 in Adj et al., "Square Root Computation over Even Extension Fields"
def Fp2sqrt():
    return 4*M + 2*Fpsqrt

# the cost of computing x(P + Q) or x(P - Q) from x(P) and x(Q)
def DblAddSub():
    return 11*M + 2*S + Fp2sqrt()

# the cost of computing the image of a point
def IsogPoint(l):
    d = (l - 1)//2
    return 4*d*M + 2*S

# the cost of KPS in an l-isogeny
def IsogKPS(l):
    d = (l - 1)//2
    return 4*(d - 1) * M + 2*(d - 1) * S

# lower bound on the cost of computing the codomain and the image of a point by sqrtVelu
# from Eq.(5) in Adj et al., "Karatsuba-based square-root Vélu’s formulas applied to two isogeny-based protocols"
def IsogCostSqrtVelu(l):
    b = int((l - 1)**0.5 / 2)
    c_kps = 24 * b
    c_eval = 2*(3 * b**log(3, 2) +  log(b, 2) * b - 5/3*b + 5/6) + 18*b + 3 * b**log(3, 2) - 2 * b
    return (c_kps + 3 * c_eval) * S

# the cost of doubling in theta coordinates
def DblDim2():
    return 6*M + 8*S

# the cost of the image of a point under a (2,2)-isogeny
def IsogPointDim2():
    return 3*M + 4*S

# the cost of the codomain of a (2,2)-isogeny
def IsogCodDim2():
    return 9*M + 8*S + I

# the cost of doubling for gluing two curves
def DblGluDim2():
    return 8*M + 4*S

# the cost of the codomain of a gluing isogeny
def IsonCodGluDim2():
    return 13*M + 8*S + I

# the cost of the image of a point under a gluing isogeny
def IsogPointGluDim2():
    return 5*M + 8*S + I

def eval_T_cost():
    n = len(TORSION_ODD_PRIMES)
    m_p = len(TORSION_PLUS_ODD_PRIMES)
    m_m = len(TORSION_MINUS_ODD_PRIMES)
    c = 0

    # cost of multiplications
    ppws = [l**e for l, e in zip(TORSION_PLUS_ODD_PRIMES, TORSION_PLUS_ODD_POWERS)]
    ppws.sort(reverse=True)
    k = prod(ppws)
    for le in ppws:
        k //= le
        c += Ladder(k)
    ppws = [l**e for l, e in zip(TORSION_PLUS_ODD_PRIMES, TORSION_PLUS_ODD_POWERS)]
    ppws.sort(reverse=True)
    k = prod(ppws)
    for le in ppws:
        k //= le
        c += Ladder(k)

    # cost of isogenies
    for i in range(n):
        e = TORSION_ODD_POWERS[i]
        l = TORSION_ODD_PRIMES[i]
        if l == 3:
            c += 3*(3*IsogPoint(l) + IsogKPS(l))
        else:
            c += min(e*(3*IsogPoint(l) + IsogKPS(l)), e*IsogCostSqrtVelu(l))
    return c

def eval_2_2_isog_cost():
    f = 139
    c = 0

    # cost of gluing
    c += DblGluDim2() * (f - 1) # doubling for making the kernel generator
    c += 2 * Ladder(2**(f-1) + 1)
    c += 2 * Ladder(2**(f-1))
    c += 4 * DblAddSub()
    c += IsonCodGluDim2()
    c += 6 * IsogPointGluDim2()

    # cost of remaining isogenies
    c += 2*strategy_cost(f-1, DblDim2(), IsogPointDim2()) + (f-1)*(IsogCodDim2() + 4*IsogPointDim2())

    return c

print(eval_T_cost())
print(eval_2_2_isog_cost())

