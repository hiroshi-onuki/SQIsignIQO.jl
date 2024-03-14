using SHA

# Sample a random ideal of prime norm D
function sample_random_ideal(D::Integer)
    @assert D % 4 == 3
    gamma, found = FullRepresentInteger(D * BigInt(2)^Log2p)
    !found && throw(ArgumentError("Could not find a random ideal"))
    a = rand(1:D-1)
    return LeftIdeal(gamma * (QOrderElem(a)  + Quaternion_i), D)
end

# Sample a random ideal of prime norm 2^e
function sample_random_ideal_2e(e::Int)
    gamma = Quaternion_1
    while norm(gamma) % BigInt(2)^e != 0
        gamma, found = FullRepresentInteger(BigInt(2)^(Log2p + e))
        !found && continue
        gamma = div(gamma, gcd(gamma))
        if gcd(gamma * (Quaternion_1 - Quaternion_i)) % 2 == 0
            gamma = div(gamma * (Quaternion_1 - Quaternion_i), 2)
        end
    end
    I = LeftIdeal(gamma, BigInt(2)^e)
    a = rand(1:BigInt(2)^(e))
    return pushforward((1 + a) * Quaternion_1 + a * Quaternion_j, I)
end

# return a random prime <= 2^KLPT_secret_key_prime_size and = 3 mod 4
function random_secret_prime()
    B = BigInt(2)^((KLPT_secret_key_prime_size) - 2) - 1
    n = rand(1:B)
    while !is_prime(4*n + 3)
        n = rand(1:B)
    end
    return 4*n + 3
end

function key_gen(cdata::CurveData)
    found = false
    counter = 0
    pk, sk = nothing, nothing
    while !found && counter < SQISIGN_keygen_attempts
        counter += 1

        # compute a secret ideal
        D_sec = random_secret_prime()
        I_sec = sample_random_ideal(D_sec)
        alpha, found = KeyGenKLPT(I_sec, D_sec)
        !found && continue
        g = gcd(alpha)
        t = g
        d = 0
        while t & 1 == 0
            t >>= 1
            d += 1
        end
        d *= 2
        alpha = div(alpha, g)
        J = ideal_transform(I_sec, alpha, D_sec)
        alpha = involution(alpha) # alpha in J

        # find m s.t. m^2 * D_sec is 2^ExponentForTorsion-good
        m = -1
        a, b = 0, 0
        found = false
        while !found
            m += 2
            a, b, found = sum_of_two_squares(BigInt(2)^ExponentForTorsion - m^2 * D_sec)
        end
        alpha = m * alpha

        # ideal to isogeny
        a24 = cdata.a24_0
        xP, xQ, xPQ = cdata.xP2e, cdata.xQ2e, cdata.xPQ2e
        M = BigInt[1 0; 0 1]
        is_first = true
        extdeg = ExtraDegree
        D = 1
        e = KLPT_keygen_length - d
        while e > ExponentForIsogeny
            n_I_d = D * extdeg * BigInt(2)^ExponentForIsogeny
            I_d = larger_ideal(J, n_I_d)
            a24, xP, xQ, xPQ, M, beta, D, found = short_ideal_to_isogeny(I_d, a24, xP, xQ, xPQ, M, D, ExponentForIsogeny, cdata, is_first, Quaternion_0, 0, 0)
            !found && break
            J = ideal_transform(J, beta, n_I_d)
            alpha = div(alpha * involution(beta), n_I_d)
            e -= ExponentForIsogeny
            is_first = false
            extdeg = 1
        end
        if found
            a24, xP, xQ, xPQ, M, beta, D, found = short_ideal_to_isogeny(J, a24, xP, xQ, xPQ, M, D, e, cdata, false, alpha, a, b)
        else
            println("Failed in e = $e")
            continue
        end
        !found && continue
        pk = Montgomery_coeff(a24)
        M = (M * m) .% BigInt(2)^ExponentFull   # M corresponds to (m * phi_I_sec^-1)^-1, so we need to multiply m
        sk = (xP, xQ, xPQ, M, I_sec)
    end
    return pk, sk, found
end

function commitment(cdata::CurveData)
    # make a random ideal of norm ExtraDegree * 2^SQISIGN_commitment_length
    I2 = sample_random_ideal_2e(SQISIGN_commitment_length)
    gamma, found = FullRepresentInteger(ExtraDegree * BigInt(2)^Log2p)
    !found && throw(ArgumentError("Could not find an ideal of norm ExtraDegree * 2^Log2p"))
    Iex = LeftIdeal(gamma, ExtraDegree)
    I = intersection(I2, Iex)

    # ideal to isogeny
    a24 = cdata.a24_0
    xP, xQ, xPQ = cdata.xP2e, cdata.xQ2e, cdata.xPQ2e
    M = BigInt[1 0; 0 1]
    found = false
    is_first = true
    extdeg = ExtraDegree
    D = 1
    e = SQISIGN_commitment_length
    while e > 0
        ed = min(e, ExponentForIsogeny)
        n_I_d = D * extdeg * BigInt(2)^ed
        I_d = larger_ideal(I, n_I_d)
        a24, xP, xQ, xPQ, M, beta, D, found = short_ideal_to_isogeny(I_d, a24, xP, xQ, xPQ, M, D, ed, cdata, is_first, Quaternion_0, 0, 0)
        !found && break
        I = ideal_transform(I, beta, n_I_d)
        e -= ed
        is_first = false
        extdeg = 1
    end
   return Montgomery_coeff(a24), xP, xQ, xPQ, M, I, found
end

# challenge is the isogeny with kernel <P + [c]Q> from a commitment curve E_com,
# where (P, Q) is a basis of E_com[2^SQISIGN_challenge_length] determined by the fixed torsion basis
function challenge(A::FqFieldElem, xP::Proj1{FqFieldElem}, xQ::Proj1{FqFieldElem}, xPQ::Proj1{FqFieldElem}, m::String, cdata::CurveData)
    h = sha3_256(string(A) * m)

    c = BigInt(0)
    for i in 1:div(SQISIGN_challenge_length,8)
        c += BigInt(h[i]) << (8*(i-1))
    end

    a24 = A_to_a24(A)
    xP = xDBLe(xP, a24, ExponentFull - SQISIGN_challenge_length)
    xQ = xDBLe(xQ, a24, ExponentFull - SQISIGN_challenge_length)
    xPQ = xDBLe(xPQ, a24, ExponentFull - SQISIGN_challenge_length)
    ker = ladder3pt(c, xP, xQ, xPQ, a24)

    a24d, im = two_e_iso(a24, ker, SQISIGN_challenge_length, [xQ])
    a24d, im = Montgomery_normalize(a24d, im)
    Ad = Montgomery_coeff(a24d)
    xK = im[1]
    xPd, xQd, xPQd = torsion_basis(a24d, SQISIGN_challenge_length)
    Pd = Point(Ad, xPd)
    Qd = Point(Ad, xQd)
    PQ = add(Pd, -Qd, Proj1(Ad))
    if !(xPQd == Proj1(PQ.X, PQ.Z))
        Qd = -Qd
    end
    s1, s2 = ec_bi_dlog_challenge(Ad, xK, Pd, Qd, cdata)
    @assert xK == linear_comb_2_e(s1, s2, xPd, xQd, xPQd, a24d, SQISIGN_challenge_length)

    ker_d = s1 % 2 == 0 ? xPd : xQd
    @assert !is_infinity(xDBLe(ker_d, a24d, SQISIGN_challenge_length - 1))
    @assert is_infinity(xDBLe(ker_d, a24d, SQISIGN_challenge_length))
    a24dd, im = two_e_iso(a24d, xK, SQISIGN_challenge_length, [ker_d])
    a24dd, im = Montgomery_normalize(a24dd, im)
    @assert !is_infinity(xDBLe(im[1], a24, SQISIGN_challenge_length - 1))
    @assert a24dd == a24
    r = ec_dlog(A, ker, im[1], xQ, cdata)
    @assert ladder(r, im[1], a24) == ker

    return c, s1, s2, r
end

function signing(pk::FqFieldElem, sk, m::String, cdata::CurveData)
    xP_A, xQ_A, xPQ_A, M_A, I_A = sk

    Acom, xP, xQ, xPQ, Mcom, Icom, found = commitment(cdata)
    !found && return nothing, nothing, false
    cha, s1, s2, r = challenge(Acom, xP, xQ, xPQ, m, cdata)

    # pull-back of the challenge ideal
    Mcom_inv = [Mcom[2,2] -Mcom[1,2]; -Mcom[2,1] Mcom[1,1]] * invmod(Mcom[1, 1] * Mcom[2, 2] - Mcom[1, 2] * Mcom[2, 1], BigInt(2)^ExponentFull)
    a, b = Mcom_inv * [1, cha]
    a, b, c, d = cdata.Matrix_2ed_inv * [b, 0, -a, 0]
    alpha = QOrderElem(a, b, c, d)
    Icha = LeftIdeal(alpha, BigInt(2)^SQISIGN_challenge_length)

    # make a left ideal I of norm I_A * 2^KLPT_signing_klpt_length
    I = intersection(Icom, Icha)
    I, found = SigningKLPT(I_A, I, norm(I_A), norm(I))
    !found && return nothing, nothing, false
    I = intersection(I_A, I)
    @assert gcd(I) == 1
    @assert norm(I) == norm(I_A) * BigInt(2)^KLPT_signing_klpt_length

    # ideal to isogeny
    a24 = A_to_a24(pk)
    xP, xQ, xPQ = xP_A, xQ_A, xPQ_A
    M = M_A
    D = norm(I_A)

    sign = Vector{BigInt}[]
    e = KLPT_signing_klpt_length
    compute_coeff = true
    while e > 0
        # compute the kernel coefficients for signature once every two times
        if compute_coeff
            ed2 = min(e, 2*ExponentForIsogeny)
            a, b = kernel_coefficients(I, M, 2, ed2, cdata.Matrices_2e)
            push!(sign, [a, b])
            compute_coeff = false
        else
            compute_coeff = true
        end
        
        ed = min(ExponentForIsogeny, e)
        n_I_d = D * BigInt(2)^ed
        I_d = larger_ideal(I, n_I_d)
        a24, xP, xQ, xPQ, M, beta, D, found = short_ideal_to_isogeny(I_d, a24, xP, xQ, xPQ, M, D, ed, cdata, false, Quaternion_0, 0, 0)
        !found && break
        I = ideal_transform(I, beta, n_I_d)
        e -= ed
    end

    # challenge ellitpic curve
    a24d = A_to_a24(pk)
    e = KLPT_signing_klpt_length
    for (a, b) in sign
        ed = min(2*ExponentForIsogeny, e)
        xP, xQ, xPQ = torsion_basis(a24d, ExponentFull)
        xP = xDBLe(xP, a24d, ExponentFull - ed)
        xQ = xDBLe(xQ, a24d, ExponentFull - ed)
        xPQ = xDBLe(xPQ, a24d, ExponentFull - ed)
        ker = linear_comb_2_e(a, b, xP, xQ, xPQ, a24d, ed)
        a24d, _ = two_e_iso(a24d, ker, ed, Proj1{FqFieldElem}[])
        a24d, _ = Montgomery_normalize(a24d, Proj1{FqFieldElem}[])
        e -= ed
    end
    @assert a24d == a24

    return sign, s1, s2, r, true
end

function verify(pk::FqFieldElem, m::String, sign::Vector{Vector{BigInt}}, s1::BigInt, s2::BigInt, r::BigInt)
    # challenge ellitpic curve
    a24 = A_to_a24(pk)
    e = KLPT_signing_klpt_length
    for (a, b) in sign
        ed = min(2*ExponentForIsogeny, e)
        xP, xQ, xPQ = torsion_basis(a24, ExponentFull)
        xP = xDBLe(xP, a24, ExponentFull - ed)
        xQ = xDBLe(xQ, a24, ExponentFull - ed)
        xPQ = xDBLe(xPQ, a24, ExponentFull - ed)
        ker = linear_comb_2_e(a, b, xP, xQ, xPQ, a24, ed)
        a24, _ = two_e_iso(a24, ker, ed, Proj1{FqFieldElem}[])
        a24, _ = Montgomery_normalize(a24, Proj1{FqFieldElem}[])
        e -= ed
    end

    # commitment elliptic curve
    xP, xQ, xPQ = torsion_basis(a24, SQISIGN_challenge_length)
    ker = linear_comb_2_e(s1, s2, xP, xQ, xPQ, a24, SQISIGN_challenge_length)
    xR = s1 % 2 == 0 ? xP : xQ
    a24com, im = two_e_iso(a24, ker, SQISIGN_challenge_length, [xR])
    a24com, im = Montgomery_normalize(a24com, im)
    Acom = Montgomery_coeff(a24com)
    xP, xQ, xPQ = torsion_basis(a24com, ExponentFull)
    xP = xDBLe(xP, a24com, ExponentFull - SQISIGN_challenge_length)
    xQ = xDBLe(xQ, a24com, ExponentFull - SQISIGN_challenge_length)
    xPQ = xDBLe(xPQ, a24com, ExponentFull - SQISIGN_challenge_length)

    # recover challenge
    h = sha3_256(string(Acom) * m)
    c = BigInt(0)
    for i in 1:div(SQISIGN_challenge_length,8)
        c += BigInt(h[i]) << (8*(i-1))
    end
    xK = ladder3pt(c, xP, xQ, xPQ, a24com)

    return xK == ladder(r, im[1], a24com)
end