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
    while !is_probable_prime(4*n + 3)
        n = rand(1:B)
    end
    return 4*n + 3
end

function key_gen(global_data::GlobalData)
    found = false
    counter = 0
    pk, sk = nothing, nothing
    E0 = global_data.E0
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

        # find m s.t. m^2 * D_sec is 2^ExponentForIsogenyDim2-good
        m = -1
        a, b = 0, 0
        found = false
        while !found
            m += 2
            a, b, found = sum_of_two_squares(BigInt(2)^ExponentForIsogenyDim2 - m^2 * D_sec)
        end
        alpha = m * alpha

        # ideal to isogeny
        a24 = E0.a24_0
        xP, xQ, xPQ = E0.xP2e, E0.xQ2e, E0.xPQ2e
        M = BigInt[1 0; 0 1]
        D = 1
        e = KLPT_keygen_length - d
        while e > ExponentForIsogenyDim1
            n_I_d = D * BigInt(2)^ExponentForIsogenyDim1
            I_d = larger_ideal(J, n_I_d)
            a24, xP, xQ, xPQ, M, beta, D = short_ideal_to_isogeny(I_d, a24, xP, xQ, xPQ, M, D, ExponentForIsogenyDim1, global_data, false, Quaternion_0, 0, 0)
            J = ideal_transform(J, beta, n_I_d)
            alpha = div(alpha * involution(beta), n_I_d)
            e -= ExponentForIsogenyDim1
        end
        a24, xP, xQ, xPQ, M, beta, D = short_ideal_to_isogeny(J, a24, xP, xQ, xPQ, M, D, e, global_data, true, alpha, a, b)
        pk = Montgomery_coeff(a24)
        M = (M * invmod(m, BigInt(2)^ExponentFull)) .% BigInt(2)^ExponentFull   # M corresponds to m * phi_I_sec, so we need to multiply m^-1
        sk = (xP, xQ, xPQ, M, I_sec)
    end
    return pk, sk, found
end

function commitment(global_data::GlobalData)
    # make a random ideal of norm 2^SQISIGN_commitment_length
    I = sample_random_ideal_2e(SQISIGN_commitment_length)

    # ideal to isogeny
    E0 = global_data.E0
    a24 = E0.a24_0
    xP, xQ, xPQ = E0.xP2e, E0.xQ2e, E0.xPQ2e
    M = BigInt[1 0; 0 1]
    D = 1
    e = SQISIGN_commitment_length
    while e > 0
        ed = min(e, ExponentForIsogenyDim1)
        is_normalized = e <= ExponentForIsogenyDim1
        n_I_d = D * BigInt(2)^ed
        I_d = larger_ideal(I, n_I_d)
        a24, xP, xQ, xPQ, M, beta, D = short_ideal_to_isogeny(I_d, a24, xP, xQ, xPQ, M, D, ed, global_data, is_normalized, Quaternion_0, 0, 0)
        I = ideal_transform(I, beta, n_I_d)
        e -= ed
    end
    return Montgomery_coeff(a24), xP, xQ, xPQ, M, I
end

# challenge is the isogeny with kernel <P + [c]Q> from a commitment curve E_com,
# where (P, Q) is a basis of E_com[2^SQISIGN_challenge_length] determined by the fixed torsion basis
function challenge(A::FqFieldElem, xP::Proj1{FqFieldElem}, xQ::Proj1{FqFieldElem}, xPQ::Proj1{FqFieldElem}, m::String, E0::E0Data)
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

    a24d, im = two_e_iso(a24, ker, SQISIGN_challenge_length, [xQ], StrategiesDim1[SQISIGN_challenge_length])
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
    s1, s2 = ec_bi_dlog_challenge(Ad, xK, Pd, Qd, E0)
    @assert xK == linear_comb_2_e(s1, s2, xPd, xQd, xPQd, a24d, SQISIGN_challenge_length)
    is_one_P = s1 % 2 != 0
    if is_one_P
        s = (s2 * invmod(s1, BigInt(2)^SQISIGN_challenge_length)) % BigInt(2)^SQISIGN_challenge_length
    else
        s = (s1 * invmod(s2, BigInt(2)^SQISIGN_challenge_length)) % BigInt(2)^SQISIGN_challenge_length
    end

    ker_d = is_one_P ? xQd : xPd
    @assert !is_infinity(xDBLe(ker_d, a24d, SQISIGN_challenge_length - 1))
    @assert is_infinity(xDBLe(ker_d, a24d, SQISIGN_challenge_length))
    a24dd, im = two_e_iso(a24d, xK, SQISIGN_challenge_length, [ker_d], StrategiesDim1[SQISIGN_challenge_length])
    a24dd, im = Montgomery_normalize(a24dd, im)
    @assert !is_infinity(xDBLe(im[1], a24, SQISIGN_challenge_length - 1))
    @assert a24dd == a24
    r = ec_dlog(A, ker, im[1], xQ, E0)
    @assert ladder(r, im[1], a24) == ker

    return c, is_one_P, s, r
end

function signing(pk::FqFieldElem, sk, m::String, global_data::GlobalData)
    xP_A, xQ_A, xPQ_A, M_A, I_A = sk

    E0 = global_data.E0
    while true
        Acom, xP, xQ, xPQ, Mcom, Icom = commitment(global_data)
        cha, is_one_P, s, r = challenge(Acom, xP, xQ, xPQ, m, E0)

        # pull-back of the challenge ideal
        Mcom_inv = [Mcom[2,2] -Mcom[1,2]; -Mcom[2,1] Mcom[1,1]] * invmod(Mcom[1, 1] * Mcom[2, 2] - Mcom[1, 2] * Mcom[2, 1], BigInt(2)^ExponentFull)
        a, b = Mcom_inv * [1, cha]
        a, b, c, d = E0.Matrix_2ed_inv * [b, 0, -a, 0]
        alpha = QOrderElem(a, b, c, d)
        Icha = LeftIdeal(alpha, BigInt(2)^SQISIGN_challenge_length)

        # make a left ideal I of norm I_A * 2^KLPT_signing_klpt_length
        Icc = intersection(Icom, Icha)
        I, found = SigningKLPT(I_A, Icc, norm(I_A), norm(Icc))
        !found && continue
        @assert norm(I) == BigInt(2)^KLPT_signing_klpt_length
        I = intersection(I_A, I)
        @assert gcd(I) == 1
        @assert norm(I) == norm(I_A) * BigInt(2)^KLPT_signing_klpt_length

        # ideal to isogeny
        a24 = A_to_a24(pk)
        xP, xQ, xPQ = xP_A, xQ_A, xPQ_A
        M = M_A
        D = norm(I_A)

        sign = Vector{UInt8}(undef, SQISIGN_sign_bytes)
        e = KLPT_signing_klpt_length
        compute_coeff = true
        idx = 1
        while e > 0
            # compute the kernel coefficients for signature once every two times
            if compute_coeff
                ed2 = min(e, 2*ExponentForIsogenyDim1)
                a, b = kernel_coefficients(I, M, 2, ed2, E0.Matrices_2e)
                if a == 1
                    sign[idx] = 0x01
                    sign[idx+1:idx+SQISIGN_sign_isogeny_bytes] = integer_to_bytes(b, SQISIGN_sign_isogeny_bytes)
                else
                    sign[idx] = 0x00
                    sign[idx+1:idx+SQISIGN_sign_isogeny_bytes] = integer_to_bytes(a, SQISIGN_sign_isogeny_bytes)
                end
                compute_coeff = false
                idx += 1 + SQISIGN_sign_isogeny_bytes
            else
                compute_coeff = true
            end

            # we do not need to compute the last 2 isogenies
            if e > 2 * ExponentForIsogenyDim1
                ed = min(ExponentForIsogenyDim1, e)
                n_I_d = D * BigInt(2)^ed
                I_d = larger_ideal(I, n_I_d)
                a24, xP, xQ, xPQ, M, beta, D = short_ideal_to_isogeny(I_d, a24, xP, xQ, xPQ, M, D, ed, global_data, compute_coeff, Quaternion_0, 0, 0)
                I = ideal_transform(I, beta, n_I_d)
            end
            e -= ExponentForIsogenyDim1
        end

        if is_one_P
            sign[idx] = 0x01
            sign[idx+1:idx+SQISIGN_challenge_bytes] = integer_to_bytes(s, SQISIGN_challenge_bytes)
        else
            sign[idx] = 0x00
            sign[idx+1:idx+SQISIGN_challenge_bytes] = integer_to_bytes(s, SQISIGN_challenge_bytes)
        end
        idx += 1 + SQISIGN_challenge_bytes
        sign[idx:idx-1+SQISIGN_challenge_bytes] = integer_to_bytes(r, SQISIGN_challenge_bytes)

        found && return sign
    end
end

function verify(pk::FqFieldElem, m::String, sign::Vector{UInt8})
    # decompress signature
    sign_coeffs = []
    idx = 1
    for i in 0:SQISIGN_signing_length-1
        bit = sign[idx]
        coeff = bytes_to_integer(sign[idx + 1:idx + SQISIGN_sign_isogeny_bytes])
        push!(sign_coeffs, (bit, coeff))
        idx += 1 + SQISIGN_sign_isogeny_bytes
    end
    bit_s = sign[idx]
    s = bytes_to_integer(sign[idx + 1:idx + SQISIGN_challenge_bytes])
    r = bytes_to_integer(sign[idx + 1 + SQISIGN_challenge_bytes:idx + 2*SQISIGN_challenge_bytes])

    # challenge ellitpic curve
    a24 = A_to_a24(pk)
    a24_d = nothing
    e = KLPT_signing_klpt_length
    for i in 1:SQISIGN_signing_length
        bit, a = sign_coeffs[i]
        ed = min(2*ExponentForIsogenyDim1, e)
        xP, xQ, xPQ = torsion_basis(a24, ExponentFull)
        xP = xDBLe(xP, a24, ExponentFull - ed)
        xQ = xDBLe(xQ, a24, ExponentFull - ed)
        xPQ = xDBLe(xPQ, a24, ExponentFull - ed)
        if bit == 1
            ker = ladder3pt(a, xP, xQ, xPQ, a24)
        else
            ker = ladder3pt(a, xQ, xP, xPQ, a24)
        end
        if i < SQISIGN_signing_length
            a24, _ = two_e_iso(a24, ker, ed, Proj1{FqFieldElem}[], StrategiesDim1[ed])
        else
            a24, _, a24_d = two_e_iso(a24, ker, ed, Proj1{FqFieldElem}[], StrategiesDim1[ed], -1)
        end
        a24, _ = Montgomery_normalize(a24, Proj1{FqFieldElem}[])
        e -= ed
    end

    # commitment elliptic curve
    xP, xQ, xPQ = torsion_basis(a24, SQISIGN_challenge_length)
    if bit_s == 1
        ker = ladder3pt(s, xP, xQ, xPQ, a24)
        xR = xQ
    else
        ker = ladder3pt(s, xQ, xP, xPQ, a24)
        xR = xP
    end
    a24com, im, a24_dd = two_e_iso(a24, ker, SQISIGN_challenge_length, [xR], StrategiesDim1[SQISIGN_challenge_length], 1)
    a24com, im = Montgomery_normalize(a24com, im)
    Acom = Montgomery_coeff(a24com)
    xP, xQ, xPQ = torsion_basis(a24com, ExponentFull)
    xP = xDBLe(xP, a24com, ExponentFull - SQISIGN_challenge_length)
    xQ = xDBLe(xQ, a24com, ExponentFull - SQISIGN_challenge_length)
    xPQ = xDBLe(xPQ, a24com, ExponentFull - SQISIGN_challenge_length)

    # check cyclicity
    j1 = jInvariant_a24(a24_d)
    j2 = jInvariant_a24(a24_dd)
    j1 == j2 && return false

    # recover challenge
    h = sha3_256(string(Acom) * m)
    c = BigInt(0)
    for i in 1:div(SQISIGN_challenge_length,8)
        c += BigInt(h[i]) << (8*(i-1))
    end
    xK = ladder3pt(c, xP, xQ, xPQ, a24com)

    return xK == ladder(r, im[1], a24com)
end