# Sample a random ideal of prime norm D
function sample_random_ideal(D::Integer)
    @assert D % 4 == 3
    gamma, found = FullRepresentInteger(D * BigInt(2)^Log2p)
    !found && throw(ArgumentError("Could not find a random ideal"))
    a = rand(1:D-1)
    return LeftIdeal(gamma * (QOrderElem(a)  + Quoternion_i), D)
end

# Sample a random ideal of prime norm 2^e
function sample_random_ideal_2e(e::Int)
    gamma = Quoternion_1
    while norm(gamma) % BigInt(2)^e != 0
        gamma, found = FullRepresentInteger(BigInt(2)^(Log2p + e))
        !found && continue
        gamma = div(gamma, gcd(gamma))
        if gcd(gamma * (Quoternion_1 - Quoternion_i)) % 2 == 0
            gamma = div(gamma * (Quoternion_1 - Quoternion_i), 2)
        end
    end
    @assert gcd(gamma * (Quoternion_1 - Quoternion_i)) % 2 == 1
    I = LeftIdeal(gamma, BigInt(2)^e)
    @assert !is_subset(I, LeftIdeal(Quoternion_1 + Quoternion_i, 2))
    a = rand(1:BigInt(2)^(e))
    J =  pushforward((1 + a) * Quoternion_1 + a * Quoternion_j, I)
    @assert !is_subset(J, LeftIdeal(Quoternion_1 + Quoternion_i, 2))
    return J
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
    while !found && counter < SQISIGN_response_attempts
        counter += 1

        # compute a secret ideal
        D_sec = random_secret_prime()
        println("D_sec: ", log(2, D_sec))
        I_sec = sample_random_ideal(D_sec)
        alpha, found = KeyGenKLPT(I_sec, D_sec)
        !found && continue
        g = gcd(alpha)
        d = 2*Int(log(2, g))
        alpha = div(alpha, g)
        J = ideal_transform(I_sec, alpha, D_sec)
        alpha = involution(alpha) # alpha in J
        println(factor(ZZ(norm(J))))
        @assert norm(J) == ExtraDegree * BigInt(2)^(KLPT_keygen_length - d)

        # find m s.t. m^2 * D_sec is 2^ExponentForTorsion-good
        m = -1
        a, b = 0, 0
        found_2e_good = false
        while !found_2e_good
            m += 2
            a, b, found_2e_good = sum_of_two_squares(BigInt(2)^ExponentForTorsion - m^2 * D_sec)
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
            a24, xP, xQ, xPQ, M, beta, D = short_ideal_to_isogeny(I_d, a24, xP, xQ, xPQ, M, D, ExponentForIsogeny, cdata, is_first, Quoternion_0, 0, 0)
            J = ideal_transform(J, beta, n_I_d)
            alpha = div(alpha * involution(beta), n_I_d)
            e -= ExponentForIsogeny
            is_first = false
            extdeg = 1
        end
        a24, xP, xQ, xPQ, M, beta, D = short_ideal_to_isogeny(J, a24, xP, xQ, xPQ, M, D, e, cdata, false, alpha, a, b)

        pk = Montgomery_coeff(a24)
        sk = (xP, xQ, xPQ, M, I_sec)
    end
    return pk, sk, found
end