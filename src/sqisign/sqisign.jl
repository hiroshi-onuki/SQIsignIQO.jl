# Sample a random ideal of prime norm D
function sample_random_ideal(D::Integer)
    @assert D % 4 == 3
    gamma, found = FullRepresentInteger(D * BigInt(2)^Log2p)
    !found && throw(ArgumentError("Could not find a random ideal"))
    a = rand(1:D-1)
    println(D)
    println(a)
    return LeftIdeal(gamma * (QOrderElem(a)  + Quoternion_i), D)
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
        D_sec = random_secret_prime()
        I_sec = sample_random_ideal(D_sec)
        alpha, found = KeyGenKLPT(I_sec, D_sec)
        !found && continue
        J = ideal_transform(I_sec, alpha, D_sec)
        g = gcd(J)
        d = 2*Int(log(2, g))
        J = div(J, g)
        println(factor(ZZ(norm(J))))
        a24, xP, xQ, xPQ, M = ideal_to_isogeny_from_O0(J, KLPT_keygen_length - d, cdata, nothing)
        pk = Montgomery_coeff(a24)
        sk = (xP, xQ, xPQ, M, I_sec)
    end
    return pk, sk
end