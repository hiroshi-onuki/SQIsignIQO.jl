# LegendreSymbol (a|q)
function quadratic_residue_symbol(a::Integer, q::Integer)
    r = powermod(a, div(q-1, 2), q)
    (r - 1) % q == 0 ? (return 1) : return -1
end

# square root of a mod q. q is a prime.
function sqrt_mod(a::Integer, q::Integer)
    q == 2 && return 1
    q % 4 == 3 && return powermod(a, div(q+1, 4), q)
    if q % 8 == 5
        x = powermod(a, div(q+3, 8), q)
        (a - x^2) % q == 0 ? (return x) : return powermod(2, div(q-1, 4), q)*x % q
    elseif q % 8 == 1
        d = 2
        while (quadratic_residue_symbol(d, q) + 1) % q != 0
            d += 1
        end
        e = 0
        t = q - 1
        while t % 2 == 0
            e += 1
            t >>= 1
        end
        A = powermod(a, t, q)
        D = powermod(d, t, q)
        m = 0
        for i in 0:e-1
            (powermod(A * powermod(D, m, p), 2^(e-1-i), q) + 1) % q == 0 && (m += 2^i)
        end
        x = powermod(a, div(t+1, 2), q) * powermod(D, div(m, 2), q)
        return x % q
    end
    error("q is not prime")
end

# Return a, b such that a^2 + b^2 = q, where q is 2 or a prime satisfing q = 1 mod 4.
function Cornacchia_Smith(q::Integer)
    q == 2 && return [1, 1]
    x = sqrt_mod(-1, q)
    a = q
    b = x
    c = integer_square_root(q)
    while b > c
        a, b = b, a % b
    end
    return b, integer_square_root(q - b^2)
end

# Return a, b such that a^2 + b^2 = n or nothing if not found.
function sum_of_two_squares(n::Integer)
    n <= 0 && return 0, 0, false
    n == 1 && return 1, 0, true
    a, b = 1, 0
    for l in SmallPrimes
        e = 0
        while n % l == 0
            n = div(n, l)
            e += 1
        end
        s = l^(div(e, 2))
        a *= s
        b *= s
        if e % 2 == 1
            l % 4 == 3 && return 0, 0, false
            s, t = Cornacchia_Smith(l)
            a, b = a*s - b*t, a*t + b*s
        end
    end
    if n % 4 == 1 && is_prime(n)
        s, t = Cornacchia_Smith(n)
        a, b = a*s - b*t, a*t + b*s
    elseif n > 1
        return 0, 0, false
    end
    return a, b, true
end
