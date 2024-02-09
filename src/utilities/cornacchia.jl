
# floor(sqrt(n))
function int_sqrt(n::T) where T <: Integer
    n == 0 && return T(0)
    x = n + 1
    y = n
    while y < x
        x = y
        y = div(x + div(n, x), 2)
    end
    return x
end

# LegendreSymbol (a|p)
function quadratic_residue_symbol(a::Integer, p::T) where T <: Integer
    r = powermod(a, T((p-1)//2), p)
    (r - 1) % p == 0 ? (return 1) : return -1
end

# square root of a mod p. p is a prime.
function sqrt_mod(a::Integer, p::Integer)
    p == 2 && return 1
    p % 4 == 3 && return powermod(a, div(p+1, 4), p)
    if p % 8 == 5
        x = powermod(a, div(p+3, 8), p)
        (a - x^2) % p == 0 ? (return x) : return powermod(2, div(p-1, 4), p)*x % p
    elseif p % 8 == 1
        d = 2
        while (quadratic_residue_symbol(d, p) + 1) % p != 0
            d += 1
        end
        e = 0
        t = p - 1
        while t % 2 == 0
            e += 1
            t >>= 1
        end
        A = powermod(a, t, p)
        D = powermod(d, t, p)
        m = 0
        for i in 0:e-1
            (powermod(A * powermod(D, m, p), 2^(e-1-i), p) + 1) % p == 0 && (m += 2^i)
        end
        x = powermod(a, div(t+1, 2), p) * powermod(D, div(m, 2), p)
        return x % p
    end
    error("p is not prime")
end

# Return a, b such that a^2 + b^2 = q, where q is 2 or a prime satisfing q = 1 mod 4.
function Cornacchia_Smith(q::Integer)
    q == 2 && return [1, 1]
    x = sqrt_mod(-1, q)
    a = q
    b = x
    c = int_sqrt(q)
    while b > c
        a, b = b, a % b
    end
    return b, int_sqrt(q - b^2)
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
            l % 4 == 3 && return 0, 1, false
            s, t = Cornacchia_Smith(l)
            a, b = a*s - b*t, a*t + b*s
        end
    end
    if n % 4 == 1 && is_prime(n)
        s, t = Cornacchia_Smith(n)
        a, b = a*s - b*t, a*t + b*s
    elseif n > 1
        return 0, 2, false
    end
    return a, b, true
end
