export square_root

# square root of x in Fp2
function square_root(x::FinFieldElem)
    Fp2 = parent(x)
    Fp = base_field(Fp2)
    i = gen(Fp2)
    p = characteristic(Fp)
    inv2 = div(p + 1, 2)
    a, b = Fp(coeff(x, 0)), Fp(coeff(x, 1))
    if b == 0
        if a^div(p - 1, 2) == 1
            return [a^div(p + 1, 4), 1]
        else
            return [(-a)^div(p + 1, 4)*i, 1]
        end
    end
    d = (a^2 + b^2)^div(p + 1, 4)
    x = ((a + d)/2)^div(p + 1, 4)
    if x^2 != (a + d) * inv2
        x = ((a - d) * inv2)^div(p + 1, 4)
    end
    y = b*inv2
    return [x^2 + y*i, x]
end
