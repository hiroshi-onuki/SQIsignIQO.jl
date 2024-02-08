# element in the maximal order <1, i, (i + j)/2, (1 + ij)/2> in B_{p, infinity}
# represented by the coefficients of the above basis
struct QOrderElem
    a::BigInt
    b::BigInt
    c::BigInt
    d::BigInt
    p::BigInt
    nj::BigInt
end

function QOrderElem(a::Integer, b::Integer, c::Integer, d::Integer, p::Integer)
    return QOrderElem(a, b, c, d, p, div(p + 1, 4))
end

function QOrderElem(a::Integer, p::Integer)
    return QOrderElem(a, 0, 0, 0, p, div(p + 1, 4))
end

function Base.getindex(x::QOrderElem, i::Integer)
    if i == 1
        return x.a
    elseif i == 2
        return x.b
    elseif i == 3
        return x.c
    elseif i == 4
        return x.d
    else
        throw(BoundsError(x, i))
    end
end

function Base.:(==)(x::QOrderElem, y::QOrderElem)
    return x.a == y.a && x.b == y.b && x.c == y.c && x.d == y.d && x.p == y.p
end

function Base.:+(x::QOrderElem, y::QOrderElem)
    return QOrderElem(x.p, x.a + y.a, x.b + y.b, x.c + y.c, x.d + y.d)
end

function Base.:-(x::QOrderElem, y::QOrderElem)
    return QOrderElem(x.p, x.a - y.a, x.b - y.b, x.c - y.c, x.d - y.d)
end

function Base.:-(x::QOrderElem)
    return QOrderElem(x.p, -x.a, -x.b, -x.c, -x.d)
end

function left_mult_matrix(x::QOrderElem)
    return [x.a -x.b -x.b-x.nj*x.c -x.nj*x.d
            x.b x.a -x.nj*x.d x.b+x.nj*x.c
            x.c x.d x.a+x.d -x.b
            x.d -x.c x.b x.a+x.d]
end

function Base.:*(x::QOrderElem, y::QOrderElem)
    a, b, c, d = left_mult_matrix(x) * [y.a, y.b, y.c, y.d]
    return QOrderElem(a, b, c, d, x.p, x.nj)
end

function involution(x::QOrderElem)
    return QOrderElem(x.a+x.d, -x.b, -x.c, -x.d, x.p, x.nj)
end

function norm(x::QOrderElem)
    return div((2*x.a + x.d)^2 + (2*x.b + x.c)^2 + x.p*(x.c^2 + x.d^2), 4)
end
