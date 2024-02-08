# element in the maximal order <1, i, (i + j)/2, (1 + ij)/2> in B_{p, infinity}
# represented by the coefficients of the above basis
struct QOrderElem
    p::ZZRingElem,
    a::ZZRingElem,
    b::ZZRingElem,
    c::ZZRingElem,
    d::ZZRingElem
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

function Base.+()
