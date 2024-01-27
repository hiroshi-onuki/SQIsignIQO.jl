abstract type ThetaLv2{T <: RingElem} end

struct ThetaPtLv2{T <: RingElem} <: ThetaLv2{T}
    a::T
    b::T
    c::T
    d::T
end

struct ThetaNullLv2{T <: RingElem} <: ThetaLv2{T}
    a::T
    b::T
    c::T
    d::T
    lams::Vector{T}
    lamds::Vector{T}
end

function ThetaNullLv2(a::T, b::T, c::T, d::T) where T <: RingElem
    ad, bd, cd, dd = Hadamard(a^2, b^2, c^2, d^2)
    inv_b, inv_c, inv_d, inv_bd, inv_cd, inv_dd = batched_inversion([b, c, d, bd, cd, dd])
    lams = [a*inv_b, a*inv_c, a*inv_d]
    lamds = [ad*inv_bd, ad*inv_cd, ad*inv_dd]
    return ThetaNullLv2(a, b, c, d, lams, lamds)
end

function ThetaPtLv2(v::Vector{T}) where T <: RingElem
    return ThetaPtLv2(v[1], v[2], v[3], v[4])
end

function Base.getindex(Theta::ThetaLv2{T}, i::Integer) where T <: RingElem
    if i == 1
        return Theta.a
    elseif i == 2
        return Theta.b
    elseif i == 3
        return Theta.c
    elseif i == 4
        return Theta.d
    else
        throw(BoundError(Theta, i))
    end
end

function Base.setindex!(Theta::ThetaLv2{T}, val::T, i::Integer) where T <: RingElem
    if i == 1
        Theta.a = val
    elseif i == 2
        Theta.b = val
    elseif i == 3
        Theta.c = val
    elseif i == 4
        Theta.d = val
    else
        throw(BoundError(Theta, i))
    end
end

function Base.:*(th1::ThetaLv2{T}, th2::ThetaLv2{T}) where T <: RingElem
    return [th1[i]*th2[i] for i in 1:4]
end

function Base.:*(c::T, th::ThetaLv2{T}) where T <: RingElem
    return ThetaPtLv2([th[i]*c for i in 1:4])
end

function square(th::ThetaLv2{T}) where T <: RingElem
    return [th[i]^2 for i in 1:4]
end

function Base.:(==)(th1::ThetaLv2{T}, th2::ThetaLv2{T}) where T <: RingElem
    i = 1
    while th1[i] == 0
        i += 1
    end
    th2[i] == 0 && return false
    return th1[1]*th2[i] == th2[1]*th1[i] && th1[2]*th2[i] == th2[2]*th1[i] && th1[3]*th2[i] == th2[3]*th1[i] && th1[4]*th2[i] == th2[4]*th1[i]
end

function Hadamard(a::T, b::T, c::T, d::T) where T <: RingElem
    ad = a + b + c + d
    bd = a - b + c - d
    cd = a + b - c - d
    dd = a - b - c + d
    return [ad, bd, cd, dd]
end

function Hadamard(v::Union{Vector{T}, ThetaLv2{T}}) where T <: RingElem
    return Hadamard(v[1], v[2], v[3], v[4])
end

function product_theta_null(t1::ThetaLv1{T}, t2::ThetaLv1{T}) where T <: RingElem
    return ThetaNullLv2(t1[1]*t2[1], t1[2]*t2[1], t1[1]*t2[2], t1[2]*t2[2])
end

function product_theta_pt(t1::ThetaLv1{T}, t2::ThetaLv1{T}) where T <: RingElem
    return ThetaPtLv2(t1[1]*t2[1], t1[2]*t2[1], t1[1]*t2[2], t1[2]*t2[2])
end
