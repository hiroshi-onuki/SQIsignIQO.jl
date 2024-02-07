
struct ThetaLv1{T <: RingElem}
    a::T
    b::T
end

function Base.getindex(Theta::ThetaLv1{T}, i::Integer) where T <: RingElem
    if i == 1
        return Theta.a
    elseif i == 2
        return Theta.b
    else
        throw(BoundError(Theta, i))
    end
end

# Montgomery coefficient to theta null
function Montgomery_to_theta(A::T) where T <: RingElem
    d = sqrt(A^2 - 4)
    a = (-A + d)/2
    return ThetaLv1(sqrt((a + 1)*(a - 1)), a - 1)
end

# theta null to Montgomery coefficient
function theta_to_Montgomery(tnull::ThetaLv1{T}) where T <: RingElem
    a, b = tnull.a, tnull.b
    a2, b2 = a^2, b^2
    T1, T2 = a2 + b2, a2 - b2
    A = -(T1^2 + T2^2) / (T1*T2)
    return A
end

# Montgomery x-coordinate to theta point
function Montgomery_point_to_theta(tnull::ThetaLv1{T}, P::Proj1{T}) where T <: RingElem
    return ThetaLv1(tnull[1]*(P.X - P.Z), tnull[2]*(P.X + P.Z))
end

# theta point to Montgomery x-coordinate
function theta_point_to_Montgomery(tnull::ThetaLv1{T}, P::ThetaLv1{T}) where T <: RingElem
    a, b = tnull.a, tnull.b
    U, V = P.a, P.b
    aV = a*V
    bU = b*U
    return Proj1(aV + bU, aV - bU)
end