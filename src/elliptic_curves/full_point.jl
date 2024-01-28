struct Point{T <: RingElem}
    X::T
    Y::T
    Z::T
end

function Point(A::T, XZ::Proj1{T}) where T <: RingElem
    X, Z = XZ.X, XZ.Z
    XZ = X*Z
    Z2 = Z^2
    Y, Yd = square_root(XZ*(X^2 + A*XZ + Z2))
    X = XZ * Yd
    Z = Z2 * Yd
    @assert Y^2*Z^2 == X^3*Z + A*X^2*Z^2 + X*Z^3
    return Point(X, Y, Z)
end