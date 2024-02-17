export Weil_pairing_2power

# Miller function f_{P}(Q)
function Miller_function(A::T, P::Point{T}, Q::Point{T}, e::Integer) where T <: RingElem
    F = parent(A)
    X, Y, Z = P.X, P.Y, P.Z
    f1, f2 = F(1), F(1)
    for i in 1:e-1
        AZ = A*Z
        QXZ = Q.X*Z
        QZX = Q.Z*X
        QZY = Q.Z*Y
        QZZ = Q.Z*Z
        lam1 = (3*X + 2*AZ) * X + Z^2
        lam12Z = lam1^2*Z
        lam2 = 2*Y*Z
        lam22 = lam2^2
        h1 = lam22 * (Q.Y*Z - QZY) - lam1*lam2 * (QXZ - QZX)
        h2 = lam22 * (QXZ + 2*QZX + A*QZZ) - lam12Z*Q.Z
        f1 = f1^2 * h1
        f2 = f2^2 * h2
        if i < e-1
            lam23 = lam2^3
            X, Y, Z = lam12Z*lam2 - lam23*(AZ + 2*X),
                lam1 * (lam22 * (3*X + AZ) - lam12Z) - Y * lam23,
                Z * lam23
        else
            X, Z = lam1^2*Z - lam22*(AZ + 2*X), Z * lam22
        end
    end
    f1 = f1^2 * (Q.X * Z - Q.Z * X)
    f2 = f2^2 * Q.Z * Z
    return f1, f2
end

# Weil pairing e_{2^e}(P, Q)
function Weil_pairing_2power(A::T, P::Point{T}, Q::Point{T}, e::Integer) where T <: RingElem
    fPQ1, fPQ2 = Miller_function(A, P, Q, e)
    if fPQ1 == 0 || fPQ2 == 0
        return parent(A)(1)
    end
    fQP1, fQP2 = Miller_function(A, Q, P, e)
    return (fPQ1*fQP2) / (fPQ2*fQP1)
end
