# double of P. Alorithm 4 in DMPR2023
function dbl(tnull::ThetaNullLv2{T}, P::ThetaPtLv2{T}) where T <: RingElem
    x, y, z, w = Hadamard(square(P))
    x2 = x^2
    y2 = tnull.lamds[1]*y^2
    z2 = tnull.lamds[2]*z^2
    w2 = tnull.lamds[3]*w^2
    xd, yd, zd, wd = Hadamard(x2, y2, z2, w2)
    return ThetaPtLv2(xd, tnull.lams[1]*yd, tnull.lams[2]*zd, tnull.lams[3]*wd)
end

# differential addition of P and Q. Algorithm 3 in DMPR2023
function diff_add(tnull::ThetaNullLv2{T}, P::ThetaPtLv2{T}, Q::ThetaPtLv2{T}, PmQ::ThetaPtLv2{T}) where T <: RingElem
    xP, yP, zP, wP = Hadamard(square(P))
    xQ, yQ, zQ, wQ = Hadamard(square(Q))
    xPQ = xP*xQ
    yPQ = tnull.lamds[1]*yP*yQ
    zPQ = tnull.lamds[2]*zP*zQ
    wPQ = tnull.lamds[3]*wP*wQ
    xPQ, yPQ, zPQ, wPQ = Hadamard(xPQ, yPQ, zPQ, wPQ)
    xyPmQ = PmQ[1]*PmQ[2]
    zwPmQ = PmQ[3]*PmQ[4]
    x = xPQ*zwPmQ*PmQ[2]
    y = yPQ*zwPmQ*PmQ[1]
    z = zPQ*xyPmQ*PmQ[4]
    w = wPQ*xyPmQ*PmQ[3]
    return ThetaPtLv2(x, y, z, w)
end

# return [m]P by Montgomey ladder
function ladder(tnull::ThetaNullLv2{T}, m::ZZRingElem, P::ThetaPtLv2{T}) where T <: RingElem
    m == 0 && return ThetaPtLv2([tnull[i] for i in 1:4])
    m == 1 && return P
    m == 2 && return dbl(tnull, P)

    t = m >> 1
    b = ZZ(1)
    while t != 1
        t >>= 1
        b <<= 1 
    end

    P0, P1 = P, dbl(tnull, P)
    while b != 0
        if m & b == 0
            P0, P1 = dbl(tnull, P0), diff_add(tnull, P0, P1, P)
        else
            P1, P0 = dbl(tnull, P1), diff_add(tnull, P0, P1, P)
        end
        b >>= 1
    end
    return P0
end