
function ec_dlog_power_of_2(xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, R::Point{T}, S::Point{T}, 
            A::T, e::Integer) where T <: RingElem
    P = Point(A, xP)
    Q = Point(A, xQ)
    if xPQ != add(P, Q, Proj1(A))
        Q = -Q
    end

    w0 = Weil_pairing_2power(A, R, S, e)
    w1 = Weil_pairing_2power(A, P, R, e)
    w2 = Weil_pairing_2power(A, P, S, e)
    w3 = Weil_pairing_2power(A, Q, R, e)
    w4 = Weil_pairing_2power(A, Q, S, e)
    
end

# return n such that x = base^n
function fq_dlog_power_of_2(x::FqFieldElem, base::FqFieldElem, e::Integer)
    n = BigInt(0)
    t = x
    for i in 1:e-1
        if t^(BigInt(2)^(e-i)) == base^(BigInt(2)^(e-i))
            n += BigInt(2)^(i-1)
        end
        t //= base^(BigInt(2)^(i-1))
    end
    return n
end