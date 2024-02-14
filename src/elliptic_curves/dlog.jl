
function ec_dlog_power_of_2(xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, R::Point{T}, S::Point{T}, 
            A::T, e::Integer) where T <: RingElem
    P = Point(A, xP)
    Q = Point(A, xQ)
    if xPQ != add(P, Q, Proj1(A))
        Q = -Q
    end

    w0 = Weil_pairing_2power(A, R, S, e)
    w1 = Weil_pairing_2power(A, P, S, e)
    w2 = Weil_pairing_2power(A, P, R, e)
    w3 = Weil_pairing_2power(A, Q, S, e)
    w4 = Weil_pairing_2power(A, Q, R, e)

    println(w0, " ", w1, " ", w2, " ", w3, " ", w4)

    n1 = fq_dlog_power_of_2(w1, w0, e)
    n2 = -fq_dlog_power_of_2(w2, w0, e)
    n3 = fq_dlog_power_of_2(w3, w0, e)
    n4 = -fq_dlog_power_of_2(w4, w0, e)

    return n1, n2, n3, n4
end

# return n such that x = base^n
function fq_dlog_power_of_2(x::FqFieldElem, base::FqFieldElem, e::Integer)
    n = BigInt(0)
    t = x
    for i in 1:e
        if t^(BigInt(2)^(e-i)) == base^(BigInt(2)^(e-1))
            n += BigInt(2)^(i-1)
            t //= base^(BigInt(2)^(i-1))
        end
    end
    return n
end