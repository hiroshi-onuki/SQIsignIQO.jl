export ec_dlog_power_of_2, make_dlog_table

# return n1, n2, n3, n4 such that P = [n1]R + [n1]S, Q = [n3]R + [n4]S
function ec_dlog_power_of_2(P::Point{T}, Q::Point{T}, R::Point{T}, S::Point{T}, 
    A::T, e::Integer) where T <: RingElem
    w0 = Weil_pairing_2power(A, R, S, e)
    w1 = Weil_pairing_2power(A, P, S, e)
    w2 = Weil_pairing_2power(A, P, R, e)
    w3 = Weil_pairing_2power(A, Q, S, e)
    w4 = Weil_pairing_2power(A, Q, R, e)

    n1 = fq_dlog_power_of_2(w1, w0, e)
    n2 = -fq_dlog_power_of_2(w2, w0, e)
    n3 = fq_dlog_power_of_2(w3, w0, e)
    n4 = -fq_dlog_power_of_2(w4, w0, e)

    return n1, n2, n3, n4
end

# return n1, n2, n3, n4 such that P = [n1]R + [n1]S, Q = [n3]R + [n4]S, where xPQ = x(P-Q)
function ec_dlog_power_of_2(xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, R::Point{T}, S::Point{T}, 
            A::T, e::Integer) where T <: RingElem
    P = Point(A, xP)
    Q = Point(A, xQ)
    PQ = add(P, -Q, Proj1(A))
    if !(xPQ == Proj1(PQ.X, PQ.Z))
        Q = -Q
    end

    return ec_dlog_power_of_2(P, Q, R, S, A, e)
end

# return n such that x = base^e
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

# make a precomputed table for dlog with base of order 2^e
function make_dlog_table(base::FqFieldElem, e::Int, window_size::Int)
    F = parent(base)
    l = 2^window_size
    f, r = divrem(e, window_size)
    T1 = [[F(1) for _ in 1:l] for _ in 1:(f+1)]
    T2 = [[F(1) for _ in 1:l] for _ in 1:f]

    T1[1][2] = 1/base
    for j in 2:l-1
        T1[1][j+1] = (T1[1][2])^j
    end
    for i in 1:f
        for j in 1:l-1
            T1[i+1][j+1] = (T1[i][j+1])^l
        end
    end

    for i in 1:f
        for j in 2:l
            T2[i][j] = T1[i][j]^(2^r)
        end
    end

    return T1, T2
end