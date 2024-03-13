# return n1, n2, n3, n4 such that P = [n1]P0 + [n1]Q0, Q = [n3]P0 + [n4]Q0, where (P0, Q0) is a fixed basis of E0[2^ExponentFull]
function ec_bi_dlog_E0(P::Point{T}, Q::Point{T}, cdata::CurveData) where T <: RingElem
    A0 = cdata.A0
    P0, Q0 = cdata.P2e, cdata.Q2e
    e = ExponentFull
    w1 = Weil_pairing_2power(A0, P, Q0, e)
    w2 = Weil_pairing_2power(A0, P, P0, e)
    w3 = Weil_pairing_2power(A0, Q, Q0, e)
    w4 = Weil_pairing_2power(A0, Q, P0, e)

    n1 = fq_dlog_power_of_2_opt(w1, cdata.dlog_data_full)
    n2 = -fq_dlog_power_of_2_opt(w2, cdata.dlog_data_full)
    n3 = fq_dlog_power_of_2_opt(w3, cdata.dlog_data_full)
    n4 = -fq_dlog_power_of_2_opt(w4, cdata.dlog_data_full)

    return n1, n2, n3, n4
end

# return n1, n2, n3, n4 such that P = [n1]P0 + [n1]Q0, Q = [n3]P0 + [n4]Q0, where (P0, Q0) is a fixed basis of E0[2^ExponentFull]
function ec_bi_dlog_E0(xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, cdata::CurveData) where T <: RingElem
    A0 = cdata.A0
    P = Point(A0, xP)
    Q = Point(A0, xQ)
    PQ = add(P, -Q, Proj1(A0))
    if !(xPQ == Proj1(PQ.X, PQ.Z))
        Q = -Q
    end

    return ec_bi_dlog_E0(P, Q, cdata)
end

# x^(2^e)
function square_e(x::FqFieldElem, e::Int)
    y = x
    for i in 1:e
        y = y^2
    end
    return y
end

# return n such that x = base^e
function fq_dlog_power_of_2(x::FqFieldElem, base::FqFieldElem, e::Integer)
    n = BigInt(0)
    t = x
    for i in 1:e
        if t^(BigInt(2)^(e-i)) == base^(BigInt(2)^(e-1))
            n += BigInt(2)^(i-1)
            t //= square_e(base, i-1)
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
            T1[i+1][j+1] = square_e(T1[i][j+1], window_size)
        end
    end

    for i in 1:f
        for j in 2:l
            T2[i][j] = square_e(T1[i][j], r)
        end
    end

    return T1, T2
end

# compute x=[x0,x1,...,x{k-1}] s.t. h = (g^(2^r*l^(e-k)))^x
function fq_dlog_subtree(e::Int, h::FqFieldElem, window_size::Int,
        strategy::Vector{Int}, table::Vector{Vector{FqFieldElem}})
    t = length(strategy)
    l = BigInt(2^window_size)
    if t == 0
        h == 1 && return [0]
        for j in 1:l-1
            if h == table[end][j+1]
                return [l - j]
            end
        end 
    end
    n = strategy[1]
    L = strategy[2:t-n+1]
    R = strategy[t-n+2:t]

    hL = square_e(h, n * window_size)
    xL = fq_dlog_subtree(e - n, hL, window_size, L, table)

    hR = h
    for i in 1:e-n
        hR *= table[end-e+i][xL[i]+1]
    end
    xR = fq_dlog_subtree(n, hR, window_size, R, table)

    return vcat(xL, xR)
end

# return n such that h = g^n, where g is a fixed base of order 2^e
function fq_dlog_power_of_2_opt(h::FqFieldElem, dlog_data::DlogData)
    e, window_size, T1, T2, strategy = dlog_data.e, dlog_data.window_size, dlog_data.T1, dlog_data.T2, dlog_data.strategy
    l = BigInt(2^window_size)
    f, r = divrem(e, window_size)
    xw = fq_dlog_subtree(f, square_e(h, r), window_size, strategy, T2)

    hr = h
    for i in 1:f
        hr *= T1[i][xw[i] + 1]
    end
    xr = 0
    for j in 0:2^r-1
        if hr == T1[end][2^r-j+1]
            xr = j
            break
        end
    end

    x = BigInt(0)
    for i in 0:f-1
        x += xw[i+1] * l^i
    end
    x += xr * l^f
    return x
end