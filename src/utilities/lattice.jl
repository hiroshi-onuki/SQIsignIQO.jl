# Return a Z-module basis from input generators gens
function get_basis(gens::Vector{Vector{T}}) where T <: Integer
    n = length(gens)
    n == 0 && return gens
    m = length(gens[1])
    prod([length(v) == m for v in gens]) || error("lengths of generatars are different")
    
    basis = deepcopy(gens)
    i = m
    j = n
    while true
        while j < 1 || prod([b[i] == 0 for b in basis[1:j]])
            i -= 1
            i == max(m - n, 0) && return filter(x->x!=zeros(m), basis)
        end
        while j > 1 && !prod([b[i] == 0 for b in basis[1:j-1]])
            k = findfirst(x->x==minimum(abs.([b[i] for b in filter(x->x[i]!=0, basis[1:j])])), abs.([b[i] for b in basis[1:j]]))
            basis[j], basis[k] = basis[k], basis[j]
            basis[j][i] < 0 && (basis[j] *= -1)
            for k in 1:j-1
                basis[k] -= div(basis[k][i], basis[j][i]) * basis[j]
            end
        end
        j -= 1
    end
end

# Algorithm 2.6.7 in H. Cohen, A Course in Computational Algebraic Number Theory
function integral_LLL(basis::Vector{Vector{T}}, quadratic_form::Function) where T <: Integer
    b = deepcopy(basis)

    # input check
    n = length(b)
    n == 0 && return b
    m = length(b[1])
    prod([length(v) == m for v in b]) || error("lengths of generatars are different")
    m < n && error("number of vectors greater than these dimension")

    q = quadratic_form
    k = 2
    kmax = 1
    d = zeros(T, n+1)
    d[1] = 1
    d[2] = q(b[1], b[1])
    H = zeros(T, n, n)
    for i in 1:n H[i,i] = 1 end
    lam = zeros(T, n, n)

    while k <= n
        if k > kmax
            kmax = k
            for j in 1:k
                u = q(b[k], b[j])
                for i in 1:j-1
                    u = div(d[i+1]*u - lam[k, i]*lam[j, i], d[i])
                end
                if j < k
                    lam[k, j] = u
                elseif j == k
                    d[k+1] = u
                end
            end
            d[k+1] == 0 && error("vectors are not linearly independent")
        end
        b, H, lam = REDI(k, k-1, b, d, H, lam)
        while d[k+1]*d[k-1] < 3//4*d[k]^2 - lam[k, k-1]^2
            b, d, H, lam = SWAPI(k, kmax, b, d, H, lam)
            k = max(2, k-1)
            b, H, lam = REDI(k, k-1, b, d, H, lam)
        end
        for l in 1:k-2
            b, H, lam = REDI(k, k-1-l, b, d, H, lam)
        end
        k += 1
    end
    return H
end

function REDI(k::Int, l::Int, b::Vector{Vector{T}}, d::Vector{T}, H::Matrix{T}, lam::Matrix{T}) where T <: Integer
    if abs(2*lam[k, l]) > d[l+1]
        q = T(round(lam[k,l]//d[l+1]))
        H[:, k] -= q*H[:, l]
        b[k] -= q*b[l]
        lam[k,l] -= q*d[l+1]
        for i in 1:l-1
            lam[k,i] -= q*lam[l,i]
        end
    end
    return b, H, lam
end

function SWAPI(k::Int, kmax::Int, b::Vector{Vector{T}}, d::Vector{T}, H::Matrix{T}, lam::Matrix{T}) where T <: Integer
    H[:, k], H[:, k-1] = H[:, k-1], H[:, k]
    b[k], b[k-1] = b[k-1], b[k]
    if k > 2
        for j in 1:k-2
            lam[k, j], lam[k-1, j] = lam[k-1, j], lam[k, j]
        end
    end
    lamd = lam[k,k-1]
    B = T((d[k-1]*d[k+1] + lamd^2) // d[k])
    for i in k+1:kmax
        t = lam[i,k]
        lam[i,k] = T((d[k+1]*lam[i,k-1] - lamd*t) // d[k])
        lam[i,k-1] = T((B*t + lamd*lam[i,k]) // d[k+1])
    end
    d[k] = B
    return b, d, H, lam
end

#=
    return a list of coefficients of short vectors in alattice.
    basis: a basis of the target lattice, M: bilinear matrix,
    C; upper bound of the quadratic forms of output vectors
    Algorithm 2.7.5 in H. Cohen, A Course in Computational Algebraic Number Theory.
=#
function short_vectors(basis::Vector{Vector{T}}, quadratic_form::Function, C::T) where T <: Integer
    # input check
    n = length(basis)
    n == 0 && return basis
    m = length(basis[1])
    prod([length(v) == m for v in basis]) || error("lengths of generatars are different")
    m < n && error("number of vectors greater than these dimension")

    # LLL reduction
    H = integral_LLL(basis, quadratic_form)
    LLLmat = hcat([b for b in basis]...) * H
    red_basis = [LLLmat[:, i] for i in 1:n]

    q = make_quadratic_form_coeffs(red_basis, M)
    S = zeros(Rational{T}, n)
    U = zeros(Rational{T}, n)
    L = zeros(T, n)
    x = zeros(T, n)
    S[n] = C
    out_vecs = []

    i = n
    tmp = T(floor(S[i] // q[i,i] * denominator(U[i])^2))
    Z = integer_square_root(tmp) // denominator(U[i])
    L[i] = T(floor(Z - U[i]))
    x[i] = T(ceil(-Z-U[i]) - 1)

    while true
        x[i] += 1
        while x[i] > L[i]
            i += 1
            x[i] += 1
        end
        if i > 1
            S[i-1] = S[i] - q[i,i]*(x[i] + U[i])^2
            i -= 1
            U[i] = sum([q[i,j]*x[j] for j in i+1:n])

            tmp = T(floor(S[i] // q[i,i] * denominator(U[i])^2))
            Z = IntegerSquareRoot(tmp) // denominator(U[i])
            L[i] = T(floor(Z - U[i]))
            x[i] = T(ceil(-Z-U[i]) - 1)
        else
            if x != zeros(T, n)
                v = sum([x[i]*red_basis[i] for i in 1:n])
                push!(out_vecs, [v, T(C - S[1] + q[1,1]*(x[1] + U[1])^2)])
            else
                return out_vecs 
            end
        end
    end
end

# Is LLL reduced?
function LLLcheck(b::Vector{Vector{T}}, M::Matrix{T}) where T <: Integer
    # input check
    n = length(b)
    n == 0 && error("no input vecotor")
    m = length(b[1])
    prod([length(v) == m for v in b]) || error("lengths of generatars are different")
    m == size(M,1) == size(M,2) || error("sizes of vector and bilinear matrix are different")
    m < n && error("number of vectors greater than these dimension")

    q(x, y) = transpose(x)*M*y
    mu = zeros(Rational{T}, n, n)
    bs = [zeros(Rational{T}, n, n)[:,i] for i in 1:n]
    for i in 1:n
        bs[i] = copy(b[i])
        for j in 1:i-1
            mu[i,j] = q(b[i],bs[j]) // q(bs[j],bs[j])
            abs(mu[i,j]) > 1//2 && return false
            bs[i] -= mu[i,j]*bs[j]
        end
        i > 1 && q(bs[i], bs[i]) < (3//4 - mu[i,i-1]^2) * q(bs[i-1],bs[i-1]) && return false
    end
    return true
end
