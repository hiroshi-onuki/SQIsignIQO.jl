# left ideal of the maximal order <1, i, (i + j)/2, (1 + ij)/2>
struct LeftIdeal
    b1::QOrderElem
    b2::QOrderElem
    b3::QOrderElem
    b4::QOrderElem
end

function LeftIdeal(basis::Vector{QOrderElem})
    return LeftIdeal(basis[1], basis[2], basis[3], basis[4])
end

function Base.:*(x::QOrderElem, I::LeftIdeal)
    return LeftIdeal(x*I.b1, x*I.b2, x*I.b3, x*I.b4)
end

function ideal_to_matrix(I::LeftIdeal)
    return hcat([[b[i] for i in 1:4] for b in [I.b1, I.b2, I.b3, I.b4]]...)
end

function norm(I::LeftIdeal)
    D = det(ideal_to_matrix(I))
    return Integer(sqrt(abs(D)))
end

# left O-ideal Ox + ON
function LeftIdeal(x::QOrderElem, N::Integer)
    basis = [QOrderElem(1,0,0,0,x.p,x.nj), QOrderElem(0,1,0,0,x.p,x.nj), QOrderElem(0,0,1,0,x.p,x.nj), QOrderElem(0,0,0,1,x.p,x.nj)]
    Ox = [[(b*x)[i] for i in 1:4] for b in basis]
    ON = [[N,0,0,0],[0,N,0,0],[0,0,N,0],[0,0,0,N]]
    basis = get_basis(vcat(Ox, ON))
    return LeftIdeal([QOrderElem(b[1], b[2], b[3], b[4], x.p, x.nj) for b in basis])
end

# smallest element in I
function small_element(I::LeftIdeal)
    p = I.b1.p
    q(x, y) = quadratic_form(QOrderElem(x, p), QOrderElem(y, p))

    # LLL reduction
    Imatrix = ideal_to_matrix(I)
    H = integral_LLL([Imatrix[:, i] for i in 1:4], q)
    LLLmat = Imatrix * H
    red_basis = [LLLmat[:, i] for i in 1:4]

    N = norm(I)
    C = BigInt(2)^137*N

    q = make_quadratic_form_coeffs(red_basis, q)
    S = zeros(Rational{Integer}, 4)
    U = zeros(Rational{Integer}, 4)
    L = zeros(Integer, 4)
    x = zeros(Integer, 4)
    S[4] = C

    i = 4
    tmp = div(S[i] * denominator(U[i])^2, q[i,i])
    Z = integer_square_root(tmp) // denominator(U[i])
    L[i] = Integer(floor(Z - U[i]))
    x[i] = Integer(ceil(-Z-U[i]) - 1)

    while true
        x[i] += 1
        while x[i] > L[i]
            i += 1
            x[i] += 1
        end
        if i > 1
            S[i-1] = S[i] - q[i,i]*(x[i] + U[i])^2
            i -= 1
            U[i] = sum([q[i,j]*x[j] for j in i+1:4])

            tmp = div(S[i] * denominator(U[i])^2, q[i,i])
            Z = integer_square_root(tmp) // denominator(U[i])
            L[i] = Integer(floor(Z - U[i]))
            x[i] = Integer(ceil(-Z - U[i]) - 1)
        else
            if x != zeros(Integer, 4)
                v = sum([x[i]*red_basis[i] for i in 1:4])
                alpha = QOrderElem(v[1], v[2], v[3], v[4], p)
                a, b, found = sum_of_two_squares(BigInt(2)^137 - div(norm(alpha), N))
                if found
                    return alpha, a, b, true
                end
            else
                return QQrderElem(0), 0, 0, false
            end
        end
    end
end

# return coefficients q_i,j s.t. Nrd(x) = sum_i q_i,i*(x_i + sum_j q_i,j*x_j)^2, where x = sum_i x_iI[i].
# See p.103 in H. Cohen, A Course in Computational Algebraic Number Theory.
function make_quadratic_form_coeffs(basis::Vector{Vector{T}}, quadratic_form::Function) where T <: Integer
    n = length(basis)
    C = zeros(Rational{T}, n, n)
    q = zeros(Rational{T}, n, n)

    for i in 1:n
        C[i, i] = quadratic_form(basis[i], basis[i])
        for j in i+1:n
            C[i, j] = quadratic_form(basis[i], basis[j])
        end
    end

    for i in 1:n
        q[i, i] = C[i, i] - (i > 1 ? sum([q[j, j] * q[j, i]^2 for j in 1:i-1]) : 0)
        for j in i+1:n
            q[i, j] = (2*C[i, j] - 2*sum([q[k,k]*q[k,i]*q[k,j] for k in 1:i])) / (2*q[i, i])
        end
    end
    return q
end
