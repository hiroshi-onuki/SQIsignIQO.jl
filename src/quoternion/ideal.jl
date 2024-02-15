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
    red_basis = [QOrderElem(LLLmat[:, i], p) for i in 1:4]

    B = 6
    N = norm(I)
    while true
        c1, c2, c3, c4 = [rand(-B:B) for i in 1:4]
        x = c1*red_basis[1] + c2*red_basis[2] + c3*red_basis[3] + c4*red_basis[4]
        if x != 0
            n = div(norm(x), N)
            _, _, found = sum_of_two_squares(BigInt(2)^137 - n)
            if found
                return x
            end
        end
    end
end