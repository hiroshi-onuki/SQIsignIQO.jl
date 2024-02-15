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

function norm(I::LeftIdeal)
    return gcd(norm(I.b1), norm(I.b2), norm(I.b3), norm(I.b4))
end

# left O-ideal Ox + ON
function LeftIdeal(x::QOrderElem, N::Integer)
    basis = [QOrderElem(1,0,0,0,x.p,x.nj), QOrderElem(0,1,0,0,x.p,x.nj), QOrderElem(0,0,1,0,x.p,x.nj), QOrderElem(0,0,0,1,x.p,x.nj)]
    Ox = [[(b*x)[i] for i in 1:4] for b in basis]
    ON = [[N,0,0,0],[0,N,0,0],[0,0,N,0],[0,0,0,N]]
    basis = get_basis(vcat(Ox, ON))
    return LeftIdeal([QOrderElem(b[1], b[2], b[3], b[4], x.p, x.nj) for b in basis])
end
