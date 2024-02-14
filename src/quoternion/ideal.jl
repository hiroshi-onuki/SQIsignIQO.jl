
struct LeftIdeal
    b1::OrderElem
    b2::OrderElem
    b3::OrderElem
    b4::OrderElem
end

function Base.(*):(x::OrderElem, I::LeftIdeal)
    return LeftIdeal(x*I.b1, x*I.b2, x*I.b3, x*I.b4)
end

function norm(I::LeftIdeal)
    return gcd(I.b1.v, I.b2.v, I.b3.v, I.b4.v)
end

