# Sample a random ideal of prime norm D
function sample_random_ideal(D::Integer)
    @assert D % 4 == 3
    gamma = FullRepresentInteger(D * BigInt(2)^Log2p)
    a = rand(1:D-1)
    return LeftIdeal(gamma * (QOrderElem(a)  + Quoternion_i), D)
end

