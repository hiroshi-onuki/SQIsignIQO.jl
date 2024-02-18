
# return (a, b) s.t. M(a, b)^t = 0 mod l^e
function kernel_coefficients(M::Matrix{T}, l::Int, e::Int) where T <: Integer
    N = T(l)^2
    if M[1, 1] % l != 0
        return (M[1, 2] *  invmod(-M[1, 1], N)) % N, 1
    elseif M[2, 1] % l != 0
        return (M[2, 2] *  invmod(-M[2, 1], N)) % N, 1
    elseif M[1, 2] % l != 0
        return 1, (M[1, 1] *  invmod(-M[1, 2], N)) % N
    else
        return 1, (M[2, 1] *  invmod(-M[2, 2], N)) % N
    end
end

# return the kernel of the l^e-isogeny determined by the kernel matrix M
function kernel_gen_power_of_prime(xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, a24::Proj1{T},
    M::Matrix{T}, l::Int, e::Int) where T <: Integer
    a, b = kernel_coefficients(M, l, e)
    if a == 1
        return ladder3pt(b, xP, xQ, xPQ, a24)
    else
        return ladder3pt(a, xQ, xP, xPQ, a24)
    end
end

function ideal_to_matrix(I::LeftIdeal)
end

function short_ideal_to_isogeny(I::LeftIdeal, a24::Proj1{T}, xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T},
    M::Matrix{BigInt}, e::Int, tdata::TorsionData) where T <: RingElem
    xPd = xDBLe(xP, a24, ExponentFull - e)
    xQd = xDBLe(xQ, a24, ExponentFull - e)
    xPQd = xDBLe(xPQ, a24, ExponentFull - e)

    # 2^e-isogeny
    alpha = primitive_element(I)
    M0 = alpha[1]*[1 0; 0 1] + alpha[2]*tdata.Matrices_2e[1] + alpha[3]*tdata.Matrices_2e[2] + alpha[4]*tdata.Matrices_2e[3]
    ker = kernel_gen_power_of_prime(xPd, xQd, xPQd, a24, M*M0, 2, e)
    a24d, images = two_e_iso(a24, ker, e, [xP, xQ, xPQ])

    beta, a, b, found = two_e_good_element(I, ExponentForIsogeny)
    !found && throw(ArgumentError("No good element found"))
    beta = involution(beta) # beta = J \bar{I}
    M_beta = beta[1]*[1 0; 0 1] + beta[2]*tdata.Matrices_2e[1] + beta[3]*tdata.Matrices_2e[2] + beta[4]*tdata.Matrices_2e[3]
    M1 = [a 0; 0 a] + b * tdata.Matrices_2e[1]

    M = ideal_to_matrix(I)
    xP, xQ, xPQ = I.basis[1], I.basis[2], I.basis[3]
    a24 = I.basis[1].a24
    return kernel_gen_power_of_prime(xP, xQ, xPQ, a24, M, l, e)
end