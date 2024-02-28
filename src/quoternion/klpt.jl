# Algorithm 10 in SQIsign documentation
# return a quoternion in O0 of norm M
function FullRepresentInteger(M::Integer)
    counter = 0
    found = false
    x, y, z, w = 0, 0, 0, 0
    while !found && counter < KLPT_repres_num_gamma_trial
        m = integer_square_root(div(4*M, p))
        z = rand(-m:m)
        md = integer_square_root(div(4*M - z^2, p))
        w = rand(-md:md)
        Md = 4*M - p*(z^2 + w^2)
        x, y, found = sum_of_two_squares(Md)
        if !found || (x - w) % 2 != 0 || (y - z) % 2 != 0
            found = false
            counter += 1
        end
    end
    if found
        return QOrderElem(div(x - w, 2), div(y - z, 2), z, w), found
    else
        return QOrderElem(0), found
    end
end

# Algorithm 8 in SQIsign documentation
# return J ~ I s.t. norm(J) ~ sqrt(p) and prime
function RandomEquivalentPrimeIdeal(I::LeftIdeal)
    counter = 0
    found = false
    J = LeftIdeal([QOrderElem(0), QOrderElem(0), QOrderElem(0), QOrderElem(0)])
    nJ = 0
    N = norm(I)

    # LLL reduction
    Imatrix = ideal_to_matrix(I)
    q(x, y) = quadratic_form(QOrderElem(x), QOrderElem(y))
    H = integral_LLL([Imatrix[:, i] for i in 1:4], q)
    LLLmat = Imatrix * H
    red_basis = [LLLmat[:, i] for i in 1:4]

    while !found && counter < KLPT_equiv_num_iter
        counter += 1
        c1, c2, c3, c4 = [rand(-KLPT_equiv_bound_coeff:KLPT_equiv_bound_coeff) for _ in 1:4]
        beta = c1 * red_basis[1] + c2 * red_basis[2] + c3 * red_basis[3] + c4 * red_basis[4]
        beta = QOrderElem(beta)
        nJ = div(norm(beta), N)
        if is_prime(nJ)
            found = true
            J = ideal_transform(I, beta, N)
        end
    end
    return J, nJ, found
end

# Algorithm 11 in SQIsign documentation
# return C, D s.t. gamma*j*(C + D*i)*delta in Z + I, where N = n(I), if divisible then N | norm(gamma).
function EichlerModConstraint(I::LeftIdeal, N::Integer, gamma::QOrderElem, delta::QOrderElem, divisible::Bool)
    M = HNFmod(ideal_to_matrix(I), N)
    v1 = to_vector(gamma * Quoternion_j)
    v2 = to_vector(-gamma * Quoternion_ij * delta)

    if divisible
        # gamma*j*(C + D*i)*delta in I
        M = hcat(v1, v2, M[:, 1:2])
    else
        # gamma*j*(C + D*i)*delta in Z + I
        M = hcat(v1, v2, [1, 0, 0, 0], M[:, 1:2])
    end
    M = Gauss_elimination_mod(M, N)
    if M[1, 1] == 0
        C, D = 1, 0
    elseif M[2, 2] == 0
        C, D = 0, 1
    elseif M[3, 3] == 0
        i = findfirst(x!=0 for x in M[1, 3:end]) + 2
        if i != nothing
            C, D = M[1, i], M[2, i]
        else
            C, D = 0, 1
        end
    elseif M[4, 4] == 0
        i = findfirst(x!=0 for x in M[1, 4:end]) + 2
        if i != nothing
            C, D = M[1, i], M[2, i]
        else
            C, D = 0, 1
        end
    else
        C, D = M[1, 5], M[2, 5]
    end
        
    return C, D
end