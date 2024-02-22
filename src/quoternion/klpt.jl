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

function RandomEquivalentPrimeIdeal(I::LeftIdeal)
    counter = 0
    found = false
    J = LeftIdeal([QOrderElem(0), QOrderElem(0), QOrderElem(0), QOrderElem(0)])
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
        if is_prime(div(norm(beta), N))
            found = true
            J = ideal_transform(I, beta, N)
        end
    end
    return J, found
end