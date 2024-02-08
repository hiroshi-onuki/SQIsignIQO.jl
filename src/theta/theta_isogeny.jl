
function two_two_isogeny_8torsion(domain::ThetaNullLv2{T}, T1::ThetaPtLv2{T}, T2::ThetaPtLv2{T},
        image_points::Vector{ThetaPtLv2{T}}, hadamard::Bool) where T <: RingElem

    F = parent(domain.a)
    xA, xB, _, _ = Hadamard(square(T1))
    zA, tB, zC, tD = Hadamard(square(T2))

    xA_inv, zA_inv, tB_inv = batched_inversion([xA, zA, tB])

    A = F(1)
    B = xB * xA_inv
    C = zC * zA_inv
    D = tD * tB_inv * B

    _, _, _, BBinv, CCinv, DDinv = precomputation!(domain)
    B_inv = BBinv * B
    C_inv = CCinv * C
    D_inv = DDinv * D

    if hadamard
        A, B, C, D = Hadamard(A, B, C, D)
    end

    # images of points
    ret = Vector{ThetaPtLv2{T}}(undef, length(image_points))
    for i in 1:length(image_points)
        x, y, z, t = Hadamard(square(image_points[i]))
        y *= B_inv
        z *= C_inv
        t *= D_inv
        if hadamard
            x, y, z, t = Hadamard(x, y, z, t)
        end
        ret[i] = ThetaPtLv2(x, y, z, t)
    end

    return ThetaNullLv2(A, B, C, D), ret
end

function two_two_isogeny_8torsion_to_product(domain::ThetaNullLv2{T}, T1::ThetaPtLv2{T}, T2::ThetaPtLv2{T},
        image_points::Vector{ThetaPtLv2{T}}) where T <: RingElem

    F = parent(domain.a)
    xA, xB, _, _ = Hadamard(square(Hadamard(T1)))
    zA, tB, zC, tD = Hadamard(square(Hadamard(T2)))

    xA_inv, zA_inv, tB_inv, xB_inv, zC_inv, tD_inv = batched_inversion(
            [xA, zA, tB, xB, zC, tD]
        )
    A = F(1)
    B = xB * xA_inv
    C = zC * zA_inv
    D = tD * tB_inv * B
    B_inv = xB_inv * xA
    C_inv = zC_inv * zA
    D_inv = tD_inv * tB * B_inv

    # images of points
    ret = Vector{ThetaPtLv2{T}}(undef, length(image_points))
    for i in 1:length(image_points)
        x, y, z, t = Hadamard(square(Hadamard(image_points[i])))
        y *= B_inv
        z *= C_inv
        t *= D_inv
        ret[i] = ThetaPtLv2(x, y, z, t)
    end

    return ThetaNullLv2(A, B, C, D), ret
end

# (2^n, 2^n)-isogeny with kernel <4*P1P2, 4*P2Q2>. P1Q1P2Q2 = (x(P1 + Q1), x(P2 + Q2))
function product_isogeny_no_strategy(a24_1::Proj1{T}, a24_2::Proj1{T},
    P1P2::CouplePoint{T}, Q1Q2::CouplePoint{T}, P1Q1P2Q2::CouplePoint{T},
    image_points::Vector{CouplePoint{T}}, n::Integer) where T <: RingElem

    push!(image_points, P1P2)
    push!(image_points, Q1Q2)

    P1P2_8 = double_iter(P1P2, a24_1, a24_2, n-1)
    Q1Q2_8 = double_iter(Q1Q2, a24_1, a24_2, n-1)

    domain, image_points = gluing_isogeny(a24_1, a24_2, P1P2_8, Q1Q2_8, P1Q1P2Q2, image_points, n)

    for k in 1:n-1
        Tp1 = double_iter(domain, image_points[end - 1], n - k - 1)
        Tp2 = double_iter(domain, image_points[end], n - k - 1)

        if k == n - 2
            domain, image_points = two_two_isogeny_8torsion(domain, Tp1, Tp2, image_points, false)
        elseif k == n - 1
            # remove kernel generators
            pop!(image_points)
            pop!(image_points)
            domain, image_points = two_two_isogeny_8torsion_to_product(domain, Tp1, Tp2, image_points)
        else
            domain, image_points = two_two_isogeny_8torsion(domain, Tp1, Tp2, image_points, true)
        end
    end

    domain, image_points = splitting_isomorphism(domain, image_points)
    return split_to_product(domain, image_points)
end

function product_isogeny(a24_1::Proj1{T}, a24_2::Proj1{T},
    P1P2::CouplePoint{T}, Q1Q2::CouplePoint{T}, P1Q1P2Q2::CouplePoint{T},
    image_points::Vector{CouplePoint{T}}, n::Integer, strategy::Vector{Int}) where T <: RingElem

    push!(image_points, P1P2)
    push!(image_points, Q1Q2)

    strategy_idx = 1
    level = Int[0]
    prev = 0

    ker1 = P1P2
    ker2 = Q1Q2
    while prev != n - 1
        push!(level, strategy[strategy_idx])

        ker1 = double_iter(ker1, a24_1, a24_2, strategy[strategy_idx])
        ker2 = double_iter(ker2, a24_1, a24_2, strategy[strategy_idx])
        push!(image_points, ker1)
        push!(image_points, ker2)

        prev += strategy[strategy_idx]
        strategy_idx += 1
    end
    pop!(image_points)
    pop!(image_points)
    pop!(level)

    domain, image_points = gluing_isogeny(a24_1, a24_2, ker1, ker2, P1Q1P2Q2, image_points, n)

    for k in 1:n-1
        prev = sum(level)
        ker1 = image_points[end - 1]
        ker2 = image_points[end]

        while prev != (n - 1 - k)
            push!(level, strategy[strategy_idx])

            ker1 = double_iter(domain, ker1, strategy[strategy_idx])
            ker2 = double_iter(domain, ker2, strategy[strategy_idx])
            push!(image_points, ker1)
            push!(image_points, ker2)

            prev += strategy[strategy_idx]
            strategy_idx += 1
        end
        pop!(image_points)
        pop!(image_points)
        pop!(level)

        if k == n - 2
            domain, image_points = two_two_isogeny_8torsion(domain, ker1, ker2, image_points, false)
        elseif k == n - 1
            # remove kernel generators
            pop!(image_points)
            pop!(image_points)
            domain, image_points = two_two_isogeny_8torsion_to_product(domain, ker1, ker2, image_points)
        else
            domain, image_points = two_two_isogeny_8torsion(domain, ker1, ker2, image_points, true)
        end
    end

    domain, image_points = splitting_isomorphism(domain, image_points)
    return split_to_product(domain, image_points)
end