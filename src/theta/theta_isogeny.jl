
function two_two_isogeny_8torsion(domain::ThetaNullLv2{T}, T1::ThetaPtLv2{T}, T2::ThetaPtLv2{T},
        is_domain_dual::Bool=false, is_codomain_dual::Bool=false) where T <: RingElem
    F = parent(domain.a)
    if is_domain_dual
        xA, xB, _, _ = Hadamard(square(hadamard(T1)))
        zA, tB, zC, tD = Hadamard(square(hadamard(T2)))
    else
        xA, xB, _, _ = Hadamard(square(T1))
        zA, tB, zC, tD = Hadamard(square(T2))
    end

    if !is_domain_dual && domain.precomputed
        xA_inv, zA_inv, tB_inv = batched_inversion([xA, zA, tB])

        A = F(1)
        B = xB * xA_inv
        C = zC * zA_inv
        D = tD * tB_inv * B

        _, _, _, BBinv, CCinv, DDinv = domain.precomputation
        B_inv = BBinv * B
        C_inv = CCinv * C
        D_inv = DDinv * D
    else
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
    end

    if is_codomain_dual
        return ThetaNullLv2(A, B, C, D)
    else
        return ThetaNullLv2(Hadamard([A, B, C, D]))
    end
end