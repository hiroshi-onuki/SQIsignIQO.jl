
function FullRepresentInteger(M::Integer, p::BigInt)
    counter = 0
    found = false
    while !found && counter < KLPT_repres_num_gamma_trial
        m = Integer(sqrt(div(4*M, p)))
        z = rand(-m:m)
        md = Integer(sqrt(div(4*M - z^2, p)))
        w = rand(-md:md)
        Md = 4*M - p*(z^2 - w^2)
        x, y, found = Cornacchia(Md)
        if found || (x - w) % 2 == 0 || (y - z) % 2 == 0
            found = false
            counter += 1
        end
    end
    if found
        return QOrderElem(div(x - w, 2), div(y - z, 2), z, w, p), found
    else
        return QOrderElem(0, p), found
    end
end