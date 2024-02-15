
# floor(sqrt(n))
function integer_square_root(n::T) where T <: Integer
    n == 0 && return T(0)
    x = n + 1
    y = n
    while y < x
        x = y
        y = div(x + div(n, x), 2)
    end
    return x
end
