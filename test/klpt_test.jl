import KaniSQIsign: FullRepresentInteger, norm

a, found = FullRepresentInteger(1000, 103)
if found
    println(a)
    println(norm(a))
end