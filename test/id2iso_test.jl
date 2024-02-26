using Nemo
using KaniSQIsign

#a = QOrderElem(168082123654050744778996387641574941677314667606066074978750034663303195240598928771394263744915916019712993146821097074055019317214177783512247070593762705522957193601412, 290475121935060369528280295437947531300840252723913376278647347304065523783601994054482835948285018743520397621340134881884658475384029520874884280915263793538215168188281, 972444482673414445000286753217595837683916086717251938159230044415087772058527255439372289073551914230201772368055038666995203906987, -2016948722110680185238295871178154499691515430257255802134678920457372233599340826676216189261523778821756564157794034320470816525733)

function test_id2iso(param::Module, n::Int)
    p = param.p
    e1 = param.ExponentForIsogeny
    _, _, cdata = param.make_field_curve_torsions()

    N = 337416411430778090000393699717 # large prime
    e = n*Int(ceil(log(2, p)))

    a, found = param.FullRepresentInteger(N*BigInt(2)^e * param.ExtraDegree)
    if !found
        println("Not found")
        return
    end
    a = div(a, gcd(a))
    println(factor(ZZ(param.norm(a))))
    I = param.LeftIdeal(a, BigInt(2)^(n*e1) * param.ExtraDegree)
    a24 = param.ideal_to_isogeny_from_O0(I, n*e1, cdata)
    println(a24)
end

test_id2iso(KaniSQIsign.Toy17, 4)