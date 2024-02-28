include("level1/prime.jl")
include("level1/klpt_constants.jl")

include("../quoternion/order.jl")
include("../quoternion/cornacchia.jl")
include("../quoternion/ideal.jl")
include("../quoternion/klpt.jl")

include("../ideal_to_isogeny/ideal_to_isogeny.jl")

const StrategyDim2 = compute_strategy(ExponentForTorsion - 4, 2, 1)

# Fp2 and values in Fp2
function make_field_curve_torsions()
    _, T = polynomial_ring(GF(p), "T")
    Fp2, Fp2_i = finite_field(T^2 + 1, "i")
    
    A0 = Fp2(0)

    # constatns from precompute/level1torsion.sage
    P2e = Point(13168888026701864374327984987573149335323956272359708223115375240549648163281*Fp2_i + 2161215833780113725019607471924727440013124162519025392101162565499276343954, 11573631817917175655308717998954552340066289950952953708354522787290655314842*Fp2_i + 9830874140031968764629465062081177133088386698799679322612340559956021219521)
    Q2e = Point(1167679482150863195523978450170412562715924357238602796020786585710888718237*Fp2_i + 195807475550494427451649608271104515847911618849141017426569917821137983194, 4862877640333236516853362674124171341900551523228923073125804599913941126908*Fp2_i + 12087566741318956310016631707446387960429051679627903688618432549189813607113)
    M_i_2e = [187804641856894237833666302824508355834507214682716838694549490552883968291 156173985651354320464277040484803017418920522788064441193166373616243357181; 168588957679671372271042170869802553932047472501894436463898332476278883566 38351782434738956352995777270585214191410724117362387945016103212571363037]
    M_ij_2e = [89525985884981099413626375992466343353718007666384484912145770930619405320 106676124532956136778423656701829964547672032461408105458899705516857934656; 86203614543351285271659433268311087941489690539722696159230979111204858683 136630438406652094773035704102627226672199931133694741727419822834835926008]
    M_1k_2e = [47293392636392707013109095476030301705868389346246377929798957148812719464 144102115607646079774418655445238838719570865336332827004198219044290149848; 197979047349103652112761733725205039071613847918896028361337091949765817183 178863031655240487173552984619063268320049549453832848709766636616642611865]
    xP79 = Proj1(17808121587231935828586247647640248406498136232625158867643836578213625687058*Fp2_i + 16286020220186406872658333272204684378108399162171061810188411060180396717158)
    xQ79 = Proj1(9666366218806412695809897752201928575473258451879582949833353562903542320936*Fp2_i + 7006434669972382365232974777910912269605319342842388545794231682129891672556)
    xPQ79 = Proj1(6177489307606325993454671866278064966740225045468462451886239598612102775407*Fp2_i + 13671566132822630425158412348447527223411962625561735380596203407515593567278)
    M_i_79 = [23 17; 6 56]
    M_ij_79 = [5 45; 17 74]
    M_1k_79 = [10 2; 34 70]
    xP3 = Proj1(7592000161467079587213592319318508362654288632426985883306372487657684985337*Fp2_i + 5871283952432172208638013510020471438703819898938916211796738350958029596793)
    xQ3 = Proj1(10274357357571942753532712008193883669393228532779273021219309419813286189574*Fp2_i + 11995073566606850132108290817491920593343697266267342692728943556512941578118)
    xPQ3 = Proj1(7592000161467079587213592319318508362654288632426985883306372487657684985337*Fp2_i + 11995073566606850132108290817491920593343697266267342692728943556512941578118)
    M_i_3 = [0 1; 2 0]
    M_ij_3 = [2 1; 0 1]
    M_1k_3 = [1 1; 1 0]
    xP5 = Proj1(16847507437177049189361159365262205020927859760801772793381914323733962762775*Fp2_i + 13605772272000940996312273595013481905001275501775595968625244296334384227937)
    xQ5 = Proj1(14903774597198409284239785752835947767246841440555037991163254909693742042546*Fp2_i + 1512147991088225013929455726144523372607638967117200841246904558202606601877)
    xPQ5 = Proj1(8423063828240058255372866588572582315449878488821245985682748700472502099684*Fp2_i + 3240355384479175514433459582131682725898481475866896002740253104610619971431)
    M_i_5 = [1 4; 2 4]
    M_ij_5 = [4 2; 3 1]
    M_1k_5 = [2 1; 0 4]

    a24_0 = A_to_a24(A0)
    xP2e = Proj1(P2e.X, P2e.Z)
    xQ2e = Proj1(Q2e.X, Q2e.Z)
    PQ2e = add(P2e, -Q2e, Proj1(A0))
    xPQ2e = Proj1(PQ2e.X, PQ2e.Z)
    xP2e_short = xDBLe(xP2e, a24_0, ExponentForIsogeny)
    xQ2e_short = xDBLe(xQ2e, a24_0, ExponentForIsogeny)
    xPQ2e_short = xDBLe(xPQ2e, a24_0, ExponentForIsogeny)
    wp_P2e_Q2e = Weil_pairing_2power(A0, P2e, Q2e, ExponentFull)

    DegreesOddTorsionBases = [3, 5, 79]
    ExponentsOddTorsionBases = [1, 1, 1]
    OddTorsionBases = [[xP3, xQ3, xPQ3], [xP5, xQ5, xPQ5], [xP79, xQ79, xPQ79]]

    Matrices_2e = [M_i_2e, M_ij_2e, M_1k_2e]
    Matrices_odd = [[M_i_3, M_ij_3, M_1k_3], [M_i_5, M_ij_5, M_1k_5], [M_i_79, M_ij_79, M_1k_79]]

    # make constants for isomorphism to the curve E_A0
    _, T = polynomial_ring(Fp2, "T")
    As = roots((256 * (T^2 - 3)^3 - 1728 * (T^2 - 4))/T^2)
    A0d = As[1]
    beta = -A0d/3
    gamma = square_root(1 / (1 - 3*beta^2))
    gamma = gamma[1]/gamma[2]
    A0dd = As[2]
    beta_d = -A0dd/3
    gamma_d = square_root(1 / (1 - 3*beta_d^2))
    gamma_d = gamma_d[1]/gamma_d[2]
    function isomorphism_to_A0(A::Proj1{FqFieldElem}, P::Proj1{FqFieldElem})
        if A == Proj1(A0)
            return P
        elseif A == Proj1(A0d)
            return Proj1(gamma*(P.X - beta*P.Z), P.Z)
        elseif A == Proj1(A0dd)
            return Proj1(gamma_d*(P.X - beta_d*P.Z), P.Z)
        else
            throw(ArgumentError("A is not A0d or A0dd"))
        end
    end

    return Fp2, Fp2_i, CurveData(A0, A0d, A0dd, a24_0, jInvariant_A(A0), P2e, Q2e, xP2e, xQ2e, xPQ2e, xP2e_short, xQ2e_short, xPQ2e_short, wp_P2e_Q2e, DegreesOddTorsionBases, ExponentsOddTorsionBases, OddTorsionBases, Matrices_2e, Matrices_odd, isomorphism_to_A0)
end