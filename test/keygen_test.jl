using KaniSQIsign

function keygen_test(param::Module)
    Fp2, Fp2_i, cdata = param.make_field_curve_torsions()
    param.key_gen(cdata)
end

keygen_test(KaniSQIsign.Level1)