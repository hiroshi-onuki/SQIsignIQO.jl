module KaniSQIsign

using Nemo
include("utilities/finite_field.jl")
include("utilities/batch_inv.jl")
include("utilities/cornacchia.jl")

include("elliptic_curves/proj1.jl")
include("elliptic_curves/full_point.jl")
include("elliptic_curves/pairing.jl")
include("elliptic_curves/montgomery.jl")
include("elliptic_curves/couple_point.jl")

include("theta/theta_structure_dim1.jl")
include("theta/theta_structure_dim2.jl")
include("theta/gluing.jl")
include("theta/theta_arithmetic.jl")
include("theta/splitting.jl")
include("theta/theta_isogeny.jl")

include("constants/constants.jl")

include("quoternion/order.jl")
include("quoternion/klpt.jl")

end # module KaniSQIsign
