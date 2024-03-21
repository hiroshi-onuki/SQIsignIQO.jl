module KaniSQIsign

using Nemo

include("utilities/integer.jl")
include("utilities/finite_field.jl")
include("utilities/batch_inv.jl")
include("utilities/lattice.jl")
include("utilities/strategy.jl")
include("utilities/matrix.jl")

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

module Level1
using Nemo, KaniSQIsign
include("parameters/level1.jl")
end # module Level1

end # module KaniSQIsign
