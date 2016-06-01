module ApApproximation

export Frequency,
BispectralSet,
cart,
camembert,
evaluate,
ap,
iap,
BesselMatrix,
cartesian2bispectral,
bispectral2pol,
pol2cart,
BispInterpolation,
rotate,
translate

const TOL=1e-5

include("bispectralset.jl")

using Grid
include("grids.jl")

include("bessel.jl")
include("ap_interpolation.jl")

end # module
