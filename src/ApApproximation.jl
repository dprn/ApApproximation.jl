module ApApproximation

export Frequency,
BispectralSet,
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

include("bispectralset.jl")

using Grid
include("grids.jl")

include("bessel.jl")
include("ap_interpolation.jl")

end # module
