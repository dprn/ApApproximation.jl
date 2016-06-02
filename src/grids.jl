##################################
# CARTESIAN/POLAR INTERPOLATIONS #
##################################

"""
Convert an image given in cartesian coordinates to polar coordinates,
assuming it to be defined on the square centered at 0 and of side [-1,1].

Arguments:
- `imgA` : 2D vector of the image (indexed by row)
- `def_rho`, `def_theta` : definition in rho and theta variables
"""
function cart2pol{T<:Real}(imgA::Array{T,2}, def_rho::Real, def_theta::Real)
	x = -1:(2./(size(imgA,1) - 1)):1.
	y = -1:(2./(size(imgA,2) - 1)):1.

	z = CoordInterpGrid((x,y), imgA, BCperiodic, InterpQuadratic)

	theta_step = 2*pi/def_theta
	eltype(imgA)[ z[rho*sin(theta), rho*cos(theta)] for rho = 0.:1./(def_rho-1):1., theta = 2*pi - theta_step:-theta_step:0.]
end

"""
Convert an image given in cartesian coordinates to polar coordinates,
assuming it to be defined on the circle centered at 0 and of radius 1.

Arguments:
- `polA` : 2D vector of the image (indexed by row)
- `def_x`, `def_y` : definition in x and y variables
"""
function pol2cart{T<:Real}(pol::Array{T,2},
                            def_x::Integer = round(Int,size(pol,1)),
                            def_y::Integer=round(Int,size(pol,1));
                            clip = false)
    coeff = clip?sqrt(2):1.
    rho_step = coeff/(size(pol,1) - 1)
    rho = 0.:rho_step:coeff
    theta_step = 2/size(pol,2)
    theta = (1/2-theta_step:-theta_step:-3/2)*pi

    z = CoordInterpGrid((rho, theta), pol, BCperiodic, InterpQuadratic)

    f(x,y) = x^2+y^2 <= 1 ? z[sqrt(x^2+y^2), atan2(y,x)] : 1.
    eltype(pol)[ f(x,y) for x = 1.:-2./(def_x-1):-1., y = -1.:2./(def_y-1):1.]
end

function pol2cart{T<:Complex}(polA::Array{T,2}, def_x::Integer = round(Int,size(pol,1)), def_y::Integer=round(Int,size(pol,1)),clip = false)
    coeff = clip?sqrt(2):1.
    rho_step = coeff/(size(pol,1) - 1)
    rho = 0.:rho_step:coeff
    theta_step = 2/size(pol,2)
    theta = (1/2-theta_step:-theta_step:-3/2)*pi

    z_real = CoordInterpGrid((rho, theta), real(polA), BCperiodic, InterpQuadratic)
    z_imag = CoordInterpGrid((rho, theta), imag(polA), BCperiodic, InterpQuadratic)

    f(z,x,y) = x^2+y^2 <= 1 ? z[sqrt(x^2+y^2), atan2(y,x)] : 1.
    eltype(polA)[complex(f(z_real, x,y),f(z_imag, x,y)) for x = 1.:-2./(def_x-1):-1., y = -1.:2./(def_y-1):1.]
end

#####################################
# INTERPOLATIONS ON BISPECTRAL SETS #
#####################################

"""
We store the interpolation vector for the image
together with its bispectral set.
"""
immutable BispInterpolation{N, T<:Number} <: AbstractArray
    f::Array{T,2}
    E::BispectralSet{N, Float64}

	function BispInterpolation(f, E)
		(size(f) != size(E)) && error("f needs to have the same dimensions as E")
		new(f,E)
	end
end

BispInterpolation{N, T<: Number}(f::Array{T,2}, E::BispectralSet{N, Float64}) = BispInterpolation{N,T}(f, E)

camembert{N, T<:Number}(::BispInterpolation{N, T}) = N

eltype{N, T<:Number}(::BispInterpolation{N,T}) = T
size(f::BispInterpolation) = size(f.f)
size(f::BispInterpolation, n) = size(f.f,n)
length(f::BispInterpolation, n) = prod(size(f))
linearindexing(::Type{BispInterpolation}) = Base.LinearFast()
getindex(f::BispInterpolation, i::Int) = f.f[i]
getindex(f::BispInterpolation, I) = [f.f[i] for i in I]
getindex(f::BispInterpolation, ::Colon) = [f.f[i] for i in 1:length(f)]
getindex(f::BispInterpolation, i::Int, j::Int) = f.f[i,j]
getindex(f::BispInterpolation, I, J) = [f.f[i,j] for i in I, j in J]

ndims(af::BispInterpolation) = length(size(af.f))

start(::BispInterpolation) = 1
next(f::BispInterpolation, state) = (f[state], state+1)
done(f::BispInterpolation, s) = s > length(f)

-{T<:Number, N}(f::BispInterpolation{N,T}, g::BispInterpolation{N,T}) =
    f.E==g.E ? BispInterpolation{N,T}(f.f-g.f, E) : error("They have to have the same bispectral set")

abs{N,T}(f::ApApproximation.BispInterpolation{N,T}) = BispInterpolation{N,Float64}(abs(f.f),f.E)

"""
Evaluates the interpolated function at the given frequency `x`.
If `x` is not in the bispectral set, it returns 0.
"""
function evaluate(f::BispInterpolation, x::Frequency)
    i,n = findin(x, f.E)
    (i != 0) ? f[i,n] : 0.
end

"""
Rotate the interpolated function (useful for tests)
"""
rotate(f::BispInterpolation, x::Real) = BispInterpolation(circshift(f.f,(0,round(Int,x*camembert(f)/(2*pi)))), f.E)

"""
Translates in the frequency space (useful for tests).
"""
function translate{T<:Complex, N}(af::BispInterpolation{N,T}, ρ::Real, θ::Real)
    af_trans = T[af[i,j]*exp(im*af.E[i].λ*ρ*cos(float(af.E[i,j].ω)-θ)) for i in 1:size(af.E,1), j in 1:size(af.E,2)]
    BispInterpolation(af_trans, af.E)
end

###########################
# INTERPOLATION FUNCTIONS #
###########################

"""
Interpolates an image, given as a 2D vector, on the given bispectral set.
"""
function cartesian2bispectral{T<:Real, R<:Real, N}(img::Array{T,2}, E::BispectralSet{N,R})
    ρs = radii(E)
    angl = camemebert_angles(E)
    def_rho = length(ρs)
    max_angles = maximum(angl)
    def_theta = N * max_angles

    imgPol = cart2pol(img, def_rho, def_theta)

    f = zeros(eltype(img), size(E))
    pos = 0
    for i in 1:length(ρs)
        for j = 1:angl[i]
            pos += 1
            for k = 1:N
                f[pos, k] = imgPol[i, j + (k-1) * max_angles]
            end
        end
    end

	BispInterpolation{N,T}(f, E)
end

"""
Interpolates the values on a bispectral set in an appropriate polar grid
"""
function bispectral2pol(f::BispInterpolation)
    ρs = [f.E[i].λ for i in 1:size(f.E,1)] |> unique
    θs = [f.E[i,j].ω for i in 1:size(f.E,1), j in 1:size(f.E,2)] |> vec |> unique
    z = [evaluate(f, Frequency(ρ,θ)) for ρ in ρs, θ in θs]
    mapslices(interpolate!, z, 2)
end

"""
Auxiliary function for `bispectral2pol`, which linearly interpolates the zero values of a vector
"""
function interpolate!(x)
    v = find(x)
    v == [] && (return x)
    if x[end]==0
        x[end] = x[v[end]]
        push!(v, length(x))
    end
    for (a,b) in zip(v, circshift(v,-1))
        for i in a+1:b-1
            dist(p,q) = q-p #min(q-p, length(x)+p-q)
            x[i] = (x[b]*dist(a,i)+x[a]*dist(i,b))/dist(a,b)
        end
    end
    x
end
