importall Base

##########
# ANGLES #
##########

"""
Encodes the angle
    2pi/N*(val + slice)
where val is always assumed to be in [0,1)
"""
immutable Angle{N, T<:Real}
    val::T          # value in the camembert slice
    slice::Int                  # camembert slice

    function Angle(val, n)
        if 0<=val<1
            new(val, mod(n, N))
        elseif val>=0
            new(val%1, mod((div(val,1)+n),N))
        else   # val<0
            new(1+val%1, mod((N+div(val,1)+n-1),N))
        end
    end
end

Angle{T<:Real}(val::T, N::Int) = Angle(val, 0, N)
Angle{T<:Real}(val::T, n::Int, N::Int) = Angle{N,T}(val, n)

camembert{N, T<:Real}(::Angle{N, T}) = N
eltype{N, T<:Real}(::Angle{N, T}) = T
value(x::Angle) = x.val
slice(x::Angle) = x.slice

rotate(x::Angle, n::Int) = Angle(value(x), (slice(x) + n)%camembert(x), camembert(x))

zero{N,T<:Real}(::Type{Angle{N,T}}) = Angle{N,T}(0., 0)
convert{T<:AbstractFloat}(::Type{T}, x::Angle) = ( convert(T, value(x)) + slice(x) )*2pi/camembert(x)

==(a::Angle, b::Angle) = camembert(a) == camembert(b) && ((value(a) == value(b) && slice(a) == slice(b)) || value(a) == value(b) == 0)
isless(a::Angle, b::Angle) = (value(a) + slice(a)) < (value(b) + slice(b))

function +(a::Angle, b::Angle)
    if camembert(a) == camembert(b)
        Angle(value(a)+value(b), slice(a) + slice(b), camembert(a))
    else
        error("Angles needs to have the same camembert slices")
    end
end

+{T<:Real}(a::Angle, b::T) = Angle(a.val+b, a.slice, camembert(a))
+{T<:Real}(b::T, a::Angle) = a + b

-(a::Angle) = Angle(-a.val, -a.slice, camembert(a))
-(a::Angle, b::Angle) = a+(-b)
-{T<:Real}(a::Angle, b::T) = -a+b
-{T<:Real}(b::T, a::Angle) = b+(-a)

*{T<:Real}(a::Angle, b::T) = Angle(a.val*b, a.slice*b, camembert(a))
*{T<:Real}(b::T, a::Angle) = a*b

for f in [:cos, :sin, :abs]
    @eval ($f)(a::Angle) = ($f)(convert(Float64, a))
end

opposite(x::Angle) = Angle(value(x)+camembert(x)/2, slice(x), camembert(x))

###############
# FREQUENCIES #
###############

"""
2D frequency in polar coordinates
"""
immutable Frequency{N, T<:Real}
    λ::T
    ω::Angle{N, T}

    Frequency(λ, ω) = λ == 0 ? new(0., 0*ω) : new(λ, ω)
end

Frequency{R<:Real}(a::R, b::Angle) = Frequency{camembert(b), typeof(a)}(a, b)
Frequency{R<:Real}(a::R, x...) = Frequency{x[end], typeof(a)}(a, Angle(x...))

λ(x::Frequency) = x.λ
ω(x::Frequency) = x.ω
camembert{N,T<:Real}(::Frequency{N,T}) = N
slice(x::Frequency) = slice( ω(x) )
eltype{N, T<:Real}(::Frequency{N,T}) = T
show(io::IO, x::Frequency) = print(io, "Frequency{$(camembert(x)), $(eltype(x))}($(λ(x)), $(slice(x.ω) + value(x.ω)))")

Base.angle(x::Frequency) = convert(Float64, x.ω)
rotate(x::Frequency, n::Int) = Frequency(λ(x), rotate(ω(x),n))

==(x::Frequency, y::Frequency) = λ(x) == λ(y) && ω(x) == ω(y)
isless(x::Frequency, y::Frequency) = λ(x) < λ(y) || (λ(x) == λ(y) && ω(x) < ω(y))

-(a::Frequency) = Frequency(λ(a), opposite(ω(a)))

*{T<:Real}(x::Frequency, b::T) = Frequency(b*λ(x), ω(x))
*{T<:Real}(b::T, x::Frequency) = x*b

function +(x::Frequency, y::Frequency)
    @assert(camembert(x) == camembert(y))
    Θ = (value(ω(x)) + slice(x)) + atan2( λ(y)sin(ω(y)-ω(x)), λ(x) + λ(y)cos(ω(y)-ω(x)) )*camembert(x)/(2pi)
    ρ = sqrt(abs( λ(x)^2 + λ(y)^2 +2(λ(x)*λ(y))*cos(ω(x)-ω(y))))
    Frequency(ρ, Θ, camembert(x))
end
-(a::Frequency, b::Frequency) = a + (-b)

"""
Product of two frequencies.
"""
composition(x::Frequency, y::Frequency) = Frequency(λ(x)*λ(y), ω(x)-ω(y))

"""
Cartesian coordinates of the frequency.
"""
cart(x::Frequency) = (λ(x)cos(ω(x)), λ(x)sin(ω(x)))
cart{T<:Real, N}(v::Vector{Frequency{N, T}}) = (T[λ(x)*cos(ω(x)) for x in v], T[λ(x)*sin(ω(x)) for x in v])

normalize(x::Frequency) = rotate(x, -slice(x))

"""
Approximate equality between frequencies
"""
approx_eq(a::Frequency) = approx_eq(a, zero(typeof(a)))
approx_eq(a::Frequency, b::Frequency) = abs(λ(a)-λ(b)) <= TOL && dist(ω(a),ω(b)) <= TOL/λ(a)

"""
Unique w.r.t. approximate equality
"""
approx_unique{N, T<:Real}(v::Vector{Frequency{N,T}}) = approx_unique(approx_eq, v)

###################
# BISPECTRAL SETS #
###################

"""
A bispectral set is a special vector of frequencies.
"""
immutable BispectralSet{N, T<:Real} <: AbstractArray{Frequency{N,T},1}
    pts::Vector{Frequency{N, T}}

    function BispectralSet(pts::Vector{Frequency{N,T}})
        @assert all(x->x<=1,[slice(x) for x in pts]) "Frequencies must be in [0,2π/N)"
        x = approx_unique(pts)
        new(issorted(x) ? x : sort(x))
    end
end

BispectralSet{N, T<: Real}(pts::Vector{Frequency{N,T}}) = BispectralSet{N,T}(pts)

camembert{N, T<:Real}(::BispectralSet{N,T}) = N
zero{N,T<:Real}(::Type{Frequency{N,T}}) = Frequency{N,T}(0., zero(Angle{N,T}) )

eltype(E::BispectralSet) = eltype(E.pts)
size(E::BispectralSet) = (length(E.pts), camembert(E))
size(E::BispectralSet, n) = size(E)[n]
function getindex(E::BispectralSet, i::Int)
    if 1<= i <= size(E,1)
        E.pts[i]
    else
        r = mod1(i, size(E,1))
        d = div(i-r, size(E,1))
        E[r, d+1]
    end
end
getindex(E::BispectralSet, i::Int, n::Int) = 1<=i<=size(E,1) && 1<=n<=size(E,2) ? rotate(E[i], n-1) : BoundsError(E, [i,n])
getindex(E::BispectralSet, ::Colon, ns) = eltype(E)[E[i,n] for i in 1:size(E,1), n in ns]
getindex(E::BispectralSet, is, ns) = eltype(E)[E[i,n] for i in is, n in ns]
getindex(E::BispectralSet, ::Colon) = vec(E[1:size(E,1), 1:size(E,2)])

start(::BispectralSet) = 1
next(E::BispectralSet, state) = (E[state], state+1)
done(E::BispectralSet, s) = s > prod(size(E))

cart(E::BispectralSet) = cart(E[:])

"""
Generates a BispectralSet given a vector (`cutoff`)
of tuples of the type (n, range), where range are the radii of the frequencies
and n the number of equispaced angles in the camembert slice for each one of them.
"""
function BispectralSet{T<:Range}(N::Int, cutoff::Vector{Tuple{Int64,T}})
    fill_freqs(x) = fill_freq_vector(N, x...)
    v = mapreduce(fill_freqs, append!, cutoff)
    BispectralSet(v)
end

"""
Auxiliary function to fill the frequency vector.
"""
function fill_freq_vector(N::Int, angles::Int, radii)
    n = length(radii)
    v = Vector{Frequency{N, Float64}}(n*angles)
    pos = 1
    for i in 1:n, j = 0:angles-1
        @inbounds v[pos] = Frequency(radii[i], j/angles, N) # j*2pi/N/angles)
        pos += 1
    end
    v
end

"""
Simple search by bisection
"""
function findin(x::Frequency, E::BispectralSet)
    y = normalize(x)
    left = 1
    approx_eq(y, E[left]) && (return left, slice(x)+1)
    right = size(E,1)
    approx_eq(y, E[right]) && (return right, slice(x)+1)
    middle = round(Int, size(E,1)/2)
    for i in 1:size(E,1)
        ((middle <= left) || middle >= right) && break
        if approx_eq(y, E[middle])
            return middle, slice(x)+1
        elseif  y > E[middle]
            left = middle
        else                    #if y<E[middle]
            right = middle
        end
        middle = round(Int, (right+left)/2)
    end
    return 0,0
end

"""
Unique w.r.t. approximate equality
"""
approx_unique{T<:Number}(v::Vector{T}) = approx_unique((x,y)->abs(x-y) <= TOL, v)
function approx_unique(eq_test::Function, v::Vector)
    x = issorted(v) ? v : sort(v)
    out = Vector{eltype(v)}()
    @inbounds for i in 1:length(v)-1
        if !eq_test(x[i], x[i+1])
            push!(out, x[i])
        end
    end
    push!(out, x[end])
    out
end



radii{N, T<:Real}(E::BispectralSet{N,T}) = T[ λ(E[i]) for i in 1:size(E,1) ] |> approx_unique
angles{N,T<:Real}(E::BispectralSet{N,T}) = T[ value( ω(E[i]) ) for i in 1:size(E,1) ] |> approx_unique

"""
Angular definition
"""
angles_def(E::BispectralSet) = angles(E) |> extract_def

"""
Radial definiton
"""
radius_def(E::BispectralSet) = radii(E) |> extract_def

"""
Extracts the minimal definition required to encode the values of the given
vecor in a regular grid (up to TOL).
"""
function extract_def(v::Vector)
    vv = (circshift(v,-1)-v)[1:end-1] |> approx_unique
    xx = round(Int, map( x->round(x, TOLexp), vv) / TOL)
    TOL*foldl(gcd, xx)
end

"""
Extracts a vector `angles` such that `angles[i]` contains the
number of elements in the camembert slice at the radius `i`.
"""
function camemebert_angles{N,T<:Real}(E::BispectralSet{N,T})
    ρs = radii(E)
    cam_angles = angles_def(E)
    angles = zeros(Int, length(ρs))
    pos_ρ = 1
    for i in 1:length(ρs)
        while ρs[i] == λ(E[pos_ρ])
            angles[i] += 1
            pos_ρ += 1
        end
    end
    angles
end

# ################### #
# PLOT BISPECTRAL SET #
# ################### #

import Plots.plot

function plot(E::BispectralSet)
    x, y = cart(E)
    ρ_max = maximum([x.λ for x in E])

    # Plot the set E
    p = scatter(x, y, marker = (:red, 2.5),
            legend = false,
            size = (315,300),
            title = "Set E")

    # Plot the lines dividing the slices
    for ang in 0:2pi/camembert(E):2pi*(1-1/camembert(E))
        plot!(t-> cos(ang)*t, t-> sin(ang)*t, 0:.1:ρ_max*(1+.1), line = (:black, :dot, 2))
    end
    p
end
