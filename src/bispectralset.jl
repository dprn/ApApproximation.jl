importall Base

typealias RatInt Union{Rational, Integer}

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
end

Angle{T<:Real}(val::T, N::Int) = Angle(val, 0, N)
function Angle{T<:Real}(val::T, n::Int, N::Int)
    if 0<=val<1
        Angle{N,T}(val, mod(n, N))
    elseif val>=0
        Angle{N,T}(val%1, mod((div(val,1)+n),N))
    else   # val<0
        Angle{N,T}(1+val%1, mod((N+div(val,1)+n-1),N))
    end
end

camembert{N, T<:Real}(::Angle{N, T}) = N
eltype{N, T<:Real}(::Angle{N, T}) = T
value(x::Angle) = x.val
slice(x::Angle) = x.slice
similar(x,y...) = Angle{camembert(x), eltype(x)}(y...)

rotate(x::Angle, n::Int) = similar(x, value(x), (slice(x) + n)%camembert(x))

convert{T<:AbstractFloat}(::Type{T}, x::Angle) = ( convert(T, value(x)) + slice(x) )*2pi/camembert(x)

==(a::Angle, b::Angle) = camembert(a) == camembert(b) && ((value(a) == value(b) && slice(a) == slice(b)) || value(a) == value(b) == 0)
isless(a::Angle, b::Angle) = (a.val + a.slice) < (b.val + b.slice)

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

sin(a::Angle) = sin(convert(Float64, a))
cos(a::Angle) = cos(convert(Float64, a))

###############
# FREQUENCIES #
###############

"""
2D frequency in polar coordinates
"""
immutable Frequency{N, T<:Real}
    λ::T
    ω::Angle{N}
end

Frequency{R<:Real}(a::R, x...) = Frequency{x[end], typeof(a)}(a, Angle(x...))

slice(x::Frequency) = slice(x.ω)
camembert{T<:Real,N}(::Frequency{N,T}) = N
Base.angle(x::Frequency) = convert(Float64, x.ω)
rotate(x::Frequency, n::Int) = Frequency(x.λ, rotate(x.ω,n))

==(x::Frequency, y::Frequency) = x.λ == y.λ && x.ω == y.ω
isless(x::Frequency, y::Frequency) = x.λ < y.λ || (x.λ == y.λ && x.ω < y.ω)

-(a::Frequency) = Frequency(a.λ, -a.ω)

"""
Product of two frequencies.
"""
composition{T<:Real,N}(x::Frequency{N,T}, y::Frequency{N,T}) = Frequency{N,T}(x.λ*y.λ, x.ω-y.ω)

"""
Cartesian coordinates of the frequency.
"""
cart(x::Frequency) = (x.λ*cos(x.ω), x.λ*sin(x.ω))
cart{T<:Real, N}(v::Vector{Frequency{N, T}}) = (T[x.λ*cos(angle(x)) for x in v], T[x.λ*sin(angle(x)) for x in v])

normalize(x::Frequency) = rotate(x, -slice(x))

###################
# BISPECTRAL SETS #
###################

"""
A bispectral set is a special vector of frequencies.
"""
immutable BispectralSet{N, T<:Real} <: AbstractArray{Frequency{N,T},1}
    pts::Vector{Frequency{N, T}}
end

function BispectralSet{T<:Real, N}(pts::Frequency{N, T})
    @assert all(x->x<=1,[slice(x) for x in pts]) "Frequencies must be in [0,2π/N)"
#     @assert all(x->x==N,[camembert(x) for x in pts]) "Number of slices must be the N for all frequencies"
    x = issorted(pts) ? pts : sort(pts)
    BispectralSet{camembert(pts[1]), T}(x)
end

# BispectralSet{T<:Real, N}(N::Int,pts::Vector{Frequency{T}}) = BispectralSet{T}(N,pts)

camembert{N, T<:Real}(::BispectralSet{N,T}) = N

Base.size(E::BispectralSet) = (length(E.pts), camembert(E))
Base.size(E::BispectralSet, n) = size(E)[n]
# Base.linearindexing(::Type{BispectralSet}) = Base.LinearSlow()
function Base.getindex(E::BispectralSet, i::Int)
    if 1<= i <= size(E,1)
        E.pts[i]
    else
        r = mod1(i, size(E,1))
        d = div(i-r, size(E,1))
        E[r, d+1]
    end
end
Base.getindex(E::BispectralSet, i::Int, n::Int) = 1<=n<=camembert(E) ? rotate(E.pts[i], n-1) : BoundsError(E, [i,n])
Base.getindex{T<:Real,N}(E::BispectralSet{N,T}, ::Colon, ns) = Frequency{N,T}[E[i,n] for i in 1:size(E,1), n in ns]
Base.getindex{T<:Real,N}(E::BispectralSet{N,T}, is, ns) = Frequency{N,T}[E[i,n] for i in is, n in ns]
Base.getindex{T<:Real,N}(E::BispectralSet{N,T}, ::Colon) = vec(E[1:size(E,1), 1:size(E,2)])

Base.start(::BispectralSet) = 1
Base.next(E::BispectralSet, state) = (E[state], state+1)
Base.done(E::BispectralSet, s) = s > prod(size(E))

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
        @inbounds v[pos] = Frequency(radii[i], j//angles, N) # j*2pi/N/angles)
        pos += 1
    end
    v
end

"""
Simple search by bisection
"""
function Base.findin(x::Frequency, E::BispectralSet)
    y = normalize(x)
    left = 1
    (y == E[left]) && (return left, slice(x)+1)
    right = size(E,1)
    (y == E[right]) && (return right, slice(x)+1)
    middle = round(Int, size(E,1)/2)
    for i in 1:size(E,1)
        ((middle <= left) || middle >= right) && break
        if y == E[middle]
            return middle, slice(x)+1
        elseif    y > E[middle]
            left = middle
        else                    #if y<E[middle]
            right = middle
        end
        middle = round(Int, (right+left)/2)
    end
    return 0,0
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
        while ρs[i] == E[pos_ρ].λ
            angles[i] += 1
            pos_ρ+=1
        end
    end
    angles
end

#############################
# EXTRACT RADIAL DEFINITION #
#############################

Base.den{T<:Integer}(x::Vector{Rational{T}}) = map(den, x)
Base.angle(E::BispectralSet) = Rational{Int}[E[i].ω.val for i in 1:size(E,1)] |> unique
radii{T<:Real,N}(E::BispectralSet{N,T}) = T[E[i].λ for i in 1:size(E,1)] |> unique
angles_def(E::BispectralSet) = angle(E) |> den |> lcm

function Base.gcd{T<:Integer}(a::Rational{T}, b::Rational{T})
    l = lcm(den(a),den(b))
    gcd(num(a)*gcd(den(a),l), num(b)*gcd(den(b),l))//l
end
Base.gcd{T<:Real}(a::T, b::T) = gcd(rationalize(Int,a), rationalize(Int, b))
Base.gcd{T<:Real, R<:Integer}(a::Rational{R}, b::T) = gcd(a, rationalize(Int, b))

radius_def(E::BispectralSet) =  radius_def(radii(E))
function radius_def{T<:Real}(x::Vector{T})
    v = (x-circshift(x,1))[2:end]
    foldl(gcd, v)
end
