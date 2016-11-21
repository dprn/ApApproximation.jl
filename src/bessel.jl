"""
Returns a Complex128 vector containing the discrete Bessel functions of parameter
α evaluated at λ for all n's.
"""
function discrete_bessel(N::Integer, α::Real, λ::Real)
    @assert(abs(α) <= 2π/N)
    x = Complex128[ exp(im * (λ * cos(α + 2π*k/N))) for k in 0:N-1 ]
    return circshift(fft(x), -1)
end

function discrete_bessel(x::Frequency)
    @assert in(slice(x), [0,camembert(x)-1])
    x = Complex128[ exp(im * (λ(x) * cos(rotate(ω(x), k)))) for k in 0:camembert(x)-1 ]
    return circshift(fft(x), -1)
end

"""
Returns the matrices J^n associated with the BispectralSets.

- E : Spatial BispectralSet
- F : Frequencies BispectralSet

ATTENTION: J[n] corresponds to J_{n-1}!
"""
discrete_bessel_matrix(E::BispectralSet) = discrete_bessel_matrix(E, E)
function discrete_bessel_matrix(E::BispectralSet, F::BispectralSet; verbose = true)
    @assert eltype(E) == eltype(F)
    verbose && (length(E) != length(F)) ? println("Non-square Jn matrices") : nothing

    J = Array(Complex128, camembert(E), size(E,1), size(F,1))

    # we store the values in a dict to avoid computing them more than once
    coeff = Dict{eltype(E),Vector{Complex128}}()
    for i = 1:size(E,1), j = 1:size(F,1)
        x = composition(E[i], F[j])
        if !haskey(coeff,x)
            coeff[x] = discrete_bessel(x)
        end
        J[:, i, j] = coeff[x]
    end

    J
end

typealias QRfact Base.LinAlg.QRCompactWY{Complex128,Array{Complex128,2}}

"""
Stores all the J^n, the weights and the precomputed QR factorization.
"""
immutable BesselMatrix
    matrix::Array{Complex128,3}
    weights::Union{Void,Vector{Float64}}
    qr::Vector{QRfact}
end

BesselMatrix(E, x...) = BesselMatrix(E,E, x...)
function BesselMatrix(E::BispectralSet, F::BispectralSet; weights = nothing)
    @assert eltype(E) == eltype(F)
    @assert(weights == nothing || length(weights) == size(F,1))

    matrix = discrete_bessel_matrix(E,F)
    qr = Vector{QRfact}(size(matrix,1))
    for i in 1:size(matrix,1)
        # J = squeeze(matrix[i,:,:],1)
        J = matrix[i,:,:]
        if weights == nothing
            qr[i] = qrfact(J)
        else
            qr[i] = qrfact(J'*J + diagm(weights))
        end
    end

    BesselMatrix(matrix, weights, qr)
end
