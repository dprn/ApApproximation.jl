import Base: fft, ifft, real
real{N,T<:Real}(f::BispInterpolation{N,Complex{T}}) = BispInterpolation{N,T}(real(f.f), f.E)


"""
Appalingly, this is faster than the equivalent:
radial_fft{N,T<:Number}(f::BispInterpolation{N,T}) = fft(f.f, 2)
"""
function radial_fft{N,T<:Number}(f::BispInterpolation{N,T})
    x = similar(f.f, Complex128)
    for i in 1:size(x,1)
        x[i,:] = fft(f.f[i,:])
    end
    x
end

"""
Same here w.r.t.
radial_ifft{T<:Number}(hat_af::Array{T, 2}, E::BispectralSet) = BispInterpolation{camembert(E), Complex128}(ifft(hat_af, 2), E)
"""
function radial_ifft{T<:Number}(hat_af::Array{T, 2}, E::BispectralSet)
    x = similar(hat_af, Complex128)
    for i in 1:size(x,1)
        x[i,:] = ifft(hat_af[i,:])
    end
    BispInterpolation{camembert(E), Complex128}(x, E)
end

"""
Actual AP approximation of a function defined on the bispectral set

- `E` : The BispectralSet of definition of the image
- `F` : The BispectralSet of the frequencies
- `J` : The precomputed Bessel matrices. If it is not given, we compute it,
        assuming the weights to be all fixed at one.
"""
ap(f, J = nothing) = ap(f, f.E, J)
function ap{R<:Real, T<:Real, N}(f::BispInterpolation{N, R}, F::BispectralSet{N,T}, J::Union{Void, BesselMatrix} = nothing)
    # Check if we are given the precomputed Bessel matrices
    J == nothing && (matrix = discrete_bessel_matrix(f.E,F))

    hat_f = radial_fft(f)

    hat_af = Array(Complex128, size(F))
    for n = 1:N
        if J == nothing
            A = squeeze(matrix[mod1(n-1,N),:,:],1)
            b = hat_f[:,n]
        else
            A = J.qr[mod1(n-1,N)]
            j = squeeze(J.matrix[mod1(n-1,N),:,:],1)
            A = (J.weights == nothing) ? j : j'*j+diagm(J.weights) #squeeze(J.matrix[mod1(n-1,N),:,:],1)
            b = (J.weights == nothing) ? hat_f[:,n] : j' * hat_f[:,n]
        end

        hat_af[:,n] = A \ b
    end

    af = radial_ifft(hat_af, F)
end

"""
Inverse AP approximation

- `af` : The values of the pseudo-periodic approximation of the image
		on the BispectralSet.
- `E` : The BispectralSet of definition of the image
- `F` : The BispectralSet of the frequencies
"""
iap(af, J = nothing) = iap(af, af.E, J)
function iap{R<:Number, T<:Real, N}(af::BispInterpolation{N, R}, E::BispectralSet{N,T}, J::Union{Void, BesselMatrix} = nothing)
    # Check if we are given the precomputed Bessel matrices
    matrix = (J == nothing) ? discrete_bessel_matrix(E,af.E) : J.matrix

    hat_af = radial_fft(af)

    hat_f = Array(Complex128, size(E))
    for n = 1:N
        hat_f[:,n] = squeeze(matrix[mod1(n-1,N),:,:],1)*hat_af[:,n]
    end

    f = radial_ifft(hat_f, E)
    real(f)
end
