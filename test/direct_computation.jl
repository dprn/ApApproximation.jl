Base.exp(x::Frequency) =  exp(im*x.λ*cos(angle(x)))

function Base.exp(x::Frequency, y::Frequency)
    composition(x,y) |> exp
end

ev(af::BispInterpolation) = ev(af, af.E)
function ev{N,T}(af::BispInterpolation, F::BispectralSet{N,T})
    f = Array{Complex128, 2}(size(F)...)
    for n in 1:size(f,2), j in 1:size(f,1)
        f[j,n] = sum([ af.f[k,m]*exp(af.E[k,m],F[j,n]) for k in 1:size(af.E,1), m in 1:size(af.E,2) ])
    end
    BispInterpolation{N,Complex128}(f, F)
end
