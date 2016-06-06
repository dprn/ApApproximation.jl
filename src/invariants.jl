"""
Evaluates the power spectrum of the given BispInterpolation.
"""
PS(af::BispInterpolation) = PS(af.f)
PS(v)  = mapslices(x->dot(x,x), v, 2) |> real |> vec

"""
Evaluates the rotational power spectrum of the given BispInterpolation.
"""
RPS(af::BispInterpolation) = RPS(af.f)
function RPS(v)
  RPS = similar(v)
  for n in 1:size(v,2)
    RPS[:, n] = mapslices(x->dot(x, circshift(x, n)), v, 2)
  end
  RPS
end

"""
Returns a dictionary associating to each couple of integers (i,j), representing a
frequency in the first slice of E, the indexes (r,s) of the frequency E[i]+E[j],
if they exist.
"""
function bisp_split(E::BispectralSet)
  d = Dict{Tuple{Int,Int}, Tuple{Int,Int}}()
  for i in 1:size(E, 1), j in 1:size(E, 1)
    x = findin(E[i]+E[j], E)
    x!=(0,0) && (d[(i,j)] = x)
  end
  d
end

"""
Evaluates the bispectrum of the given BispInterpolation.
"""
function BS(af::BispInterpolation, d = bisp_split(af.E); max_rad = 10.)
  BS = Set{ Tuple{Float64,eltype(af)} }()

  for (key, val) in d
    if (af.E[key[1]] == zero(eltype(af.E))) && (af.E[key[2]] == zero(eltype(af.E)))
    # if af.E[key[1]] == af.E[key[2]] == zero(eltype(af.E))
      ω = squeeze(af[key[1],:], 1)
      push!(BS, (0. , dot(ω, ω)))
    else
      for n in 1:camembert(af)
        ω1 = squeeze(af[key[1],:], 1)
        ω2 = circshift(squeeze(af[key[2],:], 1), n)
        ω3 = circshift(squeeze(af[val[1],:], 1), -val[2]+n)
        push!(BS, (max(λ(af.E[key[1]]), λ(af.E[key[2]])) , real(dot(ω1.*ω2, ω3))))
      end
    end
  end
  sort(collect(BS), lt = (x,y)->x[1]<=y[1])
end


"""
Evaluates the rotational bispectrum of the given BispInterpolation.

NOT WORKING CORRECTLY
"""
function RBS(af::BispInterpolation, d = bisp_split(af.E))
  RBS = Vector{eltype(af)}()

  for (key, val) in d
    for n in 1:camembert(af)
      for k in 1:camembert(af)
        ω1 = circshift(squeeze(af[key[1],:], 1), k)
        ω2 = circshift(squeeze(af[key[2],:], 1), n)
        ω3 = circshift(squeeze(af[val[1],:], 1), -val[2]+n)
        push!(RBS, dot(ω1.*ω2, ω3))
      end
    end
  end
  RBS
end
