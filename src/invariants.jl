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
