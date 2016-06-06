v = rand(Complex128, 100,10)
x = PS(v)

@assert length(x) == 100

for k in 1:length(x)
  @test x[k] == dot(squeeze(v[k,:],1),squeeze(v[k,:],1))
end

xx = RPS(v)
@assert size(xx) == size(v)
for k in 1:size(xx,1), n in 1:size(xx, 2)
  @test xx[k,n] == dot(squeeze(v[k,:],1),circshift(squeeze(v[k,:],1),n))
end
