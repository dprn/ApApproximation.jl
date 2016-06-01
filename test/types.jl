# findin test

v = Frequency{10, Float64}[Frequency(rand(), rand(), 10) for i in 1:100]
E = BispectralSet(v)

for i in 1:size(E,1), j in 1:size(E,2)
    @test findin(E[i,j],E) == (i,j)
end
