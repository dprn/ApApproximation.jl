dist(a::Angle, b::Angle)= min(abs(a-b), 2pi-abs(a-b))

function approx_in{N,T<:Real}(x::Frequency{N,T}, S::Set{Frequency{N,T}})
    for y in S
        approx_eq(x, y) && (return true)
    end
    return false
end

function addi{N,T<:Real}(a::Frequency{N,T},b::Frequency{N,T}, S::Set{Frequency{N,T}}, D::Set{Frequency{N,T}})
    temp=a.ω.val+a.ω.slice+atan2(b.λ*sin(b.ω-a.ω), a.λ+b.λ*cos(b.ω-a.ω))*N/(2*pi)
    w=Angle(temp, 0, N)
    r=sqrt(abs(a.λ^2+b.λ^2+2(a.λ*b.λ)*cos(convert(Float64,a.ω-b.ω))))
    u=Frequency(r,w)
    if approx_eq(Frequency(1.0, w),Frequency(1.0, 0., N)) || 1-w.val<TOL
        u=Frequency(r,0.,N)
    end
    v=normalize(u)
    if !approx_in(v, S) && !approx_in(v, D)
            push!(S, v)
    end
    S
end

function generate_set(N, n; verbose = false)#save_plot=false)
    D = Set{Frequency{N,Float64}}([Frequency(1.0, 0., N)])
    generate_set(D, n, verbose = verbose)
end

function generate_set{N}(D::Set{Frequency{N,Float64}}, n; verbose = false) #, save_plot=false)
  verbose && println("----------------------------------")
  verbose && println("    GENERATING BISPECTRAL SET")
  verbose && println("----------------------------------")
  AC = D
  NC = Set{Frequency{N,Float64}}()
  for iteration in 1:n
    verbose && println("Iteration $iteration started:")
    tic()
    for ac in AC, b in D
      for i in 0:N-1
        a=rotate(ac,i)
        if !approx_eq(a, b)
          addi(a, b, NC,D)
        end
      end
    end
    verbose && println("Elapsed time: $(toq())")
    AC = NC
    D = union(D, AC)
    verbose && println("Added $(length(AC)) frequencies")
    verbose && println("Total frequencies = $(length(D))")
    NC = Set{Frequency{N,Float64}}()
    # save_plot && png(plot(D),"plottings_$iteration")
    verbose && println("----------------------------------")
  end
  BispectralSet(collect(D))
end
