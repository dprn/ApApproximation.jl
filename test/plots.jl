using Plots
gadfly()

function Plots.plot{T<:Real, N}(E::BispectralSet{N,T})
    x, y = cart(E)

    x_max, y_max = maximum(x), maximum(y)

    p = scatter(x, y, marker = (:red, 4))
    ang = 0
    for i in 1:N
        if abs(ang-pi/2)>1e-10
            f(x) = tan(ang)*x
            lim = 1.3*min(x_max, y_max/abs(tan(ang)))
            plot!(f, -lim:.1:lim, line = (:black, :dot), legend = false)
        else
            plot!(zeros(100), linspace(-1.3*y_max, 1.3*y_max, 100))
        end
        ang += 2pi/N
    end
    p
end
