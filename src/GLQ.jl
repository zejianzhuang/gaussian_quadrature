include("./xw.jl")

function gaussian_quadrature!(func::Function; interval::String=missing,  a::Any=NaN, b::Any=NaN, digits::Int=10, return_maxiter::Bool=false)
    maxiter::Int = 4
    EPS = 1e-34
    epsabs = 1e-8
    Δ = epsabs * 1000
    last_inte = 0.0
    i::Int = 1
    while abs(Δ) > epsabs
        x, w = xw[i]

        if interval == "Standard" && a != NaN && b != NaN
            xp, wp = @. (b - a) * x / 2.0 + (a + b) / 2.0, (b - a) / 2.0 * w
        elseif interval == "(0, +inf)"
            xp, wp = @. (1.0 + x) / (1.0 - x + EPS), 2 / (x - 1.0 + EPS)^2 * w
        elseif interval == "(-inf, +inf)"
            xp, wp = @. (1.0 - x^2 + EPS), (1.0 + x^2) / (-1.0 + x^2 + EPS)^2 * w
        else
            print("Invalid keyword to determine the integral interval. keyword to determine an interval should be \"(0, +inf)\", \"(-inf, +inf)\" or \"Standard\".")
            break
        end

        inte = wp' * func.(xp) # Inner product
        Δ = inte - last_inte
        last_inte = inte
        maxiter += 1
        i += 1

        if length(xw) == i
            break
        end
    end

    if return_maxiter == true
        return round(last_inte, digits=digits), maxiter
    else
        return round(last_inte, digits=digits)
    end

end
