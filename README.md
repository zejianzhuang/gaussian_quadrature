# Numerical integration using gaussian quandrature
The function defined in `GLQ.jl` can be used for calculating a integral with internal $[a, b]$, $[0, \infty)$, and $(-\infty, +\infty)$.

the syntax for `gaussian_quadrature` is as follows:
```julia
function gaussian_quadrature!(func::Function; 
                                interval::String=missing,  
                                a::Any=NaN, b::Any=NaN, 
                                digits::Int=10, return_maxiter::Bool=false)
```
The `interval` must be determined by one of these: `"(0, +inf)"`, `"(-inf, +inf)"` or `"Standard"`.

Example:
```julia
julia> gaussian_quadrature!(x -> x^2, interval="Standard", a=0., b=1.) 
julia> 0.3333333
```