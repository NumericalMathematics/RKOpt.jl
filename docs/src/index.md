# RKOpt.jl

The [Julia](https://julialang.org) package
[RKOpt.jl](https://github.com/NumericalMathematics/RKOpt.jl)
is work in progress. A mature version of much more functionality
is available in the MATLAB-based package
[RK-Opt](https://github.com/ketch/RK-Opt).


## Example: Optimize a stability polynomial for imaginary axis inclusion

First, we optimize the stability polynomial for imaginary
axis inclusion using [`optimize_stability_polynomial`](@ref).

```@example polyopt-imaginary-axis
using RKOpt: optimize_stability_polynomial

spectrum = im * range(0.0, 1.0, length = 101)
accuracy_order = 1
number_of_stages = 7
dt, coefficients = optimize_stability_polynomial(accuracy_order, number_of_stages, spectrum)
```

Next, we plot the corresponding stability region.

```@example polyopt-imaginary-axis
using GLMakie
using LaTeXStrings

n = 1000
x = range(-6, 1, length = n)
y = range(-7, 7, length = 2 * n)
z = @. x + im * y'
stab = map(z -> abs(evalpoly(z, coefficients)) <= 1, z)
fig = Figure()
ax = Axis(fig[1, 1]; aspect = DataAspect(),
          xlabel = L"\mathrm{Re}(z)", ylabel = L"\mathrm{Im}(z)")
heatmap!(ax, x, y, stab, colormap = :binary, colorrange = (0.0, 2.0))
lines!(ax, [extrema(x)...], [0, 0], color = :black, linestyle = :dash)
lines!(ax, [0, 0], [extrema(y)...], color = :black, linestyle = :dash)
fig
```
