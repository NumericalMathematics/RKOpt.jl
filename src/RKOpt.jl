module RKOpt

using LinearAlgebra: eigvals

using Convex: MOI, Variable, evaluate, minimize, solve!
using ECOS: ECOS

"""
    optimize_stability_polynomial(accuracy_order,
                                  number_of_stages,
                                  spectrum;
                                  dt_min = 0.0,
                                  dt_max = 2.01 * number_of_stages^2 *
                                           maximum(abs, spectrum),
                                  tol_bisect = 1.0e-9,
                                  tol_feasible = 1.0e-9,
                                  maxiters = 1000,
                                  optimizer = ECOS.Optimizer,
                                  silent = true)

Optimize the stability polynomial of an explicit Runge-Kutta method
with a given `accuracy_order` and `number_of_stages` for a given
`spectrum` of eigenvalues (a vector of real or complex numbers) using
the algorithm of Ketcheson and Ahmadia (2012).

## Optional keyword arguments

- `dt_min = 0.0`:
  Minimum time step size for the bisection.
- `dt_max = 2.01 * number_of_stages^2 * maximum(abs, spectrum)`:
  Maximum time step size for the bisection.
- `tol_bisect = 1.0e-9`:
  Relative tolerance used to terminate the bisection.
- `tol_feasible = 1.0e-9`:
  Relative tolerance used to determine the feasibility of the optimization.
  The problem is considered feasible if the maximal absolute value of the
  stability polynomial is less than or equal to `1 + tol_feasible` at all
  eigenvalues in the `spectrum`.
- `maxiters = 1000`:
  Maximum number of iterations for the bisection.
- `optimizer = ECOS.Optimizer`:
  Convex optimization solver to be used. You can choose a solver supporting
  SOCP from the list of supported solvers in the
  [JuMP documentation](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers).
- `silent = true`:
  Whether to suppress the output of the convex optimization solver.

## References

- Ketcheson and Ahmadia (2012)
  Optimal stability polynomials for numerical integration of
  initial value problems
  [DOI: 10.2140/camcos.2012.7.247](https://doi.org/10.2140/camcos.2012.7.247)
"""
function optimize_stability_polynomial(accuracy_order,
                                       number_of_stages,
                                       spectrum;
                                       dt_min = 0.0,
                                       dt_max = 2.01 * number_of_stages^2 *
                                                maximum(abs, spectrum),
                                       tol_bisect = 1.0e-9,
                                       tol_feasible = 1.0e-9,
                                       maxiters = 1000,
                                       optimizer = ECOS.Optimizer,
                                       silent = true)
    # TODO:
    # - Allow other bases
    # - Allow function returning the spectrum for given dt

    Base.require_one_based_indexing(spectrum)
    if length(spectrum) <= number_of_stages
        throw(ArgumentError("The number of eigenvalues in `spectrum` must be greater than the `number_of_stages`."))
    end
    if accuracy_order > number_of_stages
        throw(ArgumentError("The `accuracy_order` must be less than or equal to the `number_of_stages`."))
    end

    # Pre-compute powers of the eigenvalues for efficiency
    normalized_powers = similar(spectrum, length(spectrum), number_of_stages)
    for j in 1:number_of_stages
        factorial_j = factorial(j)
        for i in eachindex(spectrum)
            normalized_powers[i, j] = spectrum[i]^j / factorial_j
        end
    end
    normalized_powers_scaled = similar(normalized_powers)

    # Initialize variable for polynomial coefficients
    free_coefficients = Variable(number_of_stages - accuracy_order)

    # Allocate memory for the evaluation of the stability polynomial
    # at the eigenvalues
    polynomial_evaluations = similar(spectrum)

    # Start bisection based on the time step size
    dt = zero(dt_max)
    for iter in 1:maxiters
        dt = 0.5 * (dt_min + dt_max)

        # Compute the stability polynomial for the current time step size
        for j in 1:number_of_stages
            dt_j = dt^j
            for i in eachindex(spectrum)
                normalized_powers_scaled[i, j] = dt_j * normalized_powers[i, j]
            end
        end

        # Construct the convex optimization problem
        ## Zeroth-order term of the polynomials
        fill!(polynomial_evaluations, 1)
        ## Terms fixed by the accuracy order
        for j in 1:accuracy_order
            for i in eachindex(spectrum)
                polynomial_evaluations[i] += normalized_powers_scaled[i, j]
            end
        end

        if number_of_stages == accuracy_order
            # If there are no free coefficients, we do not need to perform
            # the optimization but just compute the value
            max_abs_polynomial_evaluations = maximum(abs, polynomial_evaluations)
        else
            # We need to solve the convex optimization problem

            ## We cannot add the terms given by free coefficients in-place
            ## since they have types coming from Convex.jl.
            # for j in (accuracy_order + 1):number_of_stages
            #     for i in eachindex(spectrum)
            #         polynomial_evaluations[i] += free_coefficients[j - accuracy_order] * normalized_powers_scaled[i, j]
            #     end
            # end
            all_polynomial_evaluations = polynomial_evaluations
            for j in (accuracy_order + 1):number_of_stages
                all_polynomial_evaluations += free_coefficients[j - accuracy_order] *
                                              normalized_powers_scaled[:, j]
            end
            problem = minimize(maximum(abs(all_polynomial_evaluations)))

            # Solve the convex optimization problem
            solve!(problem, optimizer; silent = silent)
            max_abs_polynomial_evaluations = problem.optval
        end

        # Update the time step size based on the optimization result
        if max_abs_polynomial_evaluations <= 1 + tol_feasible
            dt_min = dt
        else
            dt_max = dt
        end

        # Terminate iteration when the relative tolerance is reached
        if (dt_max < tol_bisect) || (dt_max - dt_min < tol_bisect * dt_min)
            break
        end

        if iter == maxiters
            @warn "Maximum number of iterations reached."
        end
    end

    # Collect all polynomial coefficients after the optimization
    coefficients = ones(number_of_stages + 1)
    if number_of_stages > accuracy_order
        free_coefficients_opt = evaluate(free_coefficients)
        for i in eachindex(free_coefficients_opt)
            coefficients[accuracy_order + 1 + i] = free_coefficients_opt[i]
        end
    end

    # Due to the normalization of the powers of the eigenvalues,
    # the coefficients are effectively scaled by factorials.
    # We need to remove this scaling to obtain the actual coefficients.
    for i in eachindex(coefficients)
        coefficients[i] /= factorial(i - 1)
    end

    return dt, coefficients
end


function step_size_control_stability(stability_polynomial_coefficients_main,
                                     stability_polynomial_coefficients_embedded;
                                     beta1 = 0.6,
                                     beta2 = -0.2,
                                     radius_min = 0.1,
                                     radius_max = 2 * length(stability_polynomial_coefficients_main),
                                     tol = 1.0e-12,
                                     phi = range(π / 2, π, length = 200))
    Base.require_one_based_indexing(stability_polynomial_coefficients_main,
                                    stability_polynomial_coefficients_embedded)
    if length(stability_polynomial_coefficients_main) != length(stability_polynomial_coefficients_embedded)
        throw(DimensionMismatch("The lengths of the main and embedded stability polynomial coefficients must be the same."))
    end

    # R(z)
    coeff_main = stability_polynomial_coefficients_main
    # R'(z) z
    deriv_main = [coeff_main[i] * (i - 1) for i in eachindex(coeff_main)]
    # Rhat(z)
    coeff_embd = stability_polynomial_coefficients_embedded
    # E(z) = R(z) - Rhat(z)
    coeff_diff = coeff_main - coeff_embd
    # E'(z) z
    deriv_diff = [coeff_diff[i] * (i - 1) for i in eachindex(coeff_diff)]

    # Check order of accuracy of the stability polynomials
    p_main = length(coeff_main) - 1
    for i in eachindex(coeff_main)
        if !(coeff_main[i] * factorial(i - 1) ≈ 1)
            p_main = i - 2
            break
        end
    end
    p_embd = length(coeff_embd) - 1
    for i in eachindex(coeff_embd)
        if !(coeff_embd[i] * factorial(i - 1) ≈ 1)
            p_embd = i - 2
            break
        end
    end
    @show p_main, p_embd

    # Compute the boundary of the stability region
    # for a given angle/direction in the complex plane
    stability_polynomial(z) = evalpoly(z, coeff_main)
    stability_function_squared(r, phi) = abs2(stability_polynomial(r * cis(phi)))
    function compute_z(phi, r_min, r_max, tol)
        while r_max - r_min > tol
            r = (r_min + r_max) / 2
            if stability_function_squared(r, phi) > 1
                r_max = r
            else
                r_min = r
            end
        end
        z = r_min * cis(phi)
        return z
    end

    # Compute the spectral radius of the Jacobian matrix
    # determining step size control stability for a PI controller
    function spectral_radius(phi)
        z = compute_z(phi, radius_min, radius_max, tol)
        k = min(p_main, p_embd) + 1

        # R(z): main stability polynomial
        # u = Re( R'(z) * z / R(z))
        u = real(evalpoly(z, deriv_main) / evalpoly(z, coeff_main))
        # Rhat(z): embedded stability polynomial
        # E(z) = R(z) - Rhat(z)
        # v = Re( E'(z) * z / E(z))
        v = real(evalpoly(z, deriv_diff) / evalpoly(z, coeff_diff))
        jacobian = [1 u 0 0;
                    (-beta1 / k) (1 - v * beta1 / k) (-beta2 / k) (-v * beta2 / k);
                    1 0 0 0;
                    0 1 0 0]
        λ = eigvals(jacobian)
        return maximum(abs, λ)
    end

    rad = map(spectral_radius, phi)
    return phi, rad
end

export optimize_stability_polynomial
export step_size_control_stability

end # module RKOpt
