module RKOpt

using Convex: MOI, Variable, evaluate, minimize, solve!
using ECOS: Optimizer


"""
    optimize_stability_polynomial

TODO

## References

- Ketcheson and Ahmadia (2012)
    Optimal stability polynomials for numerical integration of
    initial value problems
    [DOI: 10.2140/camcos.2012.7.247](https://doi.org/10.2140/camcos.2012.7.247)
"""
function optimize_stability_polynomial(
        accuracy_order,
        number_of_stages,
        spectrum;
        dt_min = 0.0,
        dt_max = 2.01 * number_of_stages^2 * maximum(abs, spectrum),
        tol_bisect = 1.0e-3,
        tol_feasible = 1.0e-9,
        maxiters = 1000)
    # TODO:
    # - Allow other bases
    # - Allow other solvers
    # - Allow function returning the spectrum for given dt

    Base.require_one_based_indexing(spectrum)
    if length(spectrum) <= number_of_stages
        throw(ArgumentError("The number of eigenvalues in `spectrum` must be greater than the `number_of_stages`."))
    end

    # Pre-compute powers of the eigenvalues for efficiency
    normalized_powers = similar(spectrum,
                                length(spectrum), number_of_stages)
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
        ## Terms given by free coefficients
        for j in (accuracy_order + 1):number_of_stages
            for i in eachindex(spectrum)
                polynomial_evaluations[i] += free_coefficients[j - accuracy_order] * normalized_powers_scaled[i, j]
            end
        end
        problem = minimize(maximum(abs, polynomial_evaluations))

        # Solve the convex optimization problem
        solve!(problem,
               MOI.OptimizerWithAttributes(Optimizer,
                                           "feastol" => tol_feasible);
                                           silent = true)

        # Update the time step size based on the optimization result
        if problem.optval < 1
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
    if number_of_stages > accuracy_order + 1
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


export optimize_stability_polynomial


end # module RKOpt
