using Test
using RKOpt

using Clarabel: Clarabel
using ECOS: ECOS

# We do not use
# - COSMO.Optimizer
# - ProxSDP.Optimizer
# - SCS.Optimizer
# since they were very slow in first tests
const OPTIMIZERS = [
    Clarabel.Optimizer,
    ECOS.Optimizer
]

@testset "RKOpt" begin
    @testset "optimize_stability_polynomial" begin
        @testset "no free coefficients" begin
            @testset "real coefficients" begin
                spectrum = range(-1.0, 0.0, length = 11)
                for n in 1:10
                    accuracy_order = n
                    number_of_stages = n
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    for i in 1:n
                        @test coefficients[i] ≈ 1 / factorial(i - 1)
                    end
                end

                # other optimizers
                @testset for optimizer in OPTIMIZERS
                    for n in 1:10
                        accuracy_order = n
                        number_of_stages = n
                        dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                                   number_of_stages,
                                                                                   spectrum;
                                                                                   optimizer = optimizer)
                        for i in 1:n
                            @test coefficients[i] ≈ 1 / factorial(i - 1)
                        end
                    end
                end
            end

            @testset "complex coefficients" begin
                spectrum = im * range(-1.0, 1.0, length = 11)
                for n in 1:10
                    accuracy_order = n
                    number_of_stages = n
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    for i in 1:n
                        @test coefficients[i] ≈ 1 / factorial(i - 1)
                    end
                end

                # other optimizers
                @testset for optimizer in OPTIMIZERS
                    for n in 1:10
                        accuracy_order = n
                        number_of_stages = n
                        dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                                   number_of_stages,
                                                                                   spectrum;
                                                                                   optimizer = optimizer)
                        for i in 1:n
                            @test coefficients[i] ≈ 1 / factorial(i - 1)
                        end
                    end
                end
            end
        end

        @testset "negative real axis monomial" begin
            # Section 4.1 of Ketcheson and Ahmadia (2012), page 259
            spectrum = range(-1.0, 0.0, length = 500)
            accuracy_order = 1
            for number_of_stages in 1:8
                dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                           number_of_stages,
                                                                           spectrum)
                @test isapprox(dt, 2 * number_of_stages^2, rtol = 0.001)
            end
        end

        @testset "imaginary axis monomial" begin
            # Section 4.1 of Ketcheson and Ahmadia (2012), page 260
            @testset "order 1" begin
                spectrum = im * range(-1.0, 1.0, length = 100)
                accuracy_order = 1
                for number_of_stages in 2:9
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
                end

                # other optimizers
                @testset for optimizer in OPTIMIZERS
                    for number_of_stages in 2:9
                        dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                                   number_of_stages,
                                                                                   spectrum;
                                                                                   optimizer = optimizer)
                        @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
                    end
                end

                # Compare with optimal polynomials given by
                # Kinnmark and Gray (1984)
                # One step integration methods with maximum stability regions
                # https://doi.org/10.1016/0378-4754(84)90039-9
                # Coefficients generated using Mathematica and the code
                #     OptimalStabilityPolynomialOrder1[K_, z_] :=
                #         Expand[(-I)^
                #         K*(I*ChebyshevT[K - 1, I*z/(K - 1)] - (1 + (z/(K - 1))^2)*
                #         ChebyshevU[K - 2, I*z/(K - 1)])]
                #
                #     OptimalStabilityPolynomialOrder1[2, z]
                #     OptimalStabilityPolynomialOrder1[3, z]
                #     OptimalStabilityPolynomialOrder1[4, z]
                #     OptimalStabilityPolynomialOrder1[5, z]
                #     OptimalStabilityPolynomialOrder1[6, z]
                #     OptimalStabilityPolynomialOrder1[7, z]
                #     OptimalStabilityPolynomialOrder1[8, z]
                #     OptimalStabilityPolynomialOrder1[9, z]
                # to obtain the monomial coefficients
                @testset let number_of_stages = 2
                    spectrum = im * range(-1.0, 1.0, length = 100)
                    accuracy_order = 1
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
                    @test isapprox(coefficients, [1, 1, 1], rtol = 0.001)
                end
                @testset let number_of_stages = 3
                    spectrum = im * range(-1.0, 1.0, length = 100)
                    accuracy_order = 1
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
                    optimal_coefficients = [1, 1, 1 / 2, 1 / 4]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
                @testset let number_of_stages = 4
                    spectrum = im * range(-1.0, 1.0, length = 500)
                    accuracy_order = 1
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
                    optimal_coefficients = [1, 1, 5 / 9, 4 / 27, 4 / 81]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
                @testset let number_of_stages = 5
                    spectrum = im * range(-1.0, 1.0, length = 500)
                    accuracy_order = 1
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
                    optimal_coefficients = [1, 1, 1 / 2, 3 / 16, 1 / 32, 1 / 128]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
                @testset let number_of_stages = 6
                    spectrum = im * range(-1.0, 1.0, length = 500)
                    accuracy_order = 1
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
                    optimal_coefficients = [
                        1,
                        1,
                        13 / 25,
                        4 / 25,
                        28 / 625,
                        16 / 3125,
                        16 / 15_625
                    ]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
                @testset let number_of_stages = 7
                    spectrum = im * range(-1.0, 1.0, length = 500)
                    accuracy_order = 1
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
                    optimal_coefficients = [
                        1,
                        1,
                        1 / 2,
                        19 / 108,
                        1 / 27,
                        2 / 243,
                        1 / 1458,
                        1 / 8748
                    ]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
                @testset let number_of_stages = 8
                    spectrum = im * range(-1.0, 1.0, length = 500)
                    accuracy_order = 1
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
                    optimal_coefficients = [
                        1,
                        1,
                        25 / 49,
                        8 / 49,
                        104 / 2401,
                        16 / 2401,
                        144 / 117_649,
                        64 / 823_543,
                        64 / 5_764_801
                    ]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
                @testset let number_of_stages = 9
                    spectrum = im * range(-1.0, 1.0, length = 500)
                    accuracy_order = 1
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
                    optimal_coefficients = [
                        1,
                        1,
                        1 / 2,
                        11 / 64,
                        5 / 128,
                        17 / 2048,
                        1 / 1024,
                        5 / 32_768,
                        1 / 131_072,
                        1 / 1_048_576
                    ]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
            end

            @testset "order 2, odd number of stages" begin
                spectrum = im * range(-1.0, 1.0, length = 100)
                accuracy_order = 2
                for number_of_stages in 3:2:9
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
                end

                # other optimizers
                @testset for optimizer in OPTIMIZERS
                    for number_of_stages in 3:2:9
                        dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                                   number_of_stages,
                                                                                   spectrum;
                                                                                   optimizer = optimizer)
                        @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
                    end
                end
            end

            @testset "order 2, even number of stages" begin
                spectrum = im * range(-1.0, 1.0, length = 100)
                accuracy_order = 2
                for number_of_stages in 4:2:10
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt,
                                   sqrt(number_of_stages * (number_of_stages - 2)),
                                   rtol = 0.001)
                end

                # other optimizers
                @testset for optimizer in OPTIMIZERS
                    for number_of_stages in 4:2:10
                        dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                                   number_of_stages,
                                                                                   spectrum;
                                                                                   optimizer = optimizer)
                        @test isapprox(dt,
                                       sqrt(number_of_stages * (number_of_stages - 2)),
                                       rtol = 0.001)
                    end
                end
            end

            @testset "order 3" begin
                # Compare with optimal polynomials given by
                # Kinnmark and Gray (1984)
                # One step integration methods with maximum stability regions
                # https://doi.org/10.1016/0378-4754(84)90039-9
                # Coefficients generated using Mathematica and the code
                #     OptimalStabilityPolynomialOrder2[K_, z_] :=
                #         Expand[Module[{SI = Sqrt[(K - 1)^2 - 1]},
                #         If[OddQ[K],
                #         1/(SI^2 + 1) +
                #             I^(K - 1)*SI^2/(SI^2 + 1)*ChebyshevT[K - 1, I*z/SI] +
                #             1/(SI^2 + 1)*z +
                #             1/2*I^(K + 2)*
                #             SI/(SI^2 + 1)*((K - 2)*ChebyshevT[K, I*z/SI] -
                #             K*ChebyshevT[K - 2, I*z/SI]),
                #         I^(K + 1)*Sqrt[SI^2/(SI^2 + 1)]*ChebyshevT[K - 1, I*z/SI] +
                #             1/2*I^K*Sqrt[
                #             1/(SI^2 + 1)]*((K - 2)*ChebyshevT[K, I*z/SI] -
                #             K*ChebyshevT[K - 2, I*z/SI])]
                #         ]]
                #
                #     OptimalStabilityPolynomialOrder2[3, z]
                #     OptimalStabilityPolynomialOrder2[4, z]
                #     OptimalStabilityPolynomialOrder2[5, z]
                #     OptimalStabilityPolynomialOrder2[6, z]
                #     OptimalStabilityPolynomialOrder2[7, z]
                #     OptimalStabilityPolynomialOrder2[8, z]
                #     OptimalStabilityPolynomialOrder2[9, z]
                # to obtain the monomial coefficients
                @testset let number_of_stages = 3
                    spectrum = im * range(-1.0, 1.0, length = 100)
                    accuracy_order = 3
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, sqrt(number_of_stages * (number_of_stages - 2)),
                                   rtol = 0.001)
                    optimal_coefficients = [1, 1, 1 / 2, 1 / 6]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
                @testset let number_of_stages = 4
                    spectrum = im * range(-1.0, 1.0, length = 100)
                    accuracy_order = 3
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, sqrt(number_of_stages * (number_of_stages - 2)),
                                   rtol = 0.001)
                    optimal_coefficients = [1, 1, 1 / 2, 1 / 6, 1 / 24]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
                @testset let number_of_stages = 5
                    spectrum = im * range(-1.0, 1.0, length = 100)
                    accuracy_order = 3
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test_broken isapprox(dt,
                                          sqrt(number_of_stages * (number_of_stages - 2)),
                                          rtol = 0.001)
                    optimal_coefficients = [1, 1, 1 / 2, 1 / 6, 1 / 30, 1 / 150]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
                @testset let number_of_stages = 6
                    spectrum = im * range(-1.0, 1.0, length = 500)
                    accuracy_order = 3
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, sqrt(number_of_stages * (number_of_stages - 2)),
                                   rtol = 0.001)
                    optimal_coefficients = [1, 1, 1 / 2, 1 / 6, 1 / 24, 1 / 180, 1 / 1080]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
                @testset let number_of_stages = 7
                    spectrum = im * range(-1.0, 1.0, length = 500)
                    accuracy_order = 3
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test_broken isapprox(dt,
                                          sqrt(number_of_stages * (number_of_stages - 2)),
                                          rtol = 0.001)
                    optimal_coefficients = [
                        1,
                        1,
                        1 / 2,
                        1 / 6,
                        4 / 105,
                        4 / 525,
                        8 / 11_025,
                        8 / 77_175
                    ]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
                @testset let number_of_stages = 8
                    spectrum = im * range(-1.0, 1.0, length = 500)
                    accuracy_order = 3
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, sqrt(number_of_stages * (number_of_stages - 2)),
                                   rtol = 0.001)
                    optimal_coefficients = [
                        1,
                        1,
                        1 / 2,
                        1 / 6,
                        1 / 24,
                        1 / 144,
                        1 / 864,
                        1 / 12_096,
                        1 / 96_768
                    ]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
                @testset let number_of_stages = 9
                    spectrum = im * range(-1.0, 1.0, length = 500)
                    accuracy_order = 3
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test_broken isapprox(dt,
                                          sqrt(number_of_stages * (number_of_stages - 2)),
                                          rtol = 0.001)
                    optimal_coefficients = [
                        1,
                        1,
                        1 / 2,
                        1 / 6,
                        5 / 126,
                        1 / 126,
                        4 / 3969,
                        4 / 27_783,
                        2 / 250_047,
                        2 / 2_250_423
                    ]
                    @test isapprox(coefficients, optimal_coefficients, rtol = 0.001)
                end
            end
        end

        @testset "disk monomial" begin
            # Section 4.1 of Ketcheson and Ahmadia (2012), page 262
            @testset "order 1" begin
                angle = range(0.0, 1.0 * π, length = 500)
                spectrum = @. -1.0 + cos(angle) + im * sin(angle)
                accuracy_order = 1
                for number_of_stages in 1:9
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    if number_of_stages == 2
                        @test_broken isapprox(dt, number_of_stages, rtol = 0.01)
                    else
                        @test isapprox(dt, number_of_stages, rtol = 0.01)
                    end
                end
            end

            @testset "order 2" begin
                angle = range(0.0, 1.0 * π, length = 500)
                spectrum = @. -1.0 + cos(angle) + im * sin(angle)
                accuracy_order = 2
                for number_of_stages in 2:9
                    dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                               number_of_stages,
                                                                               spectrum)
                    @test isapprox(dt, number_of_stages - 1, rtol = 0.01)
                end

                # other optimizers
                @testset for optimizer in OPTIMIZERS
                    for number_of_stages in 2:9
                        dt, coefficients = @inferred optimize_stability_polynomial(accuracy_order,
                                                                                   number_of_stages,
                                                                                   spectrum;
                                                                                   optimizer = optimizer)
                        @test isapprox(dt, number_of_stages - 1, rtol = 0.01)
                    end
                end
            end
        end

        @testset "error handling" begin
            @test_throws ArgumentError optimize_stability_polynomial(1, 1, [-1.0])
            @test_throws ArgumentError optimize_stability_polynomial(2,
                                                                     1,
                                                                     [-1.0, -0.5, 0.0])
        end
    end

    @testset "step_size_control_stability" begin
        # BS3
        coeff_main = [1, 1, 1 / 2, 1 / 6, 0]
        coeff_embd = [1, 1, 1 / 2, 3 / 16, 1 / 48]
        phi_ref = range(π / 2, π, length = 5)
        rad_ref = [0.7994394467275345, 0.9677459749859763, 0.8021332677908201, 0.7829788340996575, 0.8029999395923118]
        coeff_main = [1, 1, 1 / 2, 1 / 6, 0]
        phi, rad = @inferred step_size_control_stability(coeff_main, coeff_embd;
                                                         beta1 = 0.6, beta2 = -0.2, phi = phi_ref)
        @test isapprox(phi, phi_ref)
        @test isapprox(rad, rad_ref)
    end
end
