using Test
using RKOpt

@testset "RKOpt" begin
    @testset "optimize_stability_polynomial" begin
        @testset "no free coefficients" begin
            # real coefficients
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

            # complex coefficients
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
            end
        end

        @testset "error handling" begin
            @test_throws ArgumentError optimize_stability_polynomial(1, 1, [-1.0])
            @test_throws ArgumentError optimize_stability_polynomial(2,
                                                                     1,
                                                                     [-1.0, -0.5, 0.0])
        end
    end
end
