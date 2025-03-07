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
                dt, coefficients = @inferred optimize_stability_polynomial(
                    accuracy_order, number_of_stages, spectrum)
                for i in 1:n
                    @test coefficients[i] ≈ 1 / factorial(i - 1)
                end
            end

            # complex coefficients
            spectrum = im * range(-1.0, 1.0, length = 11)
            for n in 1:10
                accuracy_order = n
                number_of_stages = n
                dt, coefficients = @inferred optimize_stability_polynomial(
                    accuracy_order, number_of_stages, spectrum)
                for i in 1:n
                    @test coefficients[i] ≈ 1 / factorial(i - 1)
                end
            end
        end

        @testset "real axis monomial" begin
            spectrum = range(-1.0, 0.0, length = 500)
            accuracy_order = 1
            for number_of_stages in 1:8
                dt, coefficients = @inferred optimize_stability_polynomial(
                    accuracy_order, number_of_stages, spectrum)
                @test isapprox(dt, 2 * number_of_stages^2, rtol = 0.001)
            end
        end

        @testset "imaginary axis monomial" begin
            spectrum = im * range(-1.0, 1.0, length = 100)
            accuracy_order = 1
            for number_of_stages in 2:9
                dt, coefficients = @inferred optimize_stability_polynomial(
                    accuracy_order, number_of_stages, spectrum)
                @test isapprox(dt, number_of_stages - 1, rtol = 0.001)
            end
        end
    end
end
