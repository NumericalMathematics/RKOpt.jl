using Test
using RKOpt

@testset "RKOpt" begin
    @testset "optimize_stability_polynomial" begin
        @testset "no free coefficients" begin
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
        end
    end
end
