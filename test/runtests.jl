using Toro
using Test

@testset "Toro.jl" begin

    # Test the output of CenterPressure, which returns the value of the middle pressure, P₀
    function TestCenterPressure()
        # Define our test cases (Toro, Table 4.1)
        # Order of inputs is γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ, P₀ (expected)
        test_cases = [
                      (  1.400,     1.000,  1.000,     0.000,   1.400,    0.100,  0.125,    0.000,      0.30313 ),
                      (  1.400,     0.400,  1.000,    -2.000,   1.400,    0.400,  1.000,    2.000,      0.00189 ),
                      (  1.400,  1000.000,  1.000,     0.000,   1.400,    0.010,  1.000,    0.000,    460.894   ),
                      (  1.400,     0.010,  1.000,     0.000,   1.400,  100.000,  1.000,    0.000,     46.0950  ),
                      (  1.400,   460.894,  5.99924,  19.5975,  1.400,   46.095,  5.99242, -6.19633, 1691.64    ),
                     ]
        
        absolute_tolerance = 0.01
        @testset "Center Pressure Tests" begin
            for ( γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ, ExpectedValue ) in test_cases
                # First we need to calculate the number of expected digits of precision based on the supplied answers from Toro Table 4.3
                # We can't use a single tolerance because the number of reported digits is different in each case

                # First, convert the expected value to a string, break it on the decimal place, and count the number of digits after the decimal
                ExpectedDigits = length(split(string(ExpectedValue),".")[2])
                # Then convert that number of digits to an absolute tolerance
                ExpectedTolerance = parse(Float64,"1e-$ExpectedDigits")
                # Finally, run the test
                ActualValue = CenterPressure( γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ )
                # And compare the answers using the computed tolerance
                @test isapprox( ExpectedValue, ActualValue, atol=ExpectedTolerance )
            end
        end
    end
    TestCenterPressure()

    # Test the output of CenterVelocity, which returns the value of the middle velocity, u₀
    function TestCenterVelocity()
        # Define our test cases (Toro, Table 4.1)
        # Order of inputs is P₀, γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ, u₀ (expected)
        test_cases = [
                      (    0.30313, 1.400,     1.000,  1.000,     0.000,   1.400,    0.100,   0.125,    0.000,   0.92745 ),
                      (    0.00189, 1.400,     0.400,  1.000,    -2.000,   1.400,    0.400,   1.000,    2.000,   0.00000 ),
                      (  460.894,   1.400,  1000.000,  1.000,     0.000,   1.400,    0.010,   1.000,    0.000,  19.5975  ),
                      (   46.0950,  1.400,     0.010,  1.000,     0.000,   1.400,  100.000,   1.000,    0.000,  -6.19633 ),
                      ( 1691.64,    1.400,   460.894,  5.99924,  19.5975,  1.400,   46.095,   5.99242, -6.19633, 8.6897  ), 
                     ]
        # NOTE: u₀ (expected) in the 5th test has the 5th digit removed (u₀(expected)=8.68975 in table 4.3). 
        # This was done as this value can not be reproduced beyond 4 digits, either in this code or via Wolfram Alpha, etc.
        # Indeed, running Wolfram out to >32 digits produces a value that matches Julia's result exactly
        # This could potentially indicate a typo in the table, or numerical noise if the code was 32 bit, or etc
        
        absolute_tolerance = 0.01
        @testset "Center Velocity Tests" begin
            for ( P₀, γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ, ExpectedValue ) in test_cases
                # First we need to calculate the number of expected digits of precision based on the supplied answers from Toro Table 4.3
                # We can't use a single tolerance because the number of reported digits is different in each case

                # First, convert the expected value to a string, break it on the decimal place, and count the number of digits after the decimal
                ExpectedDigits = length(split(string(ExpectedValue),".")[2])
                # Then convert that number of digits to an absolute tolerance
                ExpectedTolerance = parse(Float64,"1e-$ExpectedDigits")
                # Finally, run the test
                #P₀ = CenterPressure( γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ )
                ActualValue = CenterVelocity( P₀, γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ )
                # And compare the answers using the computed tolerance
                @test isapprox( ExpectedValue, ActualValue, atol=ExpectedTolerance )
            end
        end
    end
    TestCenterVelocity()
    #
    # Test the output of EdgeDensity, which returns the value of the density in the edge region, ρₗ (or ρᵣ)
    function TestEdgeDensity()
        # Define our test cases (Toro, Table 4.1)
        # Order of inputs is P₀, γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ, ρₗ (expected), ρᵣ (expected)
        test_cases = [
                      (    0.30313, 1.400,     1.000,  1.000,     0.000,   1.400,    0.100,   0.125,    0.000,    0.42632,  0.26557 ),
                      (    0.00189, 1.400,     0.400,  1.000,    -2.000,   1.400,    0.400,   1.000,    2.000,    0.0218,   0.0218  ),
                      (  460.894,   1.400,  1000.000,  1.000,     0.000,   1.400,    0.010,   1.000,    0.000,    0.57506,  5.99924 ),
                      (   46.0950,  1.400,     0.010,  1.000,     0.000,   1.400,  100.000,   1.000,    0.000,    5.99242,  0.57511 ),
                      ( 1691.64,    1.400,   460.894,  5.99924,  19.5975,  1.400,   46.095,   5.99242, -6.19633, 14.2823,  31.0426  ),
                     ]
        # NOTE: ρₗ(expected) and ρᵣ (expected) in the 2nd test has the 5th digit removed (ρₗ(expected)=ρᵣ(expected)=0.02185 in table 4.3). 
        # This was done as this value can not be reproduced beyond 4 digits, either in this code or via Wolfram Alpha, etc., with the expected value of P₀ supplied in the table
        # Unlike the velocity test modification above, the expected result for ρₗ or ρᵣ in table 4.3 *can* be reproduced when using P₀ calculated from CenterPressure() (or from Wolfram Alpha)
        # Thus, the failure to reproduce the tabulated result beyond 4 digits is likely due to the recorded precision of P₀
        # There are two fixes - either use the output of CenterPressure() in the test, or reduce the expected precision of the result
        # The latter was opted for in this case as we're still comparing expected results with tabulated data with this approach
        
        absolute_tolerance = 0.01
        @testset "Edge Density Tests" begin
            for ( P₀, γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ, ExpectedValueLeft, ExpectedValueRight ) in test_cases
                # First we need to calculate the number of expected digits of precision based on the supplied answers from Toro Table 4.3
                # We can't use a single tolerance because the number of reported digits is different in each case

                # Need to establish the tolerance for our answers based on the precision of the data provided in the table
                # First, convert the expected value to a string, break it on the decimal place, and count the number of digits after the decimal
                ExpectedDigitsLeft = length(split(string(ExpectedValueLeft),".")[2])
                ExpectedDigitsRight = length(split(string(ExpectedValueRight),".")[2])
                # Then convert that number of digits to an absolute tolerance
                ExpectedToleranceLeft = parse(Float64,"1e-$ExpectedDigitsLeft")
                ExpectedToleranceRight = parse(Float64,"1e-$ExpectedDigitsRight")
                # Finally, run the test
                #P₀ = CenterPressure( γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ )
                # Left edge
                ActualValueLeft = EdgeDensity( P₀, γₗ, Pₗ, ρₗ, uₗ )
                @test isapprox( ExpectedValueLeft, ActualValueLeft, atol=ExpectedToleranceLeft )
                # Right edge
                ActualValueRight = EdgeDensity( P₀, γᵣ, Pᵣ, ρᵣ, uᵣ )
                @test isapprox( ExpectedValueRight, ActualValueRight, atol=ExpectedToleranceRight )
            end
        end
    end
    TestEdgeDensity()
end
