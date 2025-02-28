# These functions implement the Riemann solver functions from Toro Ch 4
using ShockWaves

"""
    toro_constant_A( γ::T, P::T, ρ::T ) where { T <: AbstractFloat }

Computes and returns the value of the constant "A" in Eq. 4.8

# Arguments
    γ: The ratio of specific heats for the gas
    P: The pressure of the gas
    ρ: The density of the gas

# Returns
    A value of type T representing the evaluation of the constant "A" in Eq. 4.8

# Equation
    This function evaluated the equation
        2 / ( (γ+1) * ρ )
    where
        γ is the ratio of specific heats
        ρ is the density

# Notes
    The parameter P is passed to this function but is not used. This is done to keep the function signatures for the constants A and B identical.
"""
function toro_constant_A( γ::T, P::T, ρ::T ) where { T <: AbstractFloat }
    return 2.0 / ( ( γ + 1.0 ) * ρ )
end

"""
    toro_constant_B( γ::T, P::T, ρ::T ) where { T <: AbstractFloat }

Computes and returns the value of the constant "B" in Eq. 4.8

# Arguments
    γ: The ratio of specific heats for the gas
    P: The pressure of the gas
    ρ: The density of the gas

# Returns
    A value of type T representing the evaluation of the constant "B" in Eq. 4.8

# Equation
    This function evaluated the equation
        ( γ-1 )/( γ+1 ) * P
    where
        γ is the ratio of specific heats
        P is the pressure

# Notes
    The parameter ρ is passed to this function but is not used. This is done to keep the function signatures for the constants A and B identical.
"""
function toro_constant_B( γ::T, P::T, ρ::T ) where { T <: AbstractFloat }
    return ( γ - 1.0 ) * P / ( γ + 1.0 )
end

"""
    toro_fₗ( P₀::T, γ::T, P::T, ρ::T ) where { T <: AbstractFloat }
    toro_fᵣ( P₀::T, γ::T, P::T, ρ::T ) where { T <: AbstractFloat }

Computes and returns the value of the function 𝒻 (Eqs. 4.6 and 4.7)

# Arguments
    P₀: The test pressure
    γ : Ratio of specific heats for the gas
    P : The pressure of the gas
    ρ : The density of the gas

# Returns 
    A value of type T representing the evaluation of Eqs. 4.6 and/or 4.7

# Equation
    Per Proposition 4.2.1, the function 𝒻 should return
        ( P₀ - P ) * √( A / ( P + B ) )                             { if P₀ > P (Shock) }
    or
        ( 2 a ) / ( γ - 1 ) * [ ( P₀/P )^((γ-1)/2γ) - 1 ]           { if P₀ < P (Rarefaction) } 
    where
        P₀ is the test pressure
        P is the gas pressure
        γ is the ratio of specific heats
        a = √(γ P/ρ ) is the speed of sound
        A is the output of toro_constant_A( γ, P, ρ ) (Eq. 4.8)
        B is the output of toro_constant_B( γ, P, ρ ) (Eq. 4.8)

# Notes
    toro_fₗ and toro_fᵣ are identical functions. In fact, toro_fᵣ is simply aliased to toro_fₗ.
    These functions are impelmented this way to simplify function structure and make it more clear how the functions align with the equations in the text.
    In the text, there is a function f that is to be minimized, and it consists of functions fᵣ and fₗ, which is the reason for this approach.
"""
function toro_fₗ( P₀::T, γ::T, P::T, ρ::T ) where { T <: AbstractFloat }
    A = toro_constant_A( γ, P, ρ )
    B = toro_constant_B( γ, P, ρ )
    if ( P₀ > P )
        return ( P₀ - P ) * sqrt( A / ( P₀ + B ) )
    else
        return ( 2.0 * sqrt( γ * P / ρ ) ) / ( γ - 1.0 ) * ( ( P₀ / P )^( ( γ - 1.0 ) / ( 2.0 * γ ) ) - 1.0 )
    end
end

"""
    toro_fₗ( P₀::T, γ::T, P::T, ρ::T ) where { T <: AbstractFloat }
    toro_fᵣ( P₀::T, γ::T, P::T, ρ::T ) where { T <: AbstractFloat }

Computes and returns the value of the function 𝒻ₗ or fᵣ (Eqs. 4.6 and 4.7)

# Arguments
    P₀: The test pressure
    γ : Ratio of specific heats for the gas
    P : The pressure of the gas
    ρ : The density of the gas

# Returns 
    A value of type T representing the evaluation of Eqs. 4.6 and/or 4.7

# Equation
    Per Proposition 4.2.1, the function 𝒻ₗ or fᵣ should return
        ( P₀ - P ) * √( A / ( P + B ) )                             { if P₀ > P (Shock) }
    or
        ( 2 a ) / ( γ - 1 ) * [ ( P₀/P )^((γ-1)/2γ) - 1 ]           { if P₀ < P (Rarefaction) } 
    where
        P₀ is the test pressure
        P is the gas pressure
        γ is the ratio of specific heats
        a = √(γ P/ρ ) is the speed of sound
        A is the output of toro_constant_A( γ, P, ρ ) (Eq. 4.8)
        B is the output of toro_constant_B( γ, P, ρ ) (Eq. 4.8)

# Notes
    toro_fₗ and toro_fᵣ are identical functions. In fact, toro_fᵣ is simply aliased to toro_fₗ.
    These functions are impelmented this way to simplify function structure and make it more clear how the functions align with the equations in the text.
    In the text, there is a function f that is to be minimized, and it consists of functions fᵣ and fₗ, which is the reason for this approach.
"""
const toro_fᵣ = toro_fₗ

"""
    toro_f( P₀::T, γₗ::T, Pₗ::T, ρₗ::T, uₗ::T, γᵣ::T, Pᵣ::T, ρᵣ::T, uᵣ::T ) where { T <: AbstractFloat }

Computes and returns the value of the function 𝒻 (Eq. 4.5)

# Arguments
    P₀: The test pressure
    γ : Ratio of specific heats for the gas
    P : The pressure of the gas
    ρ : The density of the gas

    ⋅ Subscripts indicate whether the property is associated with the gas to the left (𝒻ₗ) or right (𝒻ᵣ) of the interface

# Returns 
    A value of type T representing the evaluation of Eq. 4.5

# Equation
    Eq. 4.5 is defined as
        𝒻( P₀, Wₗ, Wᵣ ) = fₗ( P₀, Wₗ ) + fᵣ( P₀, Wᵣ ) + ( uᵣ - uₗ )
    where
        P₀ is the test pressure
        W = [γ, P, ρ] is the state associated with the left (Wₗ) or right (Wᵣ) of the interface
"""
function toro_f( P₀::T, γₗ::T, Pₗ::T, ρₗ::T, uₗ::T, γᵣ::T, Pᵣ::T, ρᵣ::T, uᵣ::T ) where { T <: AbstractFloat }
    return toro_fₗ( P₀, γₗ, Pₗ, ρₗ ) + toro_fᵣ( P₀, γᵣ, Pᵣ, ρᵣ ) + uᵣ - uₗ
end

"""
    CenterPressure( γₗ::T, Pₗ::T, ρₗ::T, uₗ::T, γᵣ::T, Pᵣ::T, ρᵣ::T, uᵣ::T ) where { T <: AbstractFloat }

Finds the value of the center pressure, P₀ that solves Eq. 4.5 using Eqs 4.6, 4.7, 4.8, and 4.9

# Arguments
    γ: The ratio of specific heats of the gas
    P: The pressure of the gas
    ρ: The density of the gas
    u: The velocity of the gas

    ⋅ Subscripts indicate whether the property is associated with the gas to the left (𝒻ₗ) or right (𝒻ᵣ) of the interface

# Returns
    A value of type T representing the value of P₀ that satisfies Eq. 4.5

# Notes
    As an implementation detail, the root finding algorithm is currently a simple bisection method
"""
function CenterPressure( γₗ::T, Pₗ::T, ρₗ::T, uₗ::T, γᵣ::T, Pᵣ::T, ρᵣ::T, uᵣ::T ) where { T <: AbstractFloat }
    # We need to start off by getting the minimum and maximum pressures
    P₋ = min( Pₗ, Pᵣ ) # The lesser of the two pressures
    P₊ = max( Pₗ, Pᵣ ) # The greater of the two pressures
    
    # Evaluate the function f at the minimum and maximum pressures
    f₋ = toro_f( P₋, γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ )
    f₊ = toro_f( P₊, γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ )

    # Per Eq. 4.39, there are three ranges for the solution depending on the velocity difference
    # According to section 4.3.1, the function f is monotone and concave down
    # Therefore, depending on the sign of f at the minimum and maximum pressures, we can conclude where the zero crossing must lie
    # Here, we check those 3 possibilities and assign our initial lower (P₋) and upper (P₊) guesses
    if ( ( f₋ > 0.0 ) & ( f₊ > 0.0 ) ) # Both estimates are > 0. Solution is therefore below P₋, i.e. in the range [ 0, P₋ ]
        P₊ = Pₗ
        P₋ = 1e-10
    elseif ( ( f₋ <= 0.0 ) & ( f₊ >= 0.0 ) ) # f₋ less than 0 and f₊ greater than 0, so zero is between them, i.e. in the range [ P₋, P₊ ]
        P₊ = P₊
        P₋ = P₋
    elseif ( ( f₋ < 0.0 ) & ( f₊ <= 0.0 ) ) # Both estimates are less than zero, so solution must be greater than P₊, i.e. in the range [ P₊, ∞ ]
        P₋ = P₊
        P₊ = 1e10
    end

    # Now we iterate to find the solution
    # For this, a simple bisection method will be used until the upper and lower pressure estimates are within some tolerance
    ϵ = 100.0 # Initial error
    ϵₘ = 1e-10 # Maximum error
    i = 0 # Iteration number
    iₘ = 100 # Maximum number of iterations
    while ( ϵ > ϵₘ )
        # Check that we haven't exceeded our maximum iterations
        if ( i >= iₘ )
            @warn "Failed to converge after $i iterations in ToroSolve (P₋: $(P₋) P₊: $(P₊) ϵ: $(ϵ))"
            break
        else
            i=i+1
        end

        # Evaluate the function f with our current pressure estimates
        #f₋ = toro_f( P₋, γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ )
        #f₊ = toro_f( P₊, γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ )
        
        # Center pressure is an average of the left and right values
        P₀ = ( P₋ + P₊ ) / 2.0
        f₀ = toro_f( P₀, γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ )

        # Now need to decide how to move our bounds
        if ( f₀ > 0.0 ) # If f₀ > 0, the solution must lie to the left of P₀, so we reduce P₊ to P₀ and try again
            P₊ = P₀
        else # If f₀ <= 0, the solution must lie to the right of P₀, so we increase P₋ to P₀ and try again
            P₋ = P₀
        end

        ϵ = abs(f₀) # Our error will just be the value of f₀ since that's the root we want to find
    end

    # Iteration is complete
    # P₀ is the average of our left and right estimates
    return (P₊ + P₋)/2.0
end

"""
    CenterVelocity( P₀::T, γₗ::T, Pₗ::T, ρₗ::T, uₗ::T, γᵣ::T, Pᵣ::T, ρᵣ::T, uᵣ::T ) where { T <: AbstractFloat }

Finds the value of the center (contact surface) velocity using Eq. 4.9 along with the value of P₀ found from CenterPressure

# Arguments
    P₀: The pressure in the center region
    γ:  The ratio of specific heats of the gas
    P:  The pressure of the gas
    ρ:  The density of the gas
    u:  The velocity of the gas

    ⋅ Subscripts indicate whether the property is associated with the gas to the left (𝒻ₗ) or right (𝒻ᵣ) of the interface

# Returns
    A value of type T representing the value of P₀ that satisfies Eq. 4.5
"""
function CenterVelocity( P₀::T, γₗ::T, Pₗ::T, ρₗ::T, uₗ::T, γᵣ::T, Pᵣ::T, ρᵣ::T, uᵣ::T ) where { T <: AbstractFloat }
    return 0.5 * ( uₗ + uᵣ ) + 0.5 * ( toro_fᵣ( P₀, γᵣ, Pᵣ, ρᵣ ) - toro_fₗ( P₀, γₗ, Pₗ, ρₗ ) )
end

"""
    CenterVelocity( γₗ::T, Pₗ::T, ρₗ::T, uₗ::T, γᵣ::T, Pᵣ::T, ρᵣ::T, uᵣ::T ) where { T <: AbstractFloat }

Finds the value of the center (contact surface) velocity using Eq. 4.9 along with the value of P₀ found from CenterPressure

# Arguments
    γ:  The ratio of specific heats of the gas
    P:  The pressure of the gas
    ρ:  The density of the gas
    u:  The velocity of the gas

    ⋅ Subscripts indicate whether the property is associated with the gas to the left (𝒻ₗ) or right (𝒻ᵣ) of the interface

# Returns
    A value of type T representing the value of P₀ that satisfies Eq. 4.5

# Notes
    This is a convenience function to compute the center pressure, P₀, as part of the function call rather than having to compute it as a separate step. 
    Generally, it is preferable to call CenterPressure() first and then call the 9 argument version of this function with the result.
"""
function CenterVelocity( γₗ::T, Pₗ::T, ρₗ::T, uₗ::T, γᵣ::T, Pᵣ::T, ρᵣ::T, uᵣ::T ) where { T <: AbstractFloat }
    P₀ = CenterPressure( γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ )
    return CenterVelocity( P₀, γₗ, Pₗ, ρₗ, uₗ, γᵣ, Pᵣ, ρᵣ, uᵣ )
end

"""
    EdgeDensity( P₀::T, γ::T, P::T, ρ::T, u::T ) where { T <: AbstractFloat }

Finds the value of the density outside the center region

# Arguments
    P₀: The pressure in the center region
    γ : The ratio of specific heats of the gas in the edge region
    P : The initial pressure of the gas in the edge region
    ρ : The initial density of the gas in the edge region
    u : The initial velocity of the gas in the edge region

# Returns
    A value of type T representing the density in the edge region. The specific equation used depends on the ratio P₀/P. See Equations section for details.

# Equations
    The exact equation that is evaluated in this function depends on the ratio P₀/P, which determines whether the emitted wave is a shock wave or a rarefaction wave.
    As such, there are two possible cases:

    ⋅ P₀/P > 1 → Shock Wave. The equation solved in this case is
        ρₙ = ρₑ * ρₙ/ρₑ
    where
        ρₙ is the new edge density
        ρₑ is the initial edge density before the passage of the shock
    Note: ρₙ/ρₑ is calculated using the DensityJump function for a moving shock wave from the ShockWaves package.
    It is a function of
        ρₙ/ρₑ = f( γ, P₀/P )
    where
        γ  is the ratio of specific heats for the edge gas
        P₀ is the pressure in the center region, such as calculated using CenterPressure()
        P  is the pressure in the edge region, before the passage of the shock wave
    For more information, see the documentation for DensityJump in the ShockWaves package.

    ⋅ P₀/P ≤ 1 → Rarefaction wave. The equation solved in this case is Eq. 4.23:
        ρₙ = ρₑ * ( P₀/P )^( 1/γ )
    where
        ρₙ is the new edge density
        ρₑ is the initial edge density, before the passage of the rarefaction wave
        P₀ is the pressure in the center region, such as calculated using CenterPressure()
        P  is the pressure in the edge region, before the passage of the rarefaction wave

# Note
    The argument u is unused, but is kept in the function signature to be consistent with other functions in this package.
"""
function EdgeDensity( P₀::T, γ::T, P::T, ρ::T, u::T ) where { T <: AbstractFloat }
    if ( P₀ > P ) # Shock wave
        return ρ * DensityRatio( MovingShockWave(), γ, P₀/P )
    else # Rarefaction wave
        return ρ * ( P₀ / P )^( 1 / γ )
    end
end

export CenterPressure, CenterVelocity, EdgeDensity
