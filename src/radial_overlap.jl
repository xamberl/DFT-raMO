#==
This supporting file can analytically calculate the radial overlap integrals.
The functions are cross-checked with Wolfram Mathematica.
Wolfram Mathematica's solutions are simplified in an easier to code/read format, so those solutions are hardcoded in,
in case you may be wondering why there are discrepancies between the output here and in calculate_overlap().
==#

using SymPy
using SpecialFunctions
"""
    radial_overlap(n_quant::Int, l_quant::Int) -> a::Sym

Returns a SymPy object and prints the radial overlap solution at corresponding quantum numbers n, l.
"""
function radial_overlap(n_quant::Int, l_quant::Int)
    # Radial component of STOs 
    @vars n ζ r
    rad_STO = (2*ζ)^n*sqrt((2*ζ)/factorial(2*n))*r^(n-1)*exp(-ζ*r)
    
    # Spherical bessel functions
    @vars x
    sph_bes0 = sin(x)/x
    sph_bes1 = sin(x)/x^2 - cos(x)/x
    sph_bes2 = (3/x^2 - 1)(sin(x)/x) - (3*cos(x))/x^2
    
    # Return the integral
    @vars k
    l_quant == 0 ? a = integrate(rad_STO(n=>n_quant)*sph_bes0(x=>k*r)*r^2,(r,0,oo)) : nothing
    l_quant == 1 ? a = integrate(rad_STO(n=>n_quant)*sph_bes1(x=>k*r)*r^2,(r,0,oo)) : nothing
    l_quant == 2 ? a = integrate(rad_STO(n=>n_quant)*sph_bes2(x=>k*r)*r^2,(r,0,oo)) : nothing
    
    return a
end

"""
    eval_overlap(a::Sym, z::Real, k_g::Real) -> Float64

Evaluates and returns the overlap at specified ζ and |k+G|.
"""
function eval_overlap(a::Sym, z::Real, k_g::Real)
    @vars ζ k
    return a(k=>k_g, ζ=>z).evalf()
end

"""
    test_timing()

Timing test to demonstrate practical reasons to hardcode the analytical solutions.
Example used is the H 1s orbital with |k+G| = 9.92118191847175e-16
"""
function test_timing()
    z = 6.977549
    k_g = 9.92118191847175e-16
    @info "Numerical integration:"
    @time eval_overlap(radial_overlap(1,0),z,k_g)
    @info "Analytical solution:"
    @time (4*z^(5/2))/(k_g^2+z^2)^2;
end
