using Pkg, CSV, DataFrames, GLMakie, Trapz, Interpolations

function cumtrapz_int(x,y)
    # cumulative trapezoidal integration
    # input: 
    # n x 1 array of Float64 of x values
    # n x 1 array of Float64 of y values
    # output n x 1 array of Float64 of trapezoidal integral of y up to corresponding x value

    val = []

    for i in eachindex(x)

        if i == 1
            push!(val, 0)
        else
            push!(val, val[i-1] + (x[i] - x[i-1])*(y[i] + y[i-1])/2)
        end
    end

    return Float64.(val)
end

function rzi_lut(; lut_lim=[1.4e-3, 1e7], num=5000)
    # function returns look-up tables for values of the Incomplete Riemann-Zeta integral for orders p=1,2,3,4
    # lut_lim is a 2-element array of the lower and upper bounds of the look-up table domain
    # num is the number of logarithmically-spaced points in the look-up table domain

    x1, x2 = lut_lim # lower and upper bounds of the integration domain

    p_range = 1:4 # range of orders of the integral

    luts = [] # array to hold the look-up table values

    # array of x values in domain
    x = exp10.(range(log10(x1),log10(x2),num))

    for p in p_range

        # integrand values
        integrand = 1 ./(x.^(p+1).*(exp.(1 ./x).-1))

        # cumulative trapezoidal integration of all the values of the integrand
        val = cumtrapz_int(x, integrand)
        
        # return lookup table as interpolations function 
        lut = linear_interpolation(x, val, extrapolation_bc=Interpolations.Flat())

        # TODO: see if this is the fastest way
        push!(luts, lut)

        # lines!(ax1, x, val)

    end

    return luts
    
end

function rzi_lower_bound(x, p)
    # function returns the lower bound approximation of the Incomplete Riemann-Zeta integral
    # x is the lower bound of the integration domain
    # p is the order of the integral

    if p == 1
        val = exp(-1/x)
    elseif p == 2
        val = (x+1)/(x*exp(1/x))
    elseif p == 3
        val = (2*x^2+2*x+1)/(x^2*exp(1/x))
    elseif p == 4
        val = (6*x^3+6*x^2+3*x+1)/(x^3*exp(1/x))
    end
    return val
end

# function returns the upper bound approximation of the Incomplete Riemann-Zeta integral
# x is the upper bound of the integration domain
# p is the order of the integral
rzi_upper_bound(x, p) = -x^(1-p)/(1-p)

function rzi(x1, x2, p, luts; lut_lim=[1.4e-3, 1e7])
    # function returns the value of the Incomplete Riemann-Zeta integral
    # x1 is the lower bound of the integration domain
    # x2 is the upper bound of the integration domain
    # lut_lim is a 2-element array of the lower and upper bounds of the look-up table domain
    # p is the order of the integral
    # lut is the look-up tables for the Incomplete Riemann-Zeta integral created by the rzi_lut function

    x_low = lut_lim[1] # lower bound of the integration domain
    x_high = lut_lim[2] # upper bound of the integration domain

    # select the appropriate look-up table for the order of the integral
    lut = luts[p]

    # create variable to hold the value of the integral
    integral_val = 0

    # check if x1 is less than the lower bound of the integration domain
    if x1 < x_low

        # use the lower bound approximation and add to integral value

        # check if x1 is 0
        if x1 == 0
            # don't include x1 as 0 val will lead to NaN
            integral_val += rzi_lower_bound(x_low, p)
        else
            integral_val += (rzi_lower_bound(x_low, p) - rzi_lower_bound(x1, p))
        end

        # set x1 to the lower bound of the integration domain
        x1 = x_low

    end

    if x2 > x_high
        # use the upper bound approximation and add to integral value
        integral_val += rzi_upper_bound(x2, p) - rzi_upper_bound(x_high, p)

        # set x2 to the upper bound of the integration domain
        x2 = x_high

    end

    # use lut to get the value of the integral between x1 and x2 and add to integral value
    integral_val += (lut(x2) - lut(x1))

    return integral_val
    
end

function bei(; x1, x2, mu, T, luts, order)
    # function returns the value of the Bose-Einstein distribution integral
    # x1 is kT/(qe*(e_b - mu))
    # x2 is kT/(qe*(e_a - mu))
    # e_b is the upper bound of the integration
    # mu is the chemical potential
    # T is the temperature in K
    # luts is a 4 x 1 array with the look-up tables for the Incomplete Riemann-Zeta integral 
    # created by the rzi_lut function
    # order is the order of the integral, order = 2 for particle flux, order = 3 for energy flux

    # define constants
    k = 1.380649e-23 # Boltzmann's constant
    h = 6.62607015e-34 # Planck's constant
    c = 299792458 # Speed of light
    qe = 1.602176634e-19 # Electron charge

    # define prefactor
    prefactor = (2*pi*k*T)/(h^3*c^2)

    # define the values for the Bose-Einstein integral from the Riemann zeta integrals for different orders
    if order == 2
        # particle flux
        rzi0 = mu^2 * rzi(x1, x2, 1, luts)
        rzi1 = 2*k*T*mu * rzi(x1, x2, 2, luts)
        rzi2 = (k*T)^2 * rzi(x1, x2, 3, luts)
        rzi3 = 0
    elseif order == 3
        # energy flux
        rzi0 = mu^3 * rzi(x1, x2, 1, luts)
        rzi1 = 3*k*T*mu^2 * rzi(x1, x2, 2, luts)
        rzi2 = 3*(k*T)^2*mu * rzi(x1, x2, 3, luts)
        rzi3 = (k*T)^3 * rzi(x1, x2, 4, luts)
    end
    
    # calculate Bose-Einstein integral from Riemann zeta integral
    bei_val = prefactor * (rzi0 + rzi1 + rzi2 + rzi3)

    return bei_val

end

function IV(; mu, E, T, luts, a_flux_data_func)

    qe = 1.602176634e-19 # Electron charge

    # get absorbed photon flux
    a_flux = a_flux_data_func(E)

    # define x1 and x2 for the Bose-Einstein integral
    x1 = 0 # x1 = k*T/((e_b - mu)) and it is assumed e_b is infinity
    x2 = k*T/((E - mu)) # e_a is bandgap

    # emitted photon flux
    e_flux = bei(x1=x1, x2=x2, mu=mu, T=T, luts=luts, order=2)

    # get voltage
    V = mu/qe

    # get current
    J = qe*(a_flux - e_flux)

    # get power
    P = J*V

    # return current and voltage and power
    return J, V, P

end

function d_power(; mu, E, T, luts, a_flux_data_func)
    # function returns the power derivative with respect to voltage using the central difference method
    # function is used for finding the maximum power point
    # mu is the chemical potential
    # E is the photon energy
    # T is the temperature in K
    # luts is the look-up tables for the Incomplete Riemann-Zeta integral created by the rzi_lut function
    # a_flux_data_func is a function that returns the absorbed photon flux

    d_tol = 1e-3*mu # tolerance for the central difference method

    # get the power at mu + d_tol
    P_plus = IV(mu=mu+d_tol, E=E, T=T, luts=luts, a_flux_data_func=a_flux_data_func)[3]
    # get the power at mu - d_tol
    P_minus = IV(mu=mu-d_tol, E=E, T=T, luts=luts, a_flux_data_func=a_flux_data_func)[3]

    d_power = (P_plus - P_minus)/(2*d_tol)

    return d_power

end
