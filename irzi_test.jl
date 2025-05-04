using GLMakie, Pkg, CSV, DataFrames, GLMakie, Trapz, Interpolations



# function detailed_balance(spectrum_file, Eg)

#     # load the spectrum
#     spectrum = CSV.read(spectrum_file, DataFrame)

#     # wavelengths in m
#     wl = spectrum[:, 1]/1e6
#     # irradiance in W m^-3
#     irad = spectrum[:, 2]*1e6

#     # plot wavelength vs irradiance to check its ok
#     fig = Figure()
#     ax1 = Axis(fig[1, 1], limits = (0, 4e-6, nothing, nothing))
#     lines!(ax1, wl, irad)

#     # Get total power
#     p_max = trapz(wl, irad)  

#     # define constants
#     k = 1.380649e-23 # Boltzmann's constant
#     h = 6.62607015e-34 # Planck's constant
#     c = 299792458 # Speed of light
#     qe = 1.602176634e-19 # Electron charge

#     # Planckian effective temperature
#     T_ev = (h*c)./ (wl*k.*log.((2*h*pi*c^2)./(irad.*wl.^5) .+ 1))

#     # plot temp to make sure it's ok
#     ax2 = Axis(fig[1, 2], limits = (0, 4e-6, nothing, nothing))
#     lines!(ax2, wl, T_ev)



#     # Absorbed photon flux

#     # Photon flux in sunlight that has energy above bandgap
#     pflux = (2*pi/(h^3*c^2)) .* E.^2 ./ (exp.(E./(k.*T_ev)) .-1)

#     # flip pflux and photon energy array from left to right so they can be integrated correctly
#     pflux = reverse(pflux)
#     E = reverse(E)

#     ax3 = Axis(fig[1,3])
#     lines!(ax3, E/qe, pflux)

#     # integrate from bandgap to infinity 
#     # TODO replace with IRZI method

#     # energy in eV
#     E_eV = E/qe
 
#     # cumulative integrated absorbed photon flux
#     ca_flux = cumtrap_int(E, pflux)

#     ax4 = Axis(fig[1,4])
#     lines!(ax4, E/qe, ca_flux)

#     # get integrated photon flux from bandgap to infinity
#     a_flux = ca_flux[end] .- ca_flux

#     ax4 = Axis(fig[1,5])
#     lines!(ax4, E/qe, a_flux)

    
#     # Emitted photon flux

#     # Create the emitted flux array
#     e_flux = Inf*ones(length(E))

#     for i in eachindex(E)

#         # all voltages up to the current photon energy
#         V = 0:0.01:E[i]

#         # iterate through voltages to get emitted flux values for each bandgap and cell voltages
#     end
    
#     # create tuple holding outputs
#     outputs = (a_flux, E)

#     return outputs, fig

# end



# val_check, figure = detailed_balance("am0.csv", 0.1:0.001:3)

function integrand(x, p)
    # function computes the value of the Riemann zeta integrand for a given value of x and order p
    val = 1/(x^(p+1)*(exp(1/x)-1))

    return val
end

x = exp10.(range(-3,7,5000))


# trapezoidal integration function
function trap_int(x,y)

    val = 0

    for i in eachindex(x)        
        if i != 1
            val = val + (x[i] - x[i-1])*(y[i] + y[i-1])/2
        end
    end

    return val
end


function cumtrap_int(x,y)
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

fig = Figure()
ax1 = Axis(fig[1,1], xscale = log10)

function rzi_lut(; x1=1.4e-3, x2=1e7, num=5000)
    # function returns look-up tables for values of the Incomplete Riemann-Zeta integral for orders p=1,2,3,4
    # x1 is the lower bound of the integration domain
    # x2 is the upper bound of the integration domain
    # n is the number of logarithmically-spaced points in the integration domain

    p_range = 1:4 # range of orders of the integral

    luts = [] # array to hold the look-up table values

    # array of x values in domain
    x = exp10.(range(log10(x1),log10(x2),num))

    for p in p_range

        # integrand values
        integrand = 1 ./(x.^(p+1).*(exp.(1 ./x).-1))

        # cumulative trapezoidal integration of all the values of the integrand
        val = cumtrap_int(x, integrand)
        
        # return lookup table as interpolations function 

        # TODO: see if this is the fastest way
        push!(luts, linear_interpolation(x, val))

        # lines!(ax1, x, val)

    end

    return luts
    
end



# 

function lower_bound(x, p)
    # function returns the lower bound approximation of the Incomplete Riemann-Zeta integral
    # x is the lower bound of the integration domain
    # p is the order of the integral

    if p == 1
        val = e^(-1/x)
    elseif p == 2
        val = (x+1)/(x*exp(1/x))
    elseif p == 3
        val = (2*x^2+2*x+1)/(x^2*exp(1/x))
    elseif p == 4
        val = (6*x^3+6*x^2+3*x+1)/(x^3*exp(1/x))
    end
    return val
end

lower_bound(x, p) = (x+1)/(x*exp(1/x))
upper_bound(x, p) = -x^(1-p)/(1-p)

# println(lower_bound(x1))
# println(upper_bound(x2))

# define constants
k = 1.380649e-23 # Boltzmann's constant
h = 6.62607015e-34 # Planck's constant
c = 299792458 # Speed of light
qe = 1.602176634e-19 # Electron charge

T = 300
e_b = 3*qe # upper bound of the integration domain
e_a = 0.1*qe # lower bound of the integration domain
mu = 0

# construct lookup tables for the Incomplete Riemann-Zeta integral for orders p=1,2,3,4
luts = rzi_lut()

# println(lut[end,2]-pi^4/15)

function rzi(x1, x2, p, luts)
    # function returns the value of the Incomplete Riemann-Zeta integral
    # x1 is the lower bound of the integration domain
    # x2 is the upper bound of the integration domain
    # p is the order of the integral
    # lut is the look-up tables for the Incomplete Riemann-Zeta integral created by the rzi_lut function

    # select the appropriate look-up table for the order of the integral
    lut = luts[p]

    # get integration domain limit
    x_low = lut[1]
    x_high = lut[end,1]

    # create variable to hold the value of the integral
    integral_val = 0

    # check if x1 is less than the lower bound of the integration domain
    if x1 < x_low

        # use the lower bound approximation and add to integral value

        # check if x1 is 0
        if x1 == 0
            # don't include x1 as 0 val will lead to NaN
            integral_val += lower_bound(x_low, p)
        else
            integral_val += (lower_bound(x_low, p) - lower_bound(x1, p))
        end

        # set x1 to the lower bound of the integration domain
        x1 = x_low

    end

    if x2 > x_high
        # use the upper bound approximation and add to integral value
        integral_val += upper_bound(x2, p) - upper_bound(x_high, p)

        # set x2 to the upper bound of the integration domain
        x2 = x_high

    end

    # use lut to get the value of the integral between x1 and x2 and add to integral value
    integral_val += (lut(x2) - lut(x1))

    println(x_low)

    return integral_val
    
end

lut = luts[1]

# define x1 and x2 for the Bose-Einstein integral
x1 = k*T/((e_b - mu))
println(x1)
# x1 = 0
x2 = k*T/((e_a - mu))

function bei(; x1, x2, mu, T, luts)
    # function returns the value of the Bose-Einstein distribution integral
    # x1 is kT/(qe*(e_b - mu))
    # x2 is kT/(qe*(e_a - mu))
    # e_b is the upper bound of the integration
    # mu is the chemical potential
    # T is the temperature in K
    # lut is the look-up table for the Incomplete Riemann-Zeta integral created by the rzi_lut function

    # define constants
    k = 1.380649e-23 # Boltzmann's constant
    h = 6.62607015e-34 # Planck's constant
    c = 299792458 # Speed of light
    qe = 1.602176634e-19 # Electron charge

    # define prefactor
    prefactor = (2*pi*k*T)/(h^3*c^2)

    # define the values for the Bose-Einstein integral from the Riemann zeta integral for the (p-1)th order
    rzi0 = mu^3 * rzi(x1, x2, 1, luts)
    rzi1 = 3*k*T*mu^2 * rzi(x1, x2, 2, luts)
    rzi2 = 3*(k*T)^2*mu * rzi(x1, x2, 3, luts)
    rzi3 = (k*T)^3 * rzi(x1, x2, 4, luts)

    # calculate Bose-Einstein integral from Riemann zeta integral
    bei_val = prefactor * (rzi0 + rzi1 + rzi2 + rzi3)

    return bei_val

end

bei_val = bei(x1=x1, x2=x2, mu=0, T=T, luts=luts)

display(fig)

