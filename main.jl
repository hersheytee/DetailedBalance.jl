using Pkg, CSV, DataFrames, GLMakie, Trapz, Interpolations


include("irzi_functions.jl")

# define constants
k = 1.380649e-23 # Boltzmann's constant
h = 6.62607015e-34 # Planck's constant
c = 299792458 # Speed of light
qe = 1.602176634e-19 # Electron charge

function detailed_balance(; spectrum_file)

    # load the spectrum
    spectrum = CSV.read(spectrum_file, DataFrame)

    # wavelengths in m
    wl = spectrum[:, 1]/1e6
    # irradiance in W m^-3
    irad = spectrum[:, 2]*1e6

    # plot wavelength vs irradiance to check its ok
    fig = Figure()
    ax1 = Axis(fig[1, 1], limits = (0, 4e-6, nothing, nothing), xlabel = "Wavelength (m)", ylabel = "Irradiance (W/m^2)")
    lines!(ax1, wl, irad)

    # Get total power
    p_max = trapz(wl, irad)  

    # define constants
    k = 1.380649e-23 # Boltzmann's constant
    h = 6.62607015e-34 # Planck's constant
    c = 299792458 # Speed of light
    qe = 1.602176634e-19 # Electron charge

    # Planckian effective temperature
    T_ev = (h*c)./ (wl*k.*log.((2*h*pi*c^2)./(irad.*wl.^5) .+ 1))

    # plot temp to make sure it's ok
    ax2 = Axis(fig[1, 2], limits = (0, 4e-6, nothing, nothing), xlabel = "Wavelength (m)", ylabel = "PLanckian Effective Temperature (K)")
    lines!(ax2, wl, T_ev,)

    # get photon energy as fn of wavelength
    E = h*c./wl

    # Absorbed photon flux

    # Photon flux in sunlight that has energy above bandgap
    pflux = (2*pi/(h^3*c^2)) .* E.^2 ./ (exp.(E./(k.*T_ev)) .-1)

    # flip pflux and photon energy array from left to right so they can be integrated correctly
    pflux = reverse(pflux)
    E = reverse(E)

    ax3 = Axis(fig[1,3], xlabel = "Energy (eV)", ylabel = "Photon Flux (m^-2 s^-1)")
    lines!(ax3, E/qe, pflux)

    # integrate from bandgap to infinity 
    # TODO replace with IRZI method
 
    # cumulative integrated absorbed photon flux
    ca_flux = cumtrapz_int(E, pflux)

    ax4 = Axis(fig[1,4], xlabel = "Energy (eV)", ylabel = "Cumulative Absorbed Photon Flux (m^-2 s^-1)")
    lines!(ax4, E/qe, ca_flux)

    # get integrated photon flux from bandgap to infinity
    a_flux = ca_flux[end] .- ca_flux

    ax5 = Axis(fig[1,5], xlabel = "Bandgap Energy (eV)", ylabel = "Integrated Absorbed Photon Flux (m^-2 s^-1)")
    lines!(ax4, E/qe, a_flux)

    # Emitted photon flux

    # generate lookup tables
    luts = rzi_lut(lut_lim=[1.4e-3, 1e7])
    # println(luts)
    
    T = 300 # temperature of solar cell in K

    # number of voltages to iterate over
    num_voltages = 100

    # Create the emitted flux array
    e_flux = Inf*ones(length(E), num_voltages)

    # create an array to store voltages
    V_range = []

    for i in eachindex(E)

        # 1000 steps for all voltages up to the current photon energy - 0.01 to avoid 0 error
        V = range(0, E[i]-0.01*qe, length=num_voltages)

        # store voltages in the array
        V_range = push!(V_range, V)

        # iterate through voltages to get emitted flux values for each bandgap and cell voltages
        for j in eachindex(V)

            # define x1 and x2 for the Bose-Einstein integral
            x1 = 0 # x1 = k*T/((e_b - V)) and it is assumed e_b is infinity
            x2 = k*T/((E[i] - V[j])) # e_a is bandgap

            # emitted photon flux
            try
                e_flux[i,j] = bei(x1=x1, x2=x2, mu=V[j], T=T, luts=luts, order=2)
            catch
                #println("Error at i = $i, j = $j")
                #println("x1 = $x1, x2 = $x2, mu = $(V[j]), T = $T")
                break
            end

        end

    end


    
    # go through every bandgap and get currents
    for i in [1000]

        # get the voltage range
        V = V_range[i]

        ax6 = Axis(fig[1,6], limits = (0, nothing, 0, nothing), 
        xlabel = "Voltage", ylabel = "Current", title = "bandgap = $(E[i]/qe) eV")

        # get current for each voltage
        J = qe*(a_flux[i] .- e_flux[i, :])

        lines!(ax6, V/qe, J)
        
    end

    # create tuple holding outputs
    outputs = (a_flux, e_flux, E)

    return outputs, fig

end

# # trapezoidal integration function
# function trap_int(x,y)

#     val = 0

#     for i in eachindex(x)        
#         if i != 1
#             val = val + (x[i] - x[i-1])*(y[i] + y[i-1])/2
#         end
#     end

#     return val
# end

# # cumulative trapezoidal integration function
# function cumtrap_int(x,y)

#     val = []

#     for i in eachindex(x)

#         if i == 1
#             push!(val, 0)
#         else
#             push!(val, val[i-1] + (x[i] - x[i-1])*(y[i] + y[i-1])/2)
#         end
#     end

#     return Float64.(val)
# end

# function photon_emission_flux(Eg, V, T)
#     # returns photon emission flux for specified external bias, cell bandgap
#     # and temperature

#     # define constants - TODO, use constants from constants.jl
#     k = 1.380649e-23 # Boltzmann's constant
#     h = 6.62607015e-34 # Planck's constant
#     c = 299792458 # Speed of light
#     qe = 1.602176634e-19 # Electron charge

#     # bandgap in J
#     Eg = Eg*qe

#     # energies we want to integrate over
#     Ee = Eg:0.0001:5*qe

#     # photon emission flux function
#     e_flux_func = ( 2*pi/(h^3*c^2) ) * Ee.^2 ./ ( exp( (Ee-qe*V)/(k*T) ) - 1);

#     # perform trapezoidal integration
#     e_flux = trap_int(Ee, e_flux_func)

#     return e_flux

# end

val_check, figure = detailed_balance(spectrum_file = "am0.csv")

a_flux = val_check[1]
e_flux = val_check[2]

E = val_check[3]/qe


#integral_val = rzi(0, 0.05, 1, luts)

figure