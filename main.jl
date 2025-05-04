using Pkg, CSV, DataFrames, GLMakie, Trapz

function detailed_balance(spectrum_file, Eg)

    # load the spectrum
    spectrum = CSV.read(spectrum_file, DataFrame)

    # wavelengths in m
    wl = spectrum[:, 1]/1e6
    # irradiance in W m^-3
    irad = spectrum[:, 2]*1e6

    # plot wavelength vs irradiance to check its ok
    fig = Figure()
    ax1 = Axis(fig[1, 1], limits = (0, 4e-6, nothing, nothing))
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
    ax2 = Axis(fig[1, 2], limits = (0, 4e-6, nothing, nothing))
    lines!(ax2, wl, T_ev)

    # get photon energy as fn of wavelength
    E = h*c./wl


    # Absorbed photon flux

    # Photon flux in sunlight that has energy above bandgap
    pflux = (2*pi/(h^3*c^2)) .* E.^2 ./ (exp.(E./(k.*T_ev)) .-1)

    # flip pflux and photon energy array from left to right so they can be integrated correctly
    pflux = reverse(pflux)
    E = reverse(E)

    ax3 = Axis(fig[1,3])
    lines!(ax3, E/qe, pflux)

    # integrate from bandgap to infinity 
    # TODO replace with IRZI method

    # energy in eV
    E_eV = E/qe
 
    # cumulative integrated absorbed photon flux
    ca_flux = cumtrap_int(E, pflux)

    ax4 = Axis(fig[1,4])
    lines!(ax4, E/qe, ca_flux)

    # get integrated photon flux from bandgap to infinity
    a_flux = ca_flux[end] .- ca_flux

    ax4 = Axis(fig[1,5])
    lines!(ax4, E/qe, a_flux)

    
    # Emitted photon flux

    # Create the emitted flux array
    e_flux = Inf*ones(length(E))

    for i in eachindex(E)

        # all voltages up to the current photon energy
        V = 0:0.01:E[i]

        # iterate through voltages to get emitted flux values for each bandgap and cell voltages
    end
    
    # create tuple holding outputs
    outputs = (a_flux, E)

    return outputs, fig

end

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

# cumulative trapezoidal integration function
function cumtrap_int(x,y)

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

function photon_emission_flux(Eg, V, T)
    # returns photon emission flux for specified external bias, cell bandgap
    # and temperature

    # define constants - TODO, use constants from constants.jl
    k = 1.380649e-23 # Boltzmann's constant
    h = 6.62607015e-34 # Planck's constant
    c = 299792458 # Speed of light
    qe = 1.602176634e-19 # Electron charge

    # bandgap in J
    Eg = Eg*qe

    # energies we want to integrate over
    Ee = Eg:0.0001:5*qe

    # photon emission flux function
    e_flux_func = ( 2*pi/(h^3*c^2) ) * Ee.^2 ./ ( exp( (Ee-qe*V)/(k*T) ) - 1);

    # perform trapezoidal integration
    e_flux = trap_int(Ee, e_flux_func)

    return e_flux

end

val_check, figure = detailed_balance("am0.csv", 0.1:0.001:3)

figure