using Pkg, CSV, DataFrames, GLMakie, Trapz, Interpolations, Roots, ProgressBars


include("irzi_functions.jl")

# define constants
k = 1.380649e-23 # Boltzmann's constant
h = 6.62607015e-34 # Planck's constant
c = 299792458 # Speed of light
qe = 1.602176634e-19 # Electron charge



# function detailed_balance_fluxes(; spectrum_file, T=300)

#     # define constants
#     k = 1.380649e-23 # Boltzmann's constant
#     h = 6.62607015e-34 # Planck's constant
#     c = 299792458 # Speed of light
#     qe = 1.602176634e-19 # Electron charge
    
#     # load the spectrum
#     spectrum = CSV.read(spectrum_file, DataFrame)

#     # wavelengths in m
#     wl = spectrum[:, 1]/1e6
#     # irradiance in W m^-3
#     irad = spectrum[:, 2]*1e6

#     # Calculate total power in spectrum per unit area
#     p_max = trapz(wl, irad)

#     # Planckian effective temperature
#     T_ev = (h*c)./ (wl*k.*log.((2*h*pi*c^2)./(irad.*wl.^5) .+ 1))

#     # get photon energy as fn of wavelength
#     E = h*c./wl

#     # flip photon energy array from left to right so it's increasing
#     E = reverse(E)
#     # correct the effective temperature array
#     T_ev = reverse(T_ev)

#     ### Absorbed photon flux ###

#     # Photon flux in sunlight that has energy above bandgap
#     pflux = (2*pi/(h^3*c^2)) .* E.^2 ./ (exp.(E./(k.*T_ev)) .-1)

#     # integrate from bandgap to infinity 

#     # cumulative integrated absorbed photon flux
#     ca_flux = cumtrapz_int(E, pflux)

#     # get integrated photon flux from bandgap to infinity
#     a_flux = ca_flux[end] .- ca_flux

#     ### Emitted photon flux using Incomplete Riemann-Zeta integral ###

#     # generate lookup tables
#     luts = rzi_lut(lut_lim=[1.4e-3, 1e7])

#     # number of voltages to iterate over
#     num_voltages = 100

#     # Create the emitted flux array
#     e_flux = Inf*ones(length(E), num_voltages)

#     # create an array to store voltages
#     V_range = []

#     for i in eachindex(E)

#         # 1000 steps for all voltages up to the current photon energy - 0.01 to avoid 0 error
#         V = range(0, E[i]-0.01*qe, length=num_voltages)

#         # store voltages in the array
#         V_range = push!(V_range, V)

#         # iterate through voltages to get emitted flux values for each bandgap and cell voltages
#         for j in eachindex(V)

#             # define x1 and x2 for the Bose-Einstein integral
#             x1 = 0 # x1 = k*T/((e_b - V)) and it is assumed e_b is infinity
#             x2 = k*T/((E[i] - V[j])) # e_a is bandgap

#             # emitted photon flux
#             e_flux[i,j] = bei(x1=x1, x2=x2, mu=V[j], T=T, luts=luts, order=2)

#         end

#     end


#     # create tuple holding outputs
#     outputs = (a_flux=a_flux, e_flux=e_flux, E=E, V_range=V_range, )

#     return outputs

# end


function detailed_balance(; spectrum_file, T, E=-1)

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

    # Get total power to check
    p_total = trapz(wl, irad)  

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
 
    # get photon energy as fn of wavelength from data
    E_data = h*c./wl

    # flip photon energy array from left to right so it's increasing
    E_data = reverse(E_data)
    # correct the effective temperature array
    T_ev = reverse(T_ev)

    ### Absorbed photon flux ###

    # Photon flux in sunlight that has energy above bandgap
    pflux = (2*pi/(h^3*c^2)) .* E_data.^2 ./ (exp.(E_data./(k.*T_ev)) .-1)

    ax3 = Axis(fig[1,3], xlabel = "Energy (eV)", ylabel = "Photon Flux (m^-2 s^-1)")
    lines!(ax3, E_data/qe, pflux)

    # integrate from bandgap to infinity 

    # cumulative integrated absorbed photon flux
    ca_flux = cumtrapz_int(E_data, pflux)

    ax4 = Axis(fig[1,4], xlabel = "Energy (eV)", ylabel = "Cumulative Absorbed Photon Flux (m^-2 s^-1)")
    lines!(ax4, E_data/qe, ca_flux)
 
    # get integrated photon flux from bandgap to infinity
    a_flux_data = ca_flux[end] .- ca_flux

    ax5 = Axis(fig[1,5], xlabel = "Bandgap Energy (eV)", ylabel = "Integrated Absorbed Photon Flux (m^-2 s^-1)")
    lines!(ax5, E_data/qe, a_flux_data)

    # interpolate absorbed flux for bandgaps defined

    a_flux_data_func = linear_interpolation(E_data, a_flux_data)

    # if E is not defined, use the data from the spectrum
    if E == -1
        E = E_data
    end

    ### Emitted photon flux using Incomplete Riemann-Zeta Integral ###

    # generate lookup tables
    luts = rzi_lut(lut_lim=[1.4e-3, 1e7])

    # number of voltages to iterate over
    num_voltages = 100

    # Create the emitted flux array
    e_flux = Inf*ones(length(E), num_voltages)


    # create an array to store voltages, currents, and power
    V_range = []
    J_range = []
    P_range = []

    ### generate IV curve data ###

    for i in eachindex(E)

        # 1000 steps for all chemical potentials up to the current photon energy - 0.001 eV to avoid 0 error
        mu = range(0, E[i]-0.001*qe, length=num_voltages)

        # set up arrays to store IV curves
        J = zeros(num_voltages)
        V = zeros(num_voltages)
        P = zeros(num_voltages)

        # iterate through voltages to get IV curve
        for j in eachindex(mu)

            J[j], V[j], P[j] = IV(mu=mu[j], E=E[i], T=T, luts=luts, a_flux_data_func=a_flux_data_func)

        end

        # store voltage, current, and power ranges
        V_range = push!(V_range, V)
        J_range = push!(J_range, J)
        P_range = push!(P_range, P)

    end

    # store in dictionary
    IV_curves = Dict("J" => J_range, "V" => V_range, "P" => P_range)
    
    ### generate V_oc, J_sc, P_max, and eff data ###

    # set up array to store parameters for each bandgap
    V_oc_range = []
    J_sc_range = []
    P_max_range = []

    for i in tqdm(eachindex(E))

        # find and store the short circuit current
        J_sc, _ = IV(mu=0, E=E[i], T=T, luts=luts, a_flux_data_func=a_flux_data_func)
        push!(J_sc_range, J_sc)

        # check that J_sc > 0
        if J_sc < 0
            V_oc = NaN
            P_max = NaN
        else
            # get the open cicruit voltage using the roots function and store
            V_oc = find_zero(mu -> IV(mu=mu, E=E[i], T=T, luts=luts, a_flux_data_func=a_flux_data_func)[1], (0, E[i]-0.001*qe))/qe

            # find the voltage at maximum power point
            V_P_max = find_zero(mu -> d_power(mu=mu, E=E[i], T=T, luts=luts, a_flux_data_func=a_flux_data_func), (V_oc*qe/3, V_oc*qe))/qe

            # get the maximum power point
            J_P_max, _, P_max = IV(mu=V_P_max*qe, E=E[i], T=T, luts=luts, a_flux_data_func=a_flux_data_func)
            
        end

        push!(V_oc_range, V_oc)        
        push!(P_max_range, P_max)
    
    end

    # store in dictionary
    IV_parameters = Dict("V_oc" => V_oc_range, "J_sc" => J_sc_range, "P_max" => P_max_range, "eff" => P_max_range/p_total)

    # create tuple holding outputs
    outputs = (E, IV_curves, IV_parameters, p_total)

    return outputs, fig

end

outputs, figure = detailed_balance(spectrum_file = "am0.csv", T=298.15)


E = outputs[1]/qe

IV_curves = outputs[2]
IV_parameters = outputs[3]
p_total = outputs[4]
;

fig = Figure()
ax1 = Axis(fig[1, 1], limits=(0, nothing, 0, nothing))
lines!(ax1, E, IV_parameters["eff"])

;
fig
#integral_val = rzi(0, 0.05, 1, luts)
