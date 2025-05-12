using Pkg, CSV, DataFrames, GLMakie, Trapz, Interpolations, Roots, ProgressBars

function detailed_balance(; spectrum_file, T, E=-1, num_voltages=100)

    # load the spectrum
    spectrum = CSV.read(spectrum_file, DataFrame)

    # wavelengths in m
    wl = spectrum[:, 1]/1e6
    # irradiance in W m^-3
    irad = spectrum[:, 2]*1e6

    # Get total power to check
    p_total = trapz(wl, irad)

    # save spectrum parameters as dictionary
    spectrum_params = Dict("wl" => wl, "irad" => irad, "p_total" => p_total)

    # define constants
    k = 1.380649e-23 # Boltzmann's constant
    h = 6.62607015e-34 # Planck's constant
    c = 299792458 # Speed of light
    qe = 1.602176634e-19 # Electron charge

    # Planckian effective temperature
    T_ev = (h*c)./ (wl*k.*log.((2*h*pi*c^2)./(irad.*wl.^5) .+ 1))

    # get photon energy as fn of wavelength from data
    E_data = h*c./wl

    # flip photon energy array from left to right so it's increasing
    E_data = reverse(E_data)
    # correct the effective temperature array
    T_ev = reverse(T_ev)

    ### Absorbed photon flux ###

    # Photon flux in sunlight that has energy above bandgap
    pflux = (2*pi/(h^3*c^2)) .* E_data.^2 ./ (exp.(E_data./(k.*T_ev)) .-1)

    # integrate from bandgap to infinity 

    # cumulative integrated absorbed photon flux
    ca_flux = cumtrapz_int(E_data, pflux)
 
    # get integrated photon flux from bandgap to infinity
    a_flux_data = ca_flux[end] .- ca_flux

    # interpolate absorbed flux for bandgaps defined

    a_flux_data_func = linear_interpolation(E_data, a_flux_data)

    # if E is not defined, use the data from the spectrum
    if E == -1
        E = E_data
    else
        # convert E to Joules
        E = E*qe
    end

    ### Emitted photon flux using Incomplete Riemann-Zeta Integral ###

    # generate lookup tables
    luts = rzi_lut(lut_lim=[1.4e-3, 1e7])

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
            _, _, P_max = IV(mu=V_P_max*qe, E=E[i], T=T, luts=luts, a_flux_data_func=a_flux_data_func)
            
        end

        push!(V_oc_range, V_oc)        
        push!(P_max_range, P_max)
    
    end

    # store in dictionary
    IV_params = Dict("V_oc" => V_oc_range, "J_sc" => J_sc_range, "P_max" => P_max_range, "eff" => P_max_range/p_total)

    # create tuple holding outputs
    outputs = (E, spectrum_params, IV_curves, IV_params)

    return outputs

end