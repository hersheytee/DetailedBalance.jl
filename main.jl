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
    ca_flux = cumtrapz(E, pflux)

    ax4 = Axis(fig[1,4])
    lines!(ax4, E/qe, ca_flux)

    # get integrated photon flux from bandgap to infinity
    a_flux = ca_flux[end] .- ca_flux

    ax4 = Axis(fig[1,5])
    lines!(ax4, E/qe, a_flux)

    return a_flux, fig

end

# cumulative trapezoidal integration function
function cumtrapz(x,y)

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

val_check, figure = detailed_balance("am0.csv", 0.1:0.001:3)

figure