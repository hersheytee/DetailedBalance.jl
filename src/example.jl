# import the packages we need
using DetailedBalance, CairoMakie

# define the bandgaps we want to calculate the detailed balance for (in eV)
bandgaps = 0.5:0.1:2.5

# define the solar cell temperature in K
T = 300

# define the spectrum file (see the included spectrum files in the src folder for the format if you want to create your own)
spectrum_file = "src/am0.csv"

# call the detailed balance function
outputs = detailed_balance(spectrum_file=spectrum_file, T=T, E=bandgaps)

# decompose the outputs
E = outputs[1] # photon energy
spectrum_params = outputs[2] # spectrum parameters
IV_curves = outputs[3] # IV curves
IV_parameters = outputs[4] # IV parameters

# do some plotting
fig = Figure()
ax1 = Axis(fig[1, 1:2], title="T = $T K", xlabel="Bandgap (eV)", ylabel="Solar Cell Efficiency")
scatter!(ax1, bandgaps, IV_parameters["eff"])
ax2 = Axis(fig[2, 1], title="T = $T K", xlabel="Bandgap (eV)", ylabel="Open-Circuit Voltage (V)")
scatter!(ax2, bandgaps, IV_parameters["V_oc"])
ax3 = Axis(fig[2, 2], title="T = $T K", xlabel="Bandgap (eV)", ylabel="Short-Circuit Current (A/m^2)")
scatter!(ax3, bandgaps, IV_parameters["J_sc"])

# pick a specific bandgap to plot the IV curve 
bandgap = 1.2

# call the detailed balance function for this bandgap (by default, 100 points are plotted between V=0 and V=Eg for the IV curve, can be adjusted)
outputs = detailed_balance(spectrum_file=spectrum_file, T=T, E=bandgaps, num_voltages=500)

# decompose the outputs
E = outputs[1] # photon energy
spectrum_params = outputs[2] # spectrum parameters
IV_curves = outputs[3] # IV curves
IV_parameters = outputs[4] # IV parameters

ax4 = Axis(fig[3, 1:2], limits=(0, nothing, 0, maximum(IV_curves["J"][1])*1.5), title="T = $T K, Eg = $bandgap eV", xlabel="Voltage (V)", ylabel="Current (A/m^2)")
lines!(ax4, IV_curves["V"][1], IV_curves["J"][1], color=:blue, label="Bandgap = $bandgap eV")

fig