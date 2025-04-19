% Planck's constant
h = 6.62607015e-34;
% Speed of light
c = 299792458;
% Electron charge
qe = 1.602176634e-19;
% Boltzmann's constant
k = 1.380649e-23;
% Boltzmann's constant in eV
keV = k/qe;
% Cell temperature in K
T = 300;

% Bandgap in J

Eg = 1.2;
V = 0.9;

% Energies we want to integrate over
E = Eg:0.0001:5;

% Photon emission flux function
eFluxFunc = ( 2*pi/(h^3*c^2) ) * E.^2 ./ ( exp( (E-V)/(k*T) ) - 1); 

% Represent in terms of x = kT/(E-V)
integrand = 1/(exp(1/x))

plot(E, eFluxFunc);
ylim()

% Perform trapezoidal integration
% eFlux = trapz(Ee, eFluxFunc);