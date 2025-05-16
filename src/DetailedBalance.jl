module DetailedBalance

# import the packages we need
using Pkg, CSV, DataFrames, Trapz, Interpolations, Roots, ProgressBars

# export the functions we want to be accessible
export detailed_balance

include("irzi_functions.jl")
include("cell_functions.jl")

end
