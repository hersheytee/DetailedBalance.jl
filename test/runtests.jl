using DetailedBalance
using Test

@testset "DetailedBalance.jl" begin
    # Test the detailed_balance_fluxes function
    # Define the path to the spectrum file
    spectrum_file = joinpath(@__DIR__, "data", "am0.csv")
    
    # Call the function with the test spectrum file and a temperature of 300 K
    outputs = detailed_balance(spectrum_file=spectrum_file, T=300)
    
    # Check if the outputs is a tuple with the expected length
    @test typeof(outputs) == Tuple{Vector{Float64}, Dict{String, Any}, Dict{String, Vector{Any}}, Dict{String, Vector}}
    
    println(outputs[1])

end
;