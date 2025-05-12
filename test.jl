using Test

include("irzi_functions.jl")

@testset "rzi_lut" begin
    # test the rzi_lut function

    # test the output of the rzi_lut function with default lut_lim and num values
    luts = rzi_lut()
    @test size(luts) == (4,) # there should be 4 look-up tables for orders 1,2,3,4
    @testset "lut lower lim" for lut in luts # the value of the look-up tables at the lower limit should be 0
        @test lut(1.4e-3) == 0 
    end

end

@testset "rzi" begin
    # test the rzi function
end