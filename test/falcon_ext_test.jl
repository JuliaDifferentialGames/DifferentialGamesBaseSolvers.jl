# Tests for the FALCON soft dependency (ext/DifferentialGamesBaseSolversFALCONExt.jl).
# Skipped when FALCON is not installed in the test environment.
using Test
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

@testset "FALCON soft dependency" begin
    @testset "stubs error without FALCON" begin
        if Base.get_extension(DifferentialGamesBaseSolvers, :DifferentialGamesBaseSolversFALCONExt) === nothing
            @test_throws MethodError falcon()
        end
    end

    if Base.find_package("FALCON") !== nothing
        @eval using FALCON
        @testset "extension constructors" begin
            ext = Base.get_extension(DifferentialGamesBaseSolvers, :DifferentialGamesBaseSolversFALCONExt)
            @test ext !== nothing
            s = falcon(; ϵ_opt = 2e-2, ϵ_feas = 1e-4)
            @test s isa FALCON.FALCONSolver
            @test s.ϵ_opt == 2e-2 && s.ϵ_feas == 1e-4
            @test solver_capabilities(typeof(s)) == solver_capabilities(FALCON.FALCONSolver)
            @test :PD_GNEP in solver_capabilities(typeof(s))
        end
    end
end
