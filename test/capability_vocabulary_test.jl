using DifferentialGamesBase
using DifferentialGamesBaseSolvers
using Test

# ============================================================================
# Enforces that every registered GameSolver's `solver_capabilities` stays
# within DifferentialGamesBase.CANONICAL_CAPABILITIES — closing the "ad hoc
# vocabulary" gap permanently rather than just once. A future solver that
# emits a new, un-registered symbol fails this test loudly instead of
# silently drifting the vocabulary out of sync with `required_capabilities`.
# ============================================================================

const ALL_SOLVER_TYPES = [
    FNELQ, SFNELQ, CCiLQGames, iLQGames, ALGAMES, LIBR, YiPavel,
    InverseLQGames, RecedingHorizonNash,
]

@testset "Capability vocabulary is canonical" begin
    for T in ALL_SOLVER_TYPES
        caps = solver_capabilities(T)
        @testset "$T" begin
            @test validate_capabilities(caps) === nothing
        end
    end
end
