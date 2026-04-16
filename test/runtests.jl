using DifferentialGamesBaseSolvers
using Test

@testset "DifferentialGamesBaseSolvers.jl" begin
    include("../src/solvers/FNELQ/test/fnelq_tests.jl")
    include("../src/solvers/iLQGames/test/ilqgames_tests.jl")
end
