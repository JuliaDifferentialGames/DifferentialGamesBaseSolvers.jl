using DifferentialGamesBaseSolvers
using Test

@testset "DifferentialGamesBaseSolvers.jl" begin
    include("../src/solvers/FNELQ/test/fnelq_tests.jl")
    include("../src/solvers/iLQGames/test/ilqgames_tests.jl")
    include("../src/solvers/ALGAMES/test/test_algames.jl")
end
