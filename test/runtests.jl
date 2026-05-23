using DifferentialGamesBaseSolvers
using Test

@testset "DifferentialGamesBaseSolvers.jl" begin
    include("../src/solvers/FNELQ/test/fnelq_tests.jl")
    include("../src/solvers/iLQGames/test/ilqgames_tests.jl")
    include("../src/solvers/ALGAMES/test/test_algames.jl")
    include("../src/solvers/LIBR/test/libr_tests.jl")
    include("../src/solvers/YIPAVEL/test/yipavel_tests.jl")
    include("../src/solvers/INVLQ/test/invlq_tests.jl")
    include("../src/solvers/RHN/test/rhn_tests.jl")
end
