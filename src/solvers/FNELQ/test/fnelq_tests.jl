using Test
using LinearAlgebra
using SparseArrays
using DifferentialGamesBase  
using DifferentialGamesBaseSolvers

@testset "Unconstrained LQ Game" begin
    # Define a basic game
    n, n_players = 4, 2
    A = randn(n, n)
    B = [randn(n, 2) for _ in 1:n_players]
    Q = [diagm(ones(n)) for _ in 1:n_players]
    R = [diagm(0.1 * ones(2)) for _ in 1:n_players]
    Qf = [diagm(10.0 * ones(n)) for _ in 1:n_players]

    game = UnconstrainedLQGame(A, B, Q, R, Qf, zeros(n), 10.0)

    # Solve via FNELQ
    solver = FNELQ()
    sol = solve(game, solver; verbose=true)

    @test sol.converged
end