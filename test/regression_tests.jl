using DifferentialGamesBase
using DifferentialGamesBaseSolvers
using LinearAlgebra
using Test

# Regression tests for the affine (linear-term) FNELQ recursion and the iLQGames
# first-order consistency on exactly linear-quadratic problems.

@testset "FNELQ affine terms have the correct sign" begin
    # Scalar x⁺ = x + u over N = 2 steps; the second control is made prohibitively costly,
    # so u₁ minimises ½R u₁² + q_f (x₀ + u₁) with x₀ = 0, i.e. u₁ = -q_f/R.
    R, qf = 2.0, 3.0
    g = LTVLQGameProblem([ones(1,1), ones(1,1)], [[ones(1,1), ones(1,1)]],
                         [[zeros(1,1), zeros(1,1)]], [[fill(R,1,1), fill(1e9,1,1)]],
                         [zeros(1,1)], [0.0], 2.0; dt=1.0)
    obj = PlayerObjective(1, g.objectives[1].stage_cost, LQTerminalCost(zeros(1,1), [qf]))
    sol = solve(remake(g; objectives=[obj]), FNELQ())
    @test get_trajectory(sol, 1).controls[1] ≈ -qf / R atol=1e-6

    # Linear stage term r·u at the first step: u₁ = -r/R.
    sc = LQStageCost([zeros(1,1), zeros(1,1)], [fill(R,1,1), fill(1e9,1,1)],
                     [zeros(1,1), zeros(1,1)], [zeros(1), zeros(1)], [[1.0], [0.0]])
    obj2 = PlayerObjective(1, sc, LQTerminalCost(zeros(1,1)))
    sol2 = solve(remake(g; objectives=[obj2]), FNELQ())
    @test get_trajectory(sol2, 1).controls[1] ≈ -1.0 / R atol=1e-6
end

@testset "iLQGames on an exactly linear-quadratic PDGNEProblem approaches the FNELQ cost" begin
    dyn = (xi, ui, p, t) -> [xi[3]; xi[4]; ui[1]; ui[2]]
    stage    = DiagonalLQStageCost([0.0, 0.0, 1.0, 1.0], [0.1, 0.1])
    terminal = DiagonalLQTerminalCost([10.0, 10.0, 1.0, 1.0])
    p1 = PlayerSpec(1, 4, 2, [ 3.0, 0.0, 0.0, 0.0], dyn, PlayerObjective(1, stage, terminal))
    p2 = PlayerSpec(2, 4, 2, [-3.0, 0.0, 0.0, 0.0], dyn, PlayerObjective(2, stage, terminal))
    game = PDGNEProblem([p1, p2], 5.0, 0.1)
    sol  = solve(game, iLQGames(; max_iter=100, ε_conv=0.02))
    # Exact decoupled LQR value (ZOH): 7.9835 per player (see docs LQ tutorial)
    @test all(isapprox.(get_costs(sol), 7.9835; rtol=0.02))
    # costs must decrease monotonically from the zero-control initial cost
    ch = sol.solver_info[:cost_history]
    @test all(diff(ch) .< 1e-6)
end

@testset "iLQGames declares :SeparableDynamics" begin
    @test :SeparableDynamics in solver_capabilities(iLQGames)
end

@testset "LQGameProblem converting constructor" begin
    n = 4
    game = LQGameProblem(0.95 * I(n), [I(n)[:, 1:1], I(n)[:, 3:3]],
                         [diagm(ones(n)), diagm(ones(n))], [fill(0.1,1,1), fill(0.1,1,1)],
                         [10.0 * diagm(ones(n)), 10.0 * diagm(ones(n))], ones(n), 2.0; dt=0.1)
    strict = LQGameProblem(Matrix(0.95 * I(n)), [Matrix(1.0I(n))[:, 1:1], Matrix(1.0I(n))[:, 3:3]],
                         [Matrix(1.0I(n)), Matrix(1.0I(n))], [fill(0.1,1,1), fill(0.1,1,1)],
                         [Matrix(10.0I(n)), Matrix(10.0I(n))], ones(n), 2.0; dt=0.1)
    @test get_costs(solve(game, FNELQ())) ≈ get_costs(solve(strict, FNELQ()))
end

@testset "ALGAMES terminal_constraints enforces state constraints at x(N+1)" begin
    dyn = (xi, ui, p, t) -> ui
    function mkplayer(i, x0, goal)
        o = 2 * (i - 1)
        PlayerSpec(i, 2, 2, x0, dyn, PlayerObjective(i,
            NonlinearStageCost((x, u, p, t) -> 0.01 * sum(abs2, u); is_separable = true),
            NonlinearTerminalCost((x, p) -> 10.0 * ((x[o+1] - goal[1])^2 + (x[o+2] - goal[2])^2))),
            [control_bounds(i; control_offset = o, control_dim = 2, lower = [-2.0, -2.0], upper = [2.0, 2.0])])
    end
    col = SharedInequality([1, 2]; func = (x, u, p, t) -> [1.5^2 - sum(abs2, x[1:2] - x[3:4])], dim = 1)
    # The goals are only 1.0 apart, so the unconstrained terminal state violates the 1.5 separation.
    game = PDGNEProblem([mkplayer(1, [0.0, 2.0], [5.0, 0.5]), mkplayer(2, [0.0, -2.0], [5.0, -0.5])],
                        [col], 4.0, 0.1)
    dfinal(s) = norm(s.state_trajectory[1:2, end] - s.state_trajectory[3:4, end])

    s0 = solve(game, ALGAMES(; tol_opt = 1e-2))
    @test dfinal(s0) < 1.4                       # documented limitation: x(N+1) is unconstrained by default

    s1 = solve(game, ALGAMES(; tol_opt = 1e-2, terminal_constraints = true))
    @test s1.converged
    @test dfinal(s1) ≥ 1.5 - 1e-3                # terminal state now respects the constraint
    # control-dependent rows (the control bounds) are excluded at the terminal node and still hold
    @test maximum(abs, get_trajectory(s1, 1).controls) ≤ 2.0 + 1e-2
end

@testset "Receding-horizon warm-start shift handles joint-state trajectories (ALGAMES layout)" begin
    di(xi, ui, p, t) = [xi[1] + 0.1 * xi[2], xi[2] + 0.1 * ui[1]]
    stage(x, u, p, t) = 0.0
    players = [PlayerSpec(i, 2, 1, [Float64(i), 0.0], di,
                          PlayerObjective(i, NonlinearStageCost(stage), NonlinearTerminalCost((x, p) -> 0.0)))
               for i in 1:2]
    game = PDGNEProblem(players, 0.6, 0.1)
    N, n = 6, 4
    X = repeat(reshape(collect(1.0:(N + 1)), 1, :), n, 1)                 # joint state, n × (N+1)
    trajs = [Trajectory{Float64}(i, copy(X), fill(0.5, 1, N), collect(0.0:0.1:0.6), 0.0) for i in 1:2]
    shifted = DifferentialGamesBaseSolvers._shift_trajectories(trajs, game, 0.1)
    @test all(size(t.states) == (n, N + 1) for t in shifted)
    @test shifted[1].states[:, 1:N] == X[:, 2:N+1]
    # the appended node is the per-block dynamics step from the old last node
    last = X[:, N+1]
    @test shifted[1].states[1:2, N+1] ≈ [last[1] + 0.1 * last[2], last[2] + 0.1 * 0.5]
    @test shifted[1].states[3:4, N+1] ≈ [last[3] + 0.1 * last[4], last[4] + 0.1 * 0.5]
end
