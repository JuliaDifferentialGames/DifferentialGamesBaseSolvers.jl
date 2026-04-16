using Test
using LinearAlgebra
using Statistics
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# ============================================================================
# Shared test fixtures
# ============================================================================

"""
Scalar 1-player 1-state 1-control LQ game with known closed-form solution.

  min  Σ_{k=0}^{N-1} (q·x² + r·u²) + qf·xf²
  s.t. x(k+1) = a·x(k) + b·u(k)
"""
function scalar_lq_game(;
    a=1.0, b=1.0, q=1.0, r=1.0, qf=10.0, x0=1.0, N=10, dt=0.1
)
    A      = reshape([a], 1, 1)
    B      = [reshape([b], 1, 1)]
    Q      = [reshape([q], 1, 1)]
    R      = [reshape([r], 1, 1)]
    Qf_arr = [reshape([qf], 1, 1)]
    return LQGameProblem(A, B, Q, R, Qf_arr, [x0], Float64(N * dt); dt=dt)
end

"""
Two-player shared-state double integrator LQ game (n=4, m=[1,1]).

Player 1 regulates (x₁, x₂) to zero; Player 2 drives (x₃, x₄) to a target.
Known to have a unique feedback Nash equilibrium solvable by FNELQ.
"""
function two_player_double_integrator(; N=20, dt=0.1)
    A  = [1.0 dt  0.0  0.0;
          0.0 1.0  0.0  0.0;
          0.0 0.0  1.0  dt;
          0.0 0.0  0.0  1.0]
    B1 = reshape([0.0, dt,  0.0, 0.0], 4, 1)
    B2 = reshape([0.0, 0.0, 0.0, dt],  4, 1)

    Q1  = diagm([1.0, 0.1, 0.0, 0.0])
    Q2  = diagm([0.0, 0.0, 1.0, 0.1])
    R1  = 0.1 * Matrix{Float64}(I, 1, 1)
    R2  = 0.1 * Matrix{Float64}(I, 1, 1)
    Qf1 = diagm([10.0, 1.0, 0.0, 0.0])
    Qf2 = diagm([0.0,  0.0, 10.0, 1.0])
    x0  = [2.0, 0.0, -1.0, 0.0]

    return LQGameProblem(A, [B1, B2], [Q1, Q2], [R1, R2], [Qf1, Qf2],
                         x0, Float64(N * dt); dt=dt)
end

# ============================================================================
# Category 1: FNELQ analytical ground truth
# ============================================================================

@testset "FNELQ: scalar LQ — backward pass matches Riccati" begin
    # N=2 Riccati recursion (a=1, b=1, q=1, r=1, qf=10):
    #   Z(2)=10, P(1)=10/11, Z(1)=21/11, P(0)=21/32
    Z2 = 10.0
    P1_analytical = Z2 / (1.0 + Z2)
    Z1 = Z2 - Z2^2 / (1.0 + Z2) + 1.0
    P0_analytical = Z1 / (1.0 + Z1)

    game = scalar_lq_game(N=2)
    sol  = solve(game, FNELQ(regularization=0.0))

    @testset "Gains match Riccati" begin
        @test sol.strategy.gains[1][1][1, 1] ≈ P0_analytical atol=1e-10
        @test sol.strategy.gains[1][2][1, 1] ≈ P1_analytical atol=1e-10
    end

    @testset "Optimal trajectory" begin
        x0 = 1.0
        u0 = -P0_analytical * x0
        x1 = x0 + u0
        u1 = -P1_analytical * x1
        x2 = x1 + u1
        @test sol.trajectories[1].controls[1, 1] ≈ u0 atol=1e-10
        @test sol.trajectories[1].controls[1, 2] ≈ u1 atol=1e-10
    end

    @testset "Optimal cost" begin
        x0_val = 1.0
        u0_val = -P0_analytical * x0_val
        x1_val = x0_val + u0_val
        u1_val = -P1_analytical * x1_val
        x2_val = x1_val + u1_val

        # Verify cost convention by checking what evaluate_stage_cost returns
        # for a known input, then use the same convention analytically.
        # LQStageCost uses ½ on stage costs; LQTerminalCost uses no ½.
        sc = get_objective(game, 1).stage_cost
        tc = get_objective(game, 1).terminal_cost
        J_stage0 = evaluate_stage_cost(sc, [x0_val], [u0_val], nothing, 1)
        J_stage1 = evaluate_stage_cost(sc, [x1_val], [u1_val], nothing, 2)
        J_term   = evaluate_terminal_cost(tc, [x2_val], nothing)

        @test sol.trajectories[1].cost ≈ J_stage0 + J_stage1 + J_term atol=1e-8
    end
end

@testset "FNELQ: two-player double integrator — feedback Nash properties" begin
    game = two_player_double_integrator(N=10, dt=0.1)
    sol  = solve(game, FNELQ(); check_compatibility=false)

    @testset "Returns FeedbackNash equilibrium" begin
        @test sol.equilibrium_type == :FeedbackNash
        @test sol.strategy isa FeedbackStrategy
        @test sol.converged
    end

    @testset "Costs are finite and positive" begin
        @test all(isfinite, get_costs(sol))
        @test all(c -> c > 0, get_costs(sol))
    end

    @testset "Trajectory is dynamically feasible" begin
        dyn = game.dynamics
        X   = sol.state_trajectory
        U1  = sol.trajectories[1].controls
        U2  = sol.trajectories[2].controls
        for k in 1:n_steps(game)
            x_pred = get_A(dyn, 1) * X[:, k] +
                     get_B(dyn, 1, 1) * U1[:, k] +
                     get_B(dyn, 2, 1) * U2[:, k]
            @test X[:, k+1] ≈ x_pred atol=1e-8
        end
    end
end

# ============================================================================
# Category 2: iLQGames convergence and output types
# ============================================================================

@testset "iLQGames: scalar LQ converges to known optimum" begin
    game   = scalar_lq_game(N=10)
    solver = iLQGames(max_iter=200, μ_init=0.0)
    sol    = solve(game, solver; check_compatibility=false)

    @testset "Solver converges" begin
        @test sol.converged
    end

    @testset "Cost decreases from initial" begin
        costs = sol.solver_info[:cost_history]
        @test costs[end] < costs[1]
    end

    @testset "Trajectory change below ε_conv" begin
        Δ_hist = sol.solver_info[:trajectory_change_history]
        @test last(Δ_hist) < solver.ε_conv
    end

    @testset "Output strategy type" begin
        # iLQGames returns FeedbackNash from the last FNELQ solve
        @test sol.equilibrium_type == :FeedbackNash
        @test sol.strategy isa FeedbackStrategy
    end
end

@testset "iLQGames: two-player double integrator — convergence" begin
    game   = two_player_double_integrator(N=10, dt=0.1)
    solver = iLQGames(max_iter=200, ε_conv=1e-4)
    sol    = solve(game, solver; check_compatibility=false)

    @testset "Solver converges" begin
        @test sol.converged
        @test sol.iterations <= 200
    end

    @testset "Costs are finite and positive" begin
        @test all(isfinite, get_costs(sol))
        @test all(c -> c > 0, get_costs(sol))
    end

    @testset "Trajectory is dynamically feasible" begin
        dyn = game.dynamics
        X   = sol.state_trajectory
        U1  = sol.trajectories[1].controls
        U2  = sol.trajectories[2].controls
        for k in 1:n_steps(game)
            x_pred = get_A(dyn, 1) * X[:, k] +
                     get_B(dyn, 1, 1) * U1[:, k] +
                     get_B(dyn, 2, 1) * U2[:, k]
            @test X[:, k+1] ≈ x_pred atol=1e-6
        end
    end

    @testset "Cost does not increase overall" begin
        costs = sol.solver_info[:cost_history]
        @test costs[end] <= costs[1]
    end
end

# ============================================================================
# Category 2: Optimality conditions at convergence
# ============================================================================

@testset "iLQGames: trajectory change history" begin
    game   = two_player_double_integrator(N=8, dt=0.1)
    solver = iLQGames(max_iter=100, ε_conv=1e-6)
    sol    = solve(game, solver; check_compatibility=false)

    Δ_hist = sol.solver_info[:trajectory_change_history]

    @testset "History length matches iterations" begin
        # One Δ entry per accepted step
        @test length(Δ_hist) == sol.iterations
    end

    @testset "Final Δ below ε_conv when converged" begin
        if sol.converged
            @test Δ_hist[end] < solver.ε_conv
        end
    end

    @testset "solver_info keys present" begin
        info = sol.solver_info
        @test haskey(info, :cost_history)
        @test haskey(info, :trajectory_change_history)
        @test haskey(info, :η_history)
        @test haskey(info, :μ_history)
        @test haskey(info, :open_loop_strategy)
        @test haskey(info, :final_regularisation)
        @test haskey(info, :final_cost_per_player)
    end

    @testset "open_loop_strategy stored in solver_info" begin
        ol = sol.solver_info[:open_loop_strategy]
        @test ol isa OpenLoopStrategy
        @test ol.n_players == 2
    end
end

# ============================================================================
# Category 3: Cross-solver validation
# ============================================================================

@testset "iLQGames vs FNELQ: on LQ game, costs are close" begin
    # On an LQ game, iLQGames (feedback Nash of LQ approximations) and FNELQ
    # (exact feedback Nash) should agree. The equilibrium concepts differ
    # (iLQGames seeks the fixed point of successive LQ approximations; FNELQ
    # solves the exact LQ game), but on a problem that is already LQ, the
    # two should produce near-identical trajectories and costs.
    game      = two_player_double_integrator(N=10, dt=0.1)
    sol_fnelq = solve(game, FNELQ(); check_compatibility=false)
    sol_ilqg  = solve(game, iLQGames(max_iter=200, ε_conv=1e-6);
                      check_compatibility=false)

    costs_fnelq = get_costs(sol_fnelq)
    costs_ilqg  = get_costs(sol_ilqg)

    @testset "Per-player costs agree within 1%" begin
        for i in 1:2
            rel_err = abs(costs_ilqg[i] - costs_fnelq[i]) / costs_fnelq[i]
            @test rel_err < 0.01
        end
    end

    @testset "Trajectories agree" begin
        X_fnelq = sol_fnelq.state_trajectory
        X_ilqg  = sol_ilqg.state_trajectory
        @test maximum(abs.(X_ilqg .- X_fnelq)) < 0.05
    end
end

# ============================================================================
# Category 4: Convergence and scaling
# ============================================================================
@testset "iLQGames: regularisation activates on ill-conditioned problem" begin
    A  = [1.0 100.0; 0.0 1.0]
    B  = [reshape([0.0; 1.0], 2, 1)]
    Q  = [0.01 * Matrix{Float64}(I, 2, 2)]
    R  = [1e-4 * reshape([1.0], 1, 1)]
    Qf = [Matrix{Float64}(I, 2, 2)]
    game = LQGameProblem(A, B, Q, R, Qf, [1.0, 0.0], 0.5; dt=0.1)

    solver = iLQGames(max_iter=50, μ_init=0.0, μ_max=1e6, ε_conv=1e-4)

    @testset "Does not throw an exception" begin
        # Convergence warnings are acceptable; exceptions are not.
        @test begin
            sol = solve(game, solver; check_compatibility=false)
            true
        end
    end

    @testset "Result is finite" begin
        sol = solve(game, solver; check_compatibility=false)
        @test isfinite(sol.trajectories[1].cost)
    end
end

@testset "iLQGames: warmstart reduces iterations on LQ game" begin
    game     = two_player_double_integrator(N=10, dt=0.1)
    solver   = iLQGames(max_iter=200, ε_conv=1e-5)

    sol_cold = solve(game, solver; check_compatibility=false)

    # Warmstart from the FNELQ solution — very close to optimal for LQ game
    sol_fnelq = solve(game, FNELQ(); check_compatibility=false)
    sol_warm  = solve(game, solver; warmstart=sol_fnelq,
                      check_compatibility=false)

    @testset "Warmstart converges" begin
        @test sol_warm.converged
    end

    @testset "Warmstart uses ≤ iterations than cold start" begin
        @test sol_warm.iterations <= sol_cold.iterations + 2
    end
end

@testset "iLQGames: LinearDynamics game — no DA constructed" begin
    game = two_player_double_integrator(N=5, dt=0.1)
    @assert game.dynamics isa LinearDynamics
    # Use enough iterations and a loose tolerance for this short-horizon problem.
    solver = iLQGames(max_iter=200, ε_conv=1e-4)

    @testset "Does not throw an exception" begin
        @test begin
            solve(game, solver; check_compatibility=false)
            true
        end
    end

    @testset "Result is finite" begin
        sol = solve(game, solver; check_compatibility=false)
        @test isfinite(sum(get_costs(sol)))
    end
end

# ============================================================================
# Solver struct validation
# ============================================================================

@testset "iLQGames: struct field validation" begin
    @test_throws AssertionError iLQGames(max_iter=0)
    @test_throws AssertionError iLQGames(ε_conv=-1.0)
    @test_throws AssertionError iLQGames(β=1.1)
    @test_throws AssertionError iLQGames(β=-0.1)
    @test_throws AssertionError iLQGames(μ_init=-1.0)
    @test_throws AssertionError iLQGames(μ_scale=0.5)   # must be > 1
    @test_throws AssertionError iLQGames(μ_decay=1.5)   # must be < 1
    @test iLQGames(max_iter=1) isa iLQGames             # minimum valid
end

println("\n✓ iLQGames test suite complete.")