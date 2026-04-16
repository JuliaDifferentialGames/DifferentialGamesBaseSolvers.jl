using Test
using LinearAlgebra
using Statistics
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# Reach into solver internals for white-box diagnostics.
# These are intentionally unexported — we access them by qualified name.
import DifferentialGamesBaseSolvers:
    _build_da, _total_costs, _line_search, _strategy_change, _solve

# Include the example (adjust path as needed)
include(joinpath(@__DIR__, "..", "examples", "figure8_game.jl"))

# ============================================================================
# Test 1: Dynamics correctness
# ============================================================================

@testset "Figure-8 dynamics: unicycle ODE" begin
    # At φ=0, v=1: ẋ should be [1, 0, ω, a]
    x = [0.0, 0.0, 0.0, 1.0]
    u = [0.5, 0.2]
    ẋ = unicycle_ode(x, u, nothing, 0.0)
    @test ẋ[1] ≈ 1.0    # v·cos(0) = 1
    @test ẋ[2] ≈ 0.0    # v·sin(0) = 0
    @test ẋ[3] ≈ 0.5    # φ̇ = ω
    @test ẋ[4] ≈ 0.2    # v̇ = a

    # At φ=π/2, v=2: ẋ should be [0, 2, ω, a]
    x2 = [0.0, 0.0, π/2, 2.0]
    ẋ2 = unicycle_ode(x2, u, nothing, 0.0)
    @test ẋ2[1] ≈ 0.0 atol=1e-14
    @test ẋ2[2] ≈ 2.0

    # CoupledNonlinearDynamics wraps the ODE correctly
    game = figure8_game()
    @test game.dynamics isa CoupledNonlinearDynamics
    @test total_state_dim(game.dynamics) == 4
    @test total_control_dim(game.dynamics) == 2
    @test game.metadata.control_dims == [1, 1]
end

# ============================================================================
# Test 2: Cost function evaluation
# ============================================================================

@testset "Figure-8 costs: stage cost values" begin
    game = figure8_game()

    # Player 1 cost: px² + py² + ω²
    x  = [3.0, 4.0, 0.0, 1.0]   # px=3, py=4, φ=0, v=1
    u1 = [2.0]                    # ω=2
    c1 = evaluate_stage_cost(get_objective(game, 1).stage_cost, x, u1, nothing, 1)
    @test c1 ≈ 3^2 + 4^2 + 2^2   # = 9 + 16 + 4 = 29

    # Player 2 cost: (v - v_ref)² + a²  with v_ref=1
    u2 = [0.3]                    # a=0.3
    c2 = evaluate_stage_cost(get_objective(game, 2).stage_cost, x, u2, nothing, 1)
    @test c2 ≈ (1.0 - F8_VREF)^2 + 0.3^2   # = 0 + 0.09

    # At origin with zero control: P1 cost = 0
    x0 = [0.0, 0.0, 0.0, 1.0]
    @test evaluate_stage_cost(get_objective(game, 1).stage_cost, x0, [0.0], nothing, 1) ≈ 0.0

    # Terminal costs are zero (iLQGames.jl reference uses no terminal cost)
    @test evaluate_terminal_cost(get_objective(game, 1).terminal_cost, x, nothing) ≈ 0.0
    @test evaluate_terminal_cost(get_objective(game, 2).terminal_cost, x, nothing) ≈ 0.0
end

# ============================================================================
# Test 3: Discrete approximation and rollout
# ============================================================================

@testset "Figure-8 discretisation and rollout" begin
    game = figure8_game()
    N    = n_steps(game)
    dt   = game.time_horizon.dt
    t_vec = collect(range(0.0, game.time_horizon.tf, length=N+1))

    # Zero-control rollout: unicycle drifts in a straight line at v=0.5
    U_zero = zeros(2, N)
    X_zero = rollout(game.dynamics, game.initial_state, U_zero, nothing, t_vec)

    @test X_zero[:, 1] ≈ game.initial_state
    # With φ=0, v=0.5, u=0: px increases by ≈ 0.5*dt per step
    @test X_zero[1, 2] ≈ 1.0 + 0.5 * dt   atol=1e-4   # px grows
    @test X_zero[2, 2] ≈ 1.0              atol=1e-4   # py constant
    @test X_zero[3, 2] ≈ 0.0              atol=1e-6   # φ constant
    @test X_zero[4, 2] ≈ 0.5              atol=1e-6   # v constant

    # After 20s at v=0.5: px ≈ 1 + 0.5*20 = 11
    @test X_zero[1, end] ≈ 1 + 0.5 * N * dt   atol=0.1

    # DA forward step matches rollout step
    da = _build_da(game, iLQGames())
    x0 = game.initial_state
    u0 = zeros(2)
    x1_da  = da_step(da, x0, u0, nothing, 0.0)
    x1_rol = X_zero[:, 2]
    @test x1_da ≈ x1_rol atol=1e-10
end

# ============================================================================
# Test 4: Warmstart trajectory quality
# ============================================================================

@testset "Figure-8 warmstart: trajectory covers expected range" begin
    game = figure8_game()
    X_w, U_w = figure8_warmstart(game)

    # Heading oscillates over ±π/2 → px and py should span several metres
    px_range = maximum(X_w[1, :]) - minimum(X_w[1, :])
    py_range = maximum(X_w[2, :]) - minimum(X_w[2, :])

    @testset "Warmstart spans ≥ 2 m in x and y" begin
        @test px_range > 2.0
        @test py_range > 2.0
    end

    @testset "Warmstart speed ramps toward v_ref" begin
        # After 3s (30 steps) speed should be close to 1.0
        @test X_w[4, 31] ≈ F8_VREF atol=0.1
    end

    @testset "Warmstart steering is sinusoidal" begin
        # ω = A·(2π/T)·cos(2πt/T) with T=10s over 20s → 2 full cycles → 4 zero-crossings
        u1 = vec(U_w[1, :])
        sign_changes = sum(u1[k] * u1[k+1] < 0 for k in 1:length(u1)-1)
        @test sign_changes >= 4
    end

    @testset "Warmstart trajectory is dynamically consistent" begin
        N    = n_steps(game)
        dt   = game.time_horizon.dt
        t_vec = collect(range(0.0, game.time_horizon.tf, length=N+1))
        X_check = rollout(game.dynamics, game.initial_state, U_w, nothing, t_vec)
        @test X_check ≈ X_w atol=1e-10
    end
end

# ============================================================================
# Test 5: Single FNELQ subgame at warmstart nominal
# ============================================================================

@testset "Figure-8 inner LQ subgame at warmstart" begin
    game = figure8_game()
    X_w, U_w = figure8_warmstart(game)

    da   = _build_da(game, iLQGames())
    exp  = expand(game, X_w, U_w, da)
    lq   = assemble_lq_game(exp, game)

    @testset "LQ game assembled" begin
        @test lq.dynamics isa LinearDynamics
        @test n_steps(lq) == n_steps(game)
    end

    lq_sol = solve(lq, FNELQ(); check_compatibility=false)

    @testset "FNELQ converges on warmstart subgame" begin
        @test lq_sol.converged
        @test all(isfinite, get_costs(lq_sol))
    end

    @testset "FNELQ subgame cost is lower than warmstart cost" begin
        # The LQ subgame optimal should improve on the warmstart nominal
        N   = n_steps(game)
        cd  = game.metadata.control_dims
        co  = [0; cumsum(cd)[1:end-1]]
        t_vec = collect(range(0.0, game.time_horizon.tf, length=N+1))

        costs_warm = [
            total_cost(get_objective(game, i),
                       [X_w[:, k] for k in 1:N+1],
                       [U_w[co[i]+1:co[i]+cd[i], k] for k in 1:N],
                       nothing)
            for i in 1:2
        ]

        @test lq_sol.trajectories[1].cost < costs_warm[1]
    end

    @testset "FNELQ feedback gains are non-trivial" begin
        N  = n_steps(game)
        strat = lq_sol.strategy
        P1_norms = [norm(get_gain(strat, 1, k)) for k in 1:N]
        P2_norms = [norm(get_gain(strat, 2, k)) for k in 1:N]
        @test maximum(P1_norms) > 1e-4
        @test maximum(P2_norms) > 1e-4
    end
end

# ============================================================================
# Test 6: Player 2 gradient direction — does it push toward v_ref?
# ============================================================================

@testset "Figure-8 Player 2 gradient pushes v toward v_ref" begin
    game = figure8_game()
    X_w, U_w = figure8_warmstart(game)

    da  = _build_da(game, iLQGames())
    exp = expand(game, X_w, U_w, da)

    # At each step k, the gradient gx[2][k][4] (∂ℓ₂/∂v) should be negative
    # when v < v_ref (i.e., push v upward)
    n_below_ref = sum(X_w[4, k] < F8_VREF for k in 1:n_steps(game))

    if n_below_ref > 0
        # Where v < v_ref, ∂ℓ₂/∂v = 2(v - v_ref) < 0 → gradient pushes v up
        gx2_v = [exp.costs.gx[2][k][4] for k in 1:n_steps(game)
                 if X_w[4, k] < F8_VREF]
        @test all(g -> g < 0, gx2_v)
    end

    # The gradient for Player 2's control ∂ℓ₂/∂a = 2a, should be ≈ 0 when a≈0
    gu2_late = [exp.costs.gu[2][k][1] for k in 50:n_steps(game)]
    # After the ramp, a≈0 so gu should be small
    @test all(abs.(gu2_late) .< 0.1)
end

# ============================================================================
# Test 7: Line search direction at warmstart
# ============================================================================

@testset "Figure-8 line search: η=1 step from warmstart decreases cost" begin
    game = figure8_game()
    X_w, U_w = figure8_warmstart(game)

    da   = _build_da(game, iLQGames())
    exp  = expand(game, X_w, U_w, da)
    lq   = assemble_lq_game(exp, game)

    fnelq  = FNELQ(check_singularity=true)
    lq_sol = _solve(lq, fnelq, nothing, false)
    strat  = lq_sol.strategy

    # Rollout with η=1
    X1, U1 = rollout_strategy(game.dynamics, game.initial_state, strat, nothing; η=1.0)

    # Compute costs
    N  = n_steps(game)
    cd = game.metadata.control_dims
    co = [0; cumsum(cd)[1:end-1]]

    cost_warm = sum(
        total_cost(get_objective(game, i),
                   [X_w[:, k] for k in 1:N+1],
                   [U_w[co[i]+1:co[i]+cd[i], k] for k in 1:N],
                   nothing)
        for i in 1:2
    )
    cost_step = sum(
        total_cost(get_objective(game, i),
                   [X1[:, k] for k in 1:N+1],
                   [U1[co[i]+1:co[i]+cd[i], k] for k in 1:N],
                   nothing)
        for i in 1:2
    )

    @testset "η=1 step is finite" begin
        @test all(isfinite, X1)
        @test all(isfinite, U1)
    end

    @testset "η=1 step decreases total cost from warmstart" begin
        @test cost_step < cost_warm
    end
end

# ============================================================================
# Test 8: Full solve produces figure-8 shape
# ============================================================================

@testset "Figure-8 Player 2 feedforward drives positive acceleration" begin
    game = figure8_game()
    X_w, U_w = figure8_warmstart(game)

    da   = _build_da(game, iLQGames())
    exp  = expand(game, X_w, U_w, da)
    lq   = assemble_lq_game(exp, game)
    fnelq  = FNELQ(check_singularity=true)
    lq_sol = _solve(lq, fnelq, nothing, false)
    strat  = lq_sol.strategy   # FeedbackStrategy

    N = n_steps(game)

    # α for player 2 at each step — should be negative (so -η·α > 0 = acceleration)
    α2 = [get_feedforward(strat, 2, k)[1] for k in 1:N]
    @testset "Player 2 feedforward α₂ is negative (→ positive acceleration)" begin
        # The mean α₂ should be negative so that u₂ = -P₂x - η·α₂ > 0
        @test mean(α2[1:30]) < 0
    end

    # Nominal controls for player 2 from FNELQ forward pass
    u2_nom = lq_sol.trajectories[2].controls
    @testset "FNELQ nominal u₂ is positive (acceleration)" begin
        @test mean(u2_nom[1, 1:30]) > 0
    end

    # After one η=1 rollout from warmstart with this strategy,
    # speed should increase toward v_ref
    X1, U1 = rollout_strategy(game.dynamics, game.initial_state, strat, nothing; η=1.0)
    @testset "η=1 rollout: speed increases from v₀" begin
        @test mean(X1[4, 10:30]) > game.initial_state[4]
    end
end

@testset "Figure-8 full solve: trajectory shape" begin
    sol = solve_figure8(verbose=false)

    @testset "Converged" begin
        @test sol.converged
        @test sol.equilibrium_type == :OpenLoopNash
    end

    X = sol.state_trajectory

    @testset "Trajectory spans expected range (figure-8 lobes ≥ 2m)" begin
        px_range = maximum(X[1, :]) - minimum(X[1, :])
        py_range = maximum(X[2, :]) - minimum(X[2, :])
        @test px_range > 2.0
        @test py_range > 2.0
    end

    @testset "Speed converges toward v_ref = $(F8_VREF)" begin
        # After transient (first 30 steps), mean speed should be near v_ref.
        # Requires w_a << w_v so Player 2 prefers tracking over coasting.
        v_mean = mean(X[4, 30:end])
        @test v_mean ≈ F8_VREF atol=0.3
    end

    @testset "Steering oscillates (figure-8 requires repeated loops)" begin
        u1 = vec(sol.trajectories[1].controls)
        sign_changes = sum(u1[k] * u1[k+1] < 0 for k in 1:length(u1)-1)
        @test sign_changes >= 4
    end

    @testset "Player 2 cost is finite and positive" begin
        @test sol.trajectories[2].cost > 0
        @test isfinite(sol.trajectories[2].cost)
    end
end

println("\n✓ Figure-8 diagnostic tests complete.")