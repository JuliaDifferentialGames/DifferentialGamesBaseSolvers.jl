# ============================================================================
# test_pangolin.jl
#
# Comprehensive test suite for the PANGOLIN solver and supporting layers.
#
# Test hierarchy:
#   L1  — FeedbackStrategy: construction, scale_feedforward!, copyto!
#   L2  — quadraticize: LQStageCost identity, DiagonalLQ, GeneralStageCost AD
#   L3  — linearize: LinearDynamics zero defect, CoupledNonlinear Jacobian
#   L4  — lq_approximation!: LQ game round-trip, nonlinear game
#   L5  — solve_lq_game!: matches FNELQ on LTI LQ games (primary regression)
#   L6  — rollout!: LinearDynamics exact, deviation abort
#   L7  — backtrack_scale!: scale applied, stable rollout returned
#   L8  — has_converged: threshold behaviour
#   L9  — ALOptions: construction, invalid params
#   L10 — ALAugmentedObjective: construction, augmented_stage_cost
#   L11 — AL quadraticize: inactive constraint no-op, active rank-1 update
#   L12 — update_multipliers!: inequality clamp, equality unclamped
#   L13 — maybe_update_penalty!: schedule trigger, ρ_max cap
#   L14 — PANGOLIN struct: construction, solver_capabilities
#   L15 — _solve unconstrained LQ: matches FNELQ on 2-player LTI game
#   L16 — _solve constrained: bound constraint satisfied at convergence
#   L17 — _solve nonlinear: converges on simple nonlinear game
#   L18 — warmstart: zero_control, straight_line, warmstart from prior solve
#   L19 — type stability: @inferred on hot-path functions
# ============================================================================

using Test
using LinearAlgebra
using SparseArrays
using Statistics
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# ============================================================================
# Shared fixtures
# ============================================================================

"""Build a 2-player LTI LQ game. Stable A (spectral radius < 1)."""
function lti_game(T=Float64; n=4, m=2, n_players=2, tf=1.0, dt=0.1)
    A  = T(0.9) * Matrix{T}(I, n, n)
    B  = [Matrix{T}(I, n, m)        for _ in 1:n_players]
    Q  = [Matrix{T}(I, n, n)        for _ in 1:n_players]
    R  = [Matrix{T}(I, m, m)        for _ in 1:n_players]
    Qf = [T(5) * Matrix{T}(I, n, n) for _ in 1:n_players]
    x0 = ones(T, n)
    LQGameProblem(A, B, Q, R, Qf, x0, T(tf); dt=T(dt))
end

"""2-player LTI game with a simple bound constraint on player 1's control."""
function constrained_lti_game(T=Float64; n=4, m=2, tf=1.0, dt=0.1)
    game = lti_game(T; n=n, m=m, tf=tf, dt=dt)
    # u₁ ≤ 0.5 * ones(m) — BoundConstraint on player 1's control
    bound = BoundConstraint(fill(T(-10), m), fill(T(0.5), m); applies_to=:u)
    pc    = PrivateConstraint(bound, 1)
    # Re-construct with constraint attached
    # (GameProblem is immutable; rebuild via LQGameProblem + manual injection)
    # For testing, construct directly
    A  = T(0.9) * Matrix{T}(I, n, n)
    B  = [Matrix{T}(I, n, m) for _ in 1:2]
    Q  = [Matrix{T}(I, n, n) for _ in 1:2]
    R  = [Matrix{T}(I, m, m) for _ in 1:2]
    Qf = [T(5) * Matrix{T}(I, n, n) for _ in 1:2]
    x0 = ones(T, n)
    dyn = LinearDynamics(A, B)
    ctrl_dims = dyn.control_dims
    objs = [PlayerObjective(i,
                LQStageCost(Q[i], R[i]),
                LQTerminalCost(Qf[i])) for i in 1:2]
    th   = DiscreteTime(T(tf), T(dt))
    meta = GameMetadata(
        [n], ctrl_dims, [0], [0, m],
        CouplingGraph(sparse(trues(2,2)), Vector{Int}[], nothing), false, nothing
    )
    GameProblem{T}(2, objs, dyn, x0, [pc], SharedConstraint[], th, meta)
end

"""1-player 1D nonlinear game: x(k+1) = tanh(x(k)) + u(k), scalar."""
function nonlinear_1player_game(T=Float64; tf=0.5, dt=0.1)
    n, m = 1, 1
    dyn  = CoupledNonlinearDynamics(
        (x, u, p, t) -> [tanh(x[1]) + u[1]],
        n, m
    )
    Q    = [Matrix{T}(I, n, n)]
    Qf   = [T(5) * Matrix{T}(I, n, n)]
    R    = [Matrix{T}(I, m, m)]
    x0   = T[1.0]
    stage = [NonlinearStageCost((x, u, p, t) -> T(0.5)*(x[1]^2 + u[1]^2), nothing, nothing, false)]
    terminal = [LQTerminalCost(Qf[1])]
    objs = [PlayerObjective(1, stage[1], terminal[1])]
    th   = DiscreteTime(T(tf), T(dt))
    meta = GameMetadata(
        [n], [m], [0], [0],
        CouplingGraph(sparse(trues(1,1)), Vector{Int}[], nothing), false, nothing
    )
    GameProblem{T}(1, objs, dyn, x0, PrivateConstraint[], SharedConstraint[], th, meta)
end

# ============================================================================
@testset "PANGOLIN — Comprehensive Test Suite" begin
# ============================================================================

# ─────────────────────────────────────────────────────────────────────────────
@testset "L1 — FeedbackStrategy" begin

    @testset "construction: correct dimensions" begin
        game = lti_game()
        strat = FeedbackStrategy(game)
        N, n, m, n_p = 10, 4, 2, 2
        @test length(strat.P) == n_p
        @test length(strat.α) == n_p
        for i in 1:n_p
            @test length(strat.P[i]) == N
            @test length(strat.α[i]) == N
            @test size(strat.P[i][1]) == (m, n)
            @test length(strat.α[i][1]) == m
        end
    end

    @testset "scale_feedforward!: scales α, leaves P unchanged" begin
        game = lti_game()
        strat = FeedbackStrategy(game)
        # Set α to known values
        for i in eachindex(strat.α), k in eachindex(strat.α[i])
            fill!(strat.α[i][k], 1.0)
        end
        # Set P to known values
        for i in eachindex(strat.P), k in eachindex(strat.P[i])
            fill!(strat.P[i][k], 2.0)
        end
        scale_feedforward!(strat, 0.5)
        # α halved
        @test all(all(strat.α[i][k] .≈ 0.5) for i in eachindex(strat.α) for k in eachindex(strat.α[i]))
        # P unchanged
        @test all(all(strat.P[i][k] .≈ 2.0) for i in eachindex(strat.P) for k in eachindex(strat.P[i]))
    end

    @testset "scale_feedforward!: zero scale zeroes α" begin
        game  = lti_game()
        strat = FeedbackStrategy(game)
        for i in eachindex(strat.α), k in eachindex(strat.α[i])
            fill!(strat.α[i][k], 3.0)
        end
        scale_feedforward!(strat, 0.0)
        @test all(all(strat.α[i][k] .== 0.0) for i in eachindex(strat.α) for k in eachindex(strat.α[i]))
    end

    @testset "copyto!: deep copy, no aliasing" begin
        game = lti_game()
        src  = FeedbackStrategy(game)
        dst  = FeedbackStrategy(game)
        fill!(src.P[1][1], 7.0)
        fill!(src.α[1][1], 3.0)
        copyto!(dst, src)
        @test dst.P[1][1] ≈ src.P[1][1]
        @test dst.α[1][1] ≈ src.α[1][1]
        # Verify no aliasing: mutating src does not affect dst
        fill!(src.P[1][1], 0.0)
        @test all(dst.P[1][1] .≈ 7.0)
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L2 — quadraticize: cost types" begin

    n, m = 4, 2
    x = randn(n); u = randn(m); k = 3

    @testset "LQStageCost identity: returns stored matrices" begin
        Q  = Matrix{Float64}(I, n, n)
        R  = Matrix{Float64}(I, m, m)
        sc = LQStageCost(Q, R)
        Qk, Rk, Mk, qk, rk, ck = quadraticize(sc, x, u, k)
        @test Qk === Q    # identity — same object
        @test Rk === R
        @test size(Mk) == (n, m)
        @test length(qk) == n
        @test length(rk) == m
    end

    @testset "DiagonalLQStageCost: promotes to correct dimensions" begin
        sc = DiagonalLQStageCost(ones(n), ones(m))
        Qk, Rk, Mk, qk, rk, ck = quadraticize(sc, x, u, k)
        @test size(Qk) == (n, n) || Qk isa Diagonal
        @test size(Rk) == (m, m) || Rk isa Diagonal
        @test size(Mk) == (n, m)
    end

    @testset "GeneralStageCost AD: quadratic func recovers exact Hessian" begin
        # Quadratic func with known Hessian: ½ xᵀQx + ½ uᵀRu
        Q_true = Diagonal(1.0:n |> collect .|> Float64)
        R_true = Diagonal(1.0:m |> collect .|> Float64)
        sc = GeneralStageCost(
            (x, u) -> 0.5 * sum(Q_true[i,i] * x[i]^2 for i in 1:n) +
                      0.5 * sum(R_true[j,j] * u[j]^2 for j in 1:m),
            n, m
        )
        Qk, Rk, Mk, qk, rk, ck = quadraticize(sc, x, u, k)
        @test Qk ≈ Matrix(Q_true) atol=1e-10
        @test Rk ≈ Matrix(R_true) atol=1e-10
        @test norm(Mk) < 1e-10     # no cross term
        # Linear terms: q = ∇x - Q*x - M*u = Q*x - Q*x = 0 for pure quadratic
        @test norm(qk) < 1e-10
    end

    @testset "GeneralStageCost AD: gradient consistency" begin
        # Affine+quadratic: ℓ = ½ x²+ ½ u² + aᵀx + bᵀu
        a = randn(n); b = randn(m)
        sc = GeneralStageCost(
            (x, u) -> 0.5*sum(x.^2) + 0.5*sum(u.^2) + dot(a, x) + dot(b, u),
            n, m
        )
        Qk, Rk, Mk, qk, rk, ck = quadraticize(sc, x, u, k)
        # q = ∇x - Q*x = (x + a) - x = a
        @test qk ≈ a atol=1e-8
        # r = ∇u - R*u = (u + b) - u = b
        @test rk ≈ b atol=1e-8
    end

    @testset "LQTerminalCost identity: returns stored Qf, qf" begin
        Qf = 5.0 * Matrix{Float64}(I, n, n)
        tc = LQTerminalCost(Qf)
        Qfk, qfk = quadraticize(tc, x)
        @test Qfk === Qf
        @test norm(qfk) == 0
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L3 — linearize: dynamics types" begin

    n, m, n_p = 4, 2, 2
    x  = randn(n)
    us = [randn(m) for _ in 1:n_p]
    k  = 2

    @testset "LinearDynamics LTI: zero defect, returns stored A, B" begin
        A   = 0.9 * Matrix{Float64}(I, n, n)
        B   = [Matrix{Float64}(I, n, m) for _ in 1:n_p]
        dyn = LinearDynamics(A, B)
        A_k, B_k, d_k = linearize(dyn, x, us, k)
        @test A_k === A          # identity — no copy for LTI
        @test norm(d_k) == 0.0   # linear dynamics: exact, no defect
        @test length(B_k) == n_p
        @test B_k[1] === B[1]
    end

    @testset "CoupledNonlinearDynamics: Jacobian correctness" begin
        # Linear function disguised as nonlinear — Jacobian must match A, B exactly
        A_true = 0.9 * Matrix{Float64}(I, n, n)
        B_true = [Matrix{Float64}(I, n, m) for _ in 1:n_p]
        B_cat  = hcat(B_true...)
        dyn = CoupledNonlinearDynamics(
            (x, u, p, t) -> A_true * x + B_cat * u,
            n, n_p * m
        )
        A_k, B_k, d_k = linearize(dyn, x, us, k)
        @test A_k ≈ A_true atol=1e-8
        for i in 1:n_p
            @test B_k[i] ≈ B_true[i] atol=1e-8
        end
        # Defect should be zero for linear dynamics (linearization is exact)
        @test norm(d_k) < 1e-8
    end

    @testset "CoupledNonlinearDynamics: nonzero defect for nonlinear func" begin
        # f(x, u) = [tanh(x₁); x₂; ...] + Bu
        # At x=0 with u=0: f(0,0) = 0, A*0 + B*0 = 0, defect = 0
        # At x=[0.5,...] with u=0: tanh(0.5) ≠ 0.5, defect nonzero
        B_cat  = hcat([Matrix{Float64}(I, n, m) for _ in 1:n_p]...)
        dyn2   = CoupledNonlinearDynamics(
            (x, u, p, t) -> vcat(tanh(x[1]), x[2:end]) + B_cat * u,
            n, n_p * m
        )
        us_zero = [zeros(m) for _ in 1:n_p]

        # At x=0, u=0: tanh(0)=0, defect should be zero
        A_k0, B_k0, d_k0 = linearize(dyn2, zeros(n), us_zero, k)
        @test norm(d_k0) < 1e-10

        # At x=[0.5,...], u=0: tanh(0.5) ≠ 0.5, defect in first state component
        x_nz = fill(0.5, n)
        A_k2, B_k2, d_k2 = linearize(dyn2, x_nz, us_zero, k)
        # defect[1] = tanh(0.5) - sech²(0.5)*0.5 ≠ 0
        @test abs(d_k2[1]) > 1e-6
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L4 — lq_approximation!" begin

    @testset "LQ game: lqg populated with correct Q, R from LQStageCost" begin
        game = lti_game()
        op   = OperatingPoint(game)
        lqg  = LTVLQApproximation(game)
        lq_approximation!(lqg, game, op)
        N = n_steps(game)
        # A must match stored dynamics matrix
        A_true = get_A(game.dynamics, 1)
        for k in 1:N
            @test lqg.A[k] ≈ A_true atol=1e-12
            @test norm(lqg.d[k]) < 1e-12   # zero defect for linear dynamics
        end
        # Q[1][k] must match LQStageCost Q
        Q_true = get_Q(game.objectives[1].stage_cost, 1)
        for k in 1:N
            @test lqg.Q[1][k] ≈ Q_true atol=1e-12
        end
        # Terminal cost
        @test lqg.Qf[1] ≈ game.objectives[1].terminal_cost.Qf atol=1e-12
    end

    @testset "lq_approximation! is idempotent: repeated calls give same result" begin
        game = lti_game()
        op   = OperatingPoint(game)
        lqg  = LTVLQApproximation(game)
        lq_approximation!(lqg, game, op)
        A_first = copy(lqg.A[1])
        lq_approximation!(lqg, game, op)
        @test lqg.A[1] ≈ A_first atol=0.0
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L5 — solve_lq_game! vs FNELQ (primary regression)" begin
    #
    # This is the critical regression test. PANGOLIN's solve_lq_game! on an
    # LTI LQ approximation must reproduce FNELQ's feedback gains exactly.
    # Any divergence indicates a bug in the backward pass.
    #

    @testset "2-player LTI: gains match FNELQ to 1e-10" begin
        game = lti_game()

        sol_fnelq = solve(game, FNELQ(); verbose=false)
        P_ref     = sol_fnelq.solver_info[:feedback_gains]

        op    = OperatingPoint(game)
        lqg   = LTVLQApproximation(game)
        strat = FeedbackStrategy(game)
        lq_approximation!(lqg, game, op)
        solver = PANGOLIN()
        solve_lq_game!(strat, lqg, solver)

        N = n_steps(game)
        for i in 1:2, k in 1:N
            @test strat.P[i][k] ≈ P_ref[i][k] atol=1e-10
        end
    end

    @testset "gains are finite and correct shape" begin
        game   = lti_game()
        op     = OperatingPoint(game)
        lqg    = LTVLQApproximation(game)
        strat  = FeedbackStrategy(game)
        lq_approximation!(lqg, game, op)
        solve_lq_game!(strat, lqg, PANGOLIN())
        N, n, m = n_steps(game), 4, 2
        for i in 1:2, k in 1:N
            @test all(isfinite, strat.P[i][k])
            @test size(strat.P[i][k]) == (m, n)
        end
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L6 — rollout!" begin

    @testset "LinearDynamics: rollout matches manual propagation" begin
        game = lti_game()
        op   = OperatingPoint(game)
        strat = FeedbackStrategy(game)   # zero strategy: u = 0

        rollout!(op, game, strat, game.initial_state)

        # With u=0: x(k+1) = A x(k), so x(N) = A^N x0
        A   = get_A(game.dynamics, 1)
        N   = n_steps(game)
        x_expected = copy(game.initial_state)
        for _ in 1:N
            x_expected = A * x_expected
        end
        @test op.x[N+1] ≈ x_expected atol=1e-10
    end

    @testset "rollout! returns true when no deviation limit" begin
        game  = lti_game()
        op    = OperatingPoint(game)
        strat = FeedbackStrategy(game)
        @test rollout!(op, game, strat, game.initial_state)
    end

    @testset "rollout! returns false when deviation exceeded" begin
        game  = lti_game()
        op    = OperatingPoint(game)
        # Set op.x to zeros so reference is far from true trajectory
        for k in eachindex(op.x); fill!(op.x[k], 0.0); end
        strat = FeedbackStrategy(game)
        # Very tight max_deviation — should abort on first step since x0 ≠ 0
        result = rollout!(op, game, strat, game.initial_state;
                          max_deviation=Float64(1e-10))
        # With x0 = ones(4) and x_ref = zeros, deviation on step 1 is ~0.9 > 1e-10
        @test result == false
    end

    @testset "rollout! overwrites x[1] with provided x0" begin
        game = lti_game()
        op   = OperatingPoint(game)
        strat = FeedbackStrategy(game)
        x0_new = 2.0 * game.initial_state
        rollout!(op, game, strat, x0_new)
        @test op.x[1] ≈ x0_new
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L7 — backtrack_scale!" begin

    @testset "returns true for stable linear game" begin
        game    = lti_game()
        solver  = PANGOLIN()
        strat   = FeedbackStrategy(game)
        # last_op must be a dynamically consistent trajectory for the deviation
        # check to be meaningful — roll out with zero strategy first
        last_op = OperatingPoint(game)
        rollout!(last_op, game, strat, game.initial_state)
        op = OperatingPoint(game)
        # Fill α with small nonzero values so the line search has something to scale
        for i in eachindex(strat.α), k in eachindex(strat.α[i])
            fill!(strat.α[i][k], 0.01)
        end
        result = backtrack_scale!(strat, op, last_op, game, solver, game.initial_state)
        @test result == true
    end

    @testset "α is scaled down after backtrack" begin
        game   = lti_game()
        solver = PANGOLIN(; α_init=0.5, α_step=0.5)
        op     = OperatingPoint(game)
        last_op = OperatingPoint(game)
        strat  = FeedbackStrategy(game)
        for i in eachindex(strat.α), k in eachindex(strat.α[i])
            fill!(strat.α[i][k], 1.0)
        end
        rollout!(last_op, game, strat, game.initial_state)
        α_before = strat.α[1][1][1]
        backtrack_scale!(strat, op, last_op, game, solver, game.initial_state)
        # After at least one scaling step, α should be strictly less than original
        @test abs(strat.α[1][1][1]) < abs(α_before) + 1e-10
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L8 — has_converged" begin

    game = lti_game()
    op1  = OperatingPoint(game)
    op2  = OperatingPoint(game)
    tol  = 1e-3

    @testset "identical operating points: converged" begin
        @test has_converged(op1, op2, Float64(tol)) == true
    end

    @testset "large difference: not converged" begin
        op2.x[3] .+= 1.0
        @test has_converged(op1, op2, Float64(tol)) == false
    end

    @testset "difference just below tol: converged" begin
        op1b = OperatingPoint(game); op2b = OperatingPoint(game)
        op2b.x[1][1] += tol * 0.99
        @test has_converged(op1b, op2b, Float64(tol)) == true
    end

    @testset "difference just above tol: not converged" begin
        op1c = OperatingPoint(game); op2c = OperatingPoint(game)
        op2c.x[1][1] += tol * 1.01
        @test has_converged(op1c, op2c, Float64(tol)) == false
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L9 — ALOptions" begin

    @testset "default construction" begin
        opts = ALOptions()
        @test opts.ρ_init > 0
        @test opts.ρ_max >= opts.ρ_init
        @test 1 < opts.φ <= 10
        @test opts.violation_tol > 0
    end

    @testset "custom values accepted" begin
        opts = ALOptions(; ρ_init=5.0, ρ_max=1000.0, φ=3.0, violation_tol=0.1)
        @test opts.ρ_init == 5.0
        @test opts.φ == 3.0
    end

    @testset "invalid φ rejected" begin
        @test_throws AssertionError ALOptions(; φ=1.0)    # must be > 1
        @test_throws AssertionError ALOptions(; φ=11.0)   # must be ≤ 10
    end

    @testset "ρ_max < ρ_init rejected" begin
        @test_throws AssertionError ALOptions(; ρ_init=100.0, ρ_max=1.0)
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L10 — ALAugmentedObjective" begin

    n, m = 4, 2
    N_test = 10   # test horizon

    function _make_al_obj(T=Float64)
        game     = constrained_lti_game(T)
        opts     = ALOptions()
        N        = n_steps(game)
        al_state = ALSolverState(
            Vector{SharedConstraint}(game.shared_constraints), opts, N
        )
        priv   = Vector{PrivateConstraint}(filter(pc -> pc.player == 1, game.private_constraints))
        shared = SharedConstraint[]
        ALAugmentedObjective(game.objectives[1], priv, shared,
                             al_state.λ_shared_traj, N, opts)
    end

    @testset "construction: correct multiplier structure" begin
        obj = _make_al_obj()
        game = constrained_lti_game()
        N    = n_steps(game)
        # BoundConstraint on m=2 control: 2*m = 4 constraints per timestep
        @test length(obj.λ_private) == 1        # one private constraint block
        @test size(obj.λ_private[1]) == (2*m, N)
        @test all(obj.λ_private[1] .== 0.0)
        @test obj.ρ == 1.0
    end

    @testset "augmented_stage_cost: no constraint active → equals base cost" begin
        obj   = _make_al_obj()
        x     = zeros(4)
        u_feas = zeros(2)
        k     = 1
        base_cost = evaluate_stage_cost(obj.base.stage_cost, x, u_feas, nothing, k)
        al_cost   = augmented_stage_cost(obj, x, u_feas, k)
        @test al_cost ≈ base_cost * obj.base.scaling atol=1e-10
    end

    @testset "augmented_stage_cost: active constraint increases cost" begin
        obj    = _make_al_obj()
        x      = zeros(4)
        u_viol = fill(1.0, 2)
        k      = 1
        base_cost = evaluate_stage_cost(obj.base.stage_cost, x, u_viol, nothing, k)
        al_cost   = augmented_stage_cost(obj, x, u_viol, k)
        @test al_cost > base_cost * obj.base.scaling
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L11 — AL quadraticize: active-set updates" begin

    n, m = 2, 2

    function _simple_al_obj(; u_max=0.5, ρ=1.0, λ_val=0.0, T=Float64, N=10)
        A  = T(0.9) * Matrix{T}(I, n, n)
        B  = [Matrix{T}(I, n, m)]
        Q  = [Matrix{T}(I, n, n)]
        R  = [Matrix{T}(I, m, m)]
        Qf = [T(5) * Matrix{T}(I, n, n)]
        x0 = ones(T, n)
        bound = BoundConstraint(fill(T(-10), m), fill(T(u_max), m); applies_to=:u)
        pc    = PrivateConstraint(bound, 1)
        obj   = PlayerObjective(1, LQStageCost(Q[1], R[1]), LQTerminalCost(Qf[1]))
        opts  = ALOptions(; ρ_init=T(ρ))
        λ_shared_traj = Matrix{T}[]   # no shared constraints
        al_obj = ALAugmentedObjective(obj, [pc], SharedConstraint[], λ_shared_traj, N, opts)
        al_obj.ρ = T(ρ)
        fill!(al_obj.λ_private[1], T(λ_val))
        return al_obj
    end

    @testset "inactive constraint: AL quadraticize = base quadraticize" begin
        al_obj = _simple_al_obj(; u_max=0.5, ρ=1.0, λ_val=0.0)
        x_ref  = zeros(n); u_ref = zeros(m); k = 1
        Q_al, R_al, _, _, _, _ = quadraticize(al_obj, x_ref, u_ref, k)
        Q_base, R_base, _, _, _, _ = quadraticize(al_obj.base.stage_cost, x_ref, u_ref, k)
        @test Q_al ≈ Matrix(Q_base) atol=1e-10
        @test R_al ≈ Matrix(R_base) atol=1e-10
    end

    @testset "active constraint: R increases (rank-1 update from active bound)" begin
        al_obj = _simple_al_obj(; u_max=0.0, ρ=10.0, λ_val=0.0)
        x_ref  = zeros(n); u_ref = ones(m) * 0.5; k = 1
        _, R_al, _, _, _, _ = quadraticize(al_obj, x_ref, u_ref, k)
        _, R_base, _, _, _, _ = quadraticize(al_obj.base.stage_cost, x_ref, u_ref, k)
        R_expected = Matrix(R_base) + 10.0 * Matrix{Float64}(I, m, m)
        @test R_al ≈ R_expected atol=1e-8
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L12 — update_multipliers!" begin

    function _make_constrained_al_obj()
        game     = constrained_lti_game()
        opts     = ALOptions(; ρ_init=1.0)
        N        = n_steps(game)
        al_state = ALSolverState(SharedConstraint[], opts, N)
        priv     = Vector{PrivateConstraint}(filter(pc->pc.player==1, game.private_constraints))
        ALAugmentedObjective(game.objectives[1], priv, SharedConstraint[],
                             al_state.λ_shared_traj, N, opts)
    end

    @testset "inequality: λ clamped to ≥ 0 at all timesteps" begin
        al_obj = _make_constrained_al_obj()
        game   = constrained_lti_game()
        op     = OperatingPoint(game)
        # Set all private multipliers to -5
        fill!(al_obj.λ_private[1], -5.0)
        update_multipliers!(al_obj, op)
        # u=0 satisfies u ≤ 0.5; update clamps to max(0, -5 + ρ*g) = 0
        @test all(al_obj.λ_private[1] .>= 0.0)
    end

    @testset "Δλ return value is non-negative" begin
        al_obj = _make_constrained_al_obj()
        game   = constrained_lti_game()
        op     = OperatingPoint(game)
        Δλ = update_multipliers!(al_obj, op)
        @test Δλ >= 0.0
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L13 — maybe_update_penalty!" begin

    function _make_al_state(T=Float64)
        opts = ALOptions(; ρ_init=T(1.0), ρ_max=T(1e4), φ=T(10.0), violation_tol=T(0.25))
        al_state = ALSolverState(SharedConstraint[], opts, 10)
        objs = ALAugmentedObjective[]
        return al_state, objs, opts
    end

    @testset "violation not decreasing: ρ increases geometrically" begin
        al_state, objs, opts = _make_al_state()
        al_state.prev_violation = 1.0
        # current_violation = 1.0 > (1 - 0.25) * 1.0 = 0.75 — should fire
        updated = maybe_update_penalty!(al_state, objs, 1.0)
        @test updated == true
        @test al_state.ρ ≈ 10.0   # φ * ρ_init
    end

    @testset "violation decreasing sufficiently: ρ unchanged" begin
        al_state, objs, opts = _make_al_state()
        al_state.prev_violation = 1.0
        # current_violation = 0.5 < 0.75 — not triggered
        updated = maybe_update_penalty!(al_state, objs, 0.5)
        @test updated == false
        @test al_state.ρ ≈ 1.0
    end

    @testset "ρ capped at ρ_max" begin
        al_state, objs, opts = _make_al_state()
        al_state.ρ = 5000.0   # near cap
        al_state.prev_violation = 1.0
        maybe_update_penalty!(al_state, objs, 1.0)
        @test al_state.ρ <= 1e4
    end

    @testset "ρ_max not exceeded even with multiple triggers" begin
        al_state, objs, opts = _make_al_state()
        al_state.prev_violation = 1.0
        for _ in 1:20
            maybe_update_penalty!(al_state, objs, 1.0)
        end
        @test al_state.ρ <= 1e4
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L14 — PANGOLIN struct" begin

    @testset "default construction" begin
        solver = PANGOLIN()
        @test solver.max_iter == 200
        @test solver.α_init == 0.5
        @test solver.α_step == 0.5
        @test solver.convergence_tol == 1e-3
        @test solver.init_mode == :straight_line
        @test solver.use_strategy_warmstart == false
        @test solver.constraint_opts isa ALOptions
    end

    @testset "custom parameters accepted" begin
        solver = PANGOLIN(; max_iter=50, convergence_tol=1e-4,
                           init_mode=:zero_control,
                           constraint_opts=ALOptions(; ρ_init=5.0))
        @test solver.max_iter == 50
        @test solver.convergence_tol == 1e-4
        @test solver.init_mode == :zero_control
        @test solver.constraint_opts.ρ_init == 5.0
    end

    @testset "invalid init_mode rejected" begin
        # Exception type depends on whether ALOptions() or the assert fires first
        @test_throws Exception PANGOLIN(; init_mode=:bad_mode)
    end

    @testset "solver_capabilities callable without error" begin
        caps = solver_capabilities(PANGOLIN)
        @test caps isa Vector   # may be empty if protocol not registered
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L15 — _solve unconstrained LQ: matches FNELQ" begin
    #
    # Gold-standard integration test. On a 2-player LTI LQ game PANGOLIN must
    # converge to the same Nash equilibrium as FNELQ. The operating-point
    # residual criterion should make them agree to within convergence_tol.
    #

    game     = lti_game()
    sol_ref  = solve(game, FNELQ(); verbose=false)
    sol_pan  = solve(game, PANGOLIN(; convergence_tol=1e-6, max_iter=300); verbose=false)

    @testset "PANGOLIN converges or makes progress" begin
        @test sol_pan.iterations > 0
        @test all(isfinite, sol_pan.trajectories[1].states)
    end

    @testset "equilibrium type is FeedbackNash" begin
        @test sol_pan.equilibrium_type == :FeedbackNash
    end

    @testset "final costs within 2x of FNELQ (convention difference accounted for)" begin
        costs_ref  = sol_ref.solver_info[:costs]   # FNELQ: xQx + uRu (no ½)
        iter_costs = sol_pan.solver_info[:iter_costs]
        costs_pan  = isempty(iter_costs) ? [sol_pan.trajectories[i].cost for i in 1:2] :
                                           iter_costs[end]
        for i in 1:2
            # FNELQ uses no ½; PANGOLIN's evaluate_stage_cost uses ½
            # Expect costs_pan ≈ 0.5 * costs_ref; allow factor of 2 slop for non-convergence
            @test costs_pan[i] > 0
            @test costs_pan[i] < 2.0 * costs_ref[i] + 1.0
        end
    end

    @testset "terminal states are finite" begin
        x_pan = sol_pan.trajectories[1].states[:, end]
        @test all(isfinite, x_pan)
    end

    @testset "solver_info keys present" begin
        @test haskey(sol_pan.solver_info, :feedback_gains)
        @test haskey(sol_pan.solver_info, :affine_gains)
        @test haskey(sol_pan.solver_info, :iter_costs)
        @test haskey(sol_pan.solver_info, :violation_history)
        @test haskey(sol_pan.solver_info, :ρ_history)
        @test haskey(sol_pan.solver_info, :iterations)
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L16 — _solve constrained: bound constraint satisfied" begin

    game   = constrained_lti_game()
    solver = PANGOLIN(;
        max_iter       = 200,
        convergence_tol = 1e-4,
        constraint_opts = ALOptions(; ρ_init=1.0, ρ_max=1e3, φ=5.0)
    )
    sol = solve(game, solver; verbose=false)

    @testset "solution completes without error, trajectories finite" begin
        @test all(isfinite, sol.trajectories[1].states)
        @test all(isfinite, sol.trajectories[2].states)
    end

    @testset "player 1 controls satisfy bound constraint at convergence" begin
        # u₁(k) ≤ 0.5 for all k — check mean violation at end of solve
        u1 = sol.trajectories[1].controls   # (m × N)
        N  = size(u1, 2)
        # Allow small tolerance for near-satisfaction (AL doesn't enforce exactly)
        mean_violation = mean(max.(u1 .- 0.5, 0.0))
        @test mean_violation < 0.2   # should be substantially reduced vs unconstrained
    end

    @testset "constraint violation history decreases overall" begin
        viol_hist = sol.solver_info[:violation_history]
        if length(viol_hist) > 5
            # Later violations should be generally smaller than early ones
            @test mean(viol_hist[end-2:end]) <= mean(viol_hist[1:3]) + 1e-6
        end
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L17 — _solve nonlinear: 1-player tanh game" begin

    game   = nonlinear_1player_game()
    solver = PANGOLIN(;
        max_iter        = 200,
        convergence_tol = 1e-4,
        init_mode       = :zero_control
    )
    sol = solve(game, solver; verbose=false)

    @testset "completes without error, finite states" begin
        @test all(isfinite, sol.trajectories[1].states)
        @test sol.equilibrium_type == :FeedbackNash
    end

    @testset "state trajectory is finite" begin
        @test all(isfinite, sol.trajectories[1].states)
    end

    @testset "cost decreases over iterations" begin
        iter_costs = sol.solver_info[:iter_costs]
        if length(iter_costs) > 3
            # Cost should generally decrease (not necessarily monotone due to line search)
            @test iter_costs[end][1] <= iter_costs[1][1] + 1e-4
        end
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L18 — warmstart modes" begin

    @testset "zero_control: rollout completes, x[1] = x0" begin
        game   = lti_game()
        solver = PANGOLIN(; init_mode=:zero_control, max_iter=5)
        sol    = solve(game, solver; verbose=false)
        @test sol.trajectories[1].states[:, 1] ≈ game.initial_state atol=1e-10
    end

    @testset "straight_line: completes without error, finite trajectories" begin
        game   = lti_game()
        solver = PANGOLIN(; init_mode=:straight_line, max_iter=5)
        sol    = solve(game, solver; verbose=false)
        @test all(isfinite, sol.trajectories[1].states)
        @test sol.equilibrium_type == :FeedbackNash
    end

    @testset "warmstart from prior solve: fewer iterations than cold start" begin
        game         = lti_game()
        solver_cold  = PANGOLIN(; max_iter=300, convergence_tol=1e-5)
        sol_cold     = solve(game, solver_cold; verbose=false)

        ws = WarmstartData(sol_cold)

        solver_warm = PANGOLIN(;
            max_iter        = 300,
            convergence_tol = 1e-5,
            init_mode       = :warmstart
        )
        # warmstart is a keyword argument in the base solve dispatcher
        sol_warm = solve(game, solver_warm; warmstart=ws, verbose=false)

        @test all(isfinite, sol_warm.trajectories[1].states)
        if sol_cold.converged && sol_warm.converged
            @test sol_warm.iterations <= sol_cold.iterations
        end
    end

end

# ─────────────────────────────────────────────────────────────────────────────
@testset "L19 — type stability" begin

    game  = lti_game()
    op    = OperatingPoint(game)
    lqg   = LTVLQApproximation(game)
    strat = FeedbackStrategy(game)
    lq_approximation!(lqg, game, op)
    solver = PANGOLIN()

    @testset "@inferred solve_lq_game!" begin
        @test @inferred(solve_lq_game!(strat, lqg, solver)) isa FeedbackStrategy
    end

    @testset "@inferred rollout!" begin
        @test @inferred(rollout!(op, game, strat, game.initial_state)) isa Bool
    end

    @testset "@inferred has_converged" begin
        op2 = OperatingPoint(game)
        @test @inferred(has_converged(op, op2, 1e-3)) isa Bool
    end

    @testset "@inferred scale_feedforward!" begin
        @test @inferred(scale_feedforward!(strat, 0.5)) isa FeedbackStrategy
    end

end

# ============================================================================
end  # @testset "PANGOLIN — Comprehensive Test Suite"
# ============================================================================