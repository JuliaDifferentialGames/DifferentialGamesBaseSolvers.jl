# ============================================================================
# test_algames.jl
#
# Test suite for ALGAMES (algames.jl + utils.jl).
# Run via:  include("test_algames.jl")  from a session that has DGB on LOAD_PATH.
#
# Categories (per CLAUDE.md guidelines):
#   Cat 0 — Unit:        workspace allocation, packing, index helpers
#   Cat 1 — Ground truth: single-player LQR, 2-player LQ game
#   Cat 2 — Optimality:  stationarity, dynamics, dual feasibility
#   Cat 3 — Cross-check: gradient finite-difference test at solution
#   Cat 4 — Convergence: converged flag, warm-start, constraint satisfaction
# ============================================================================

using Test
using LinearAlgebra
using ForwardDiff
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# ============================================================================
# Problem factories
# ============================================================================

"""
    make_double_integrator_game(; N, dt, np, Q_scale, R_scale, Qf_scale)

N-step, np-player shared double-integrator game.
State: x = [pos; vel], n=2.  Each player controls a scalar acceleration.
Dynamics (exact discrete): A=[1 dt; 0 1], Bᵢ = [dt²/2; dt] * eᵢ (unit vector).
"""
function make_double_integrator_game(;
    N       ::Int     = 10,
    dt      ::Float64 = 0.1,
    np      ::Int     = 2,
    Q_scale ::Float64 = 1.0,
    R_scale ::Float64 = 1.0,
    Qf_scale::Float64 = 5.0,
)
    T  = Float64
    tf = dt * N
    n  = 2
    A  = [1.0 dt; 0.0 1.0]
    # Each player gets the same B column (scalar control, affects velocity only).
    # This makes the symmetric game have u¹=u² at equilibrium.
    B  = [reshape([0.5*dt^2; dt], n, 1) for _ in 1:np]
    Q  = [Q_scale  * Matrix{T}(I, n, n) for _ in 1:np]
    R  = [R_scale  * Matrix{T}(I, 1, 1) for _ in 1:np]
    Qf = [Qf_scale * Matrix{T}(I, n, n) for _ in 1:np]
    x0 = [1.0, 0.0]
    return LQGameProblem(A, B, Q, R, Qf, x0, T(tf); dt=T(dt))
end

"""Evaluate the total trajectory cost for player i given X (n×N+1) and Ui (mᵢ×N)."""
function traj_cost(game::GameProblem{T}, i::Int,
                   X::Matrix{T}, Ui::Matrix{T}) where {T}
    obj = get_objective(game, i)
    N   = size(Ui, 2)
    s   = sum(evaluate_stage_cost(obj.stage_cost, X[:,k], Ui[:,k], nothing, k)
              for k in 1:N)
    s  += evaluate_terminal_cost(obj.terminal_cost, X[:,N+1], nothing)
    return s
end

# ============================================================================
# Cat 0 — Unit tests: workspace, packing, index helpers
# ============================================================================

@testset "Cat 0: workspace allocation dimensions" begin
    game   = make_double_integrator_game(N=8, np=2)
    solver = ALGAMES()
    ws     = ALGAMESWorkspace(game, solver, nothing)

    @test ws.n  == 2
    @test ws.np == 2
    @test ws.N  == 8
    @test size(ws.X)        == (2, 9)     # n × N+1
    @test size(ws.U[1])     == (1, 8)     # m₁ × N
    @test size(ws.U[2])     == (1, 8)
    @test size(ws.Λ_dyn[1]) == (2, 8)    # n × N
    @test size(ws.Λ_dyn[2]) == (2, 8)
    @test ws.X[:, 1] ≈ game.initial_state
    # Unconstrained → no constraint duals
    @test length(ws.λ) == 0
    @test length(ws.ρ) == 0
end

@testset "Cat 0: _y_dim equals _G_dim" begin
    game = make_double_integrator_game(N=5, np=3)
    ws   = ALGAMESWorkspace(game, ALGAMES(), nothing)
    @test _y_dim(ws) == _G_dim(ws)
end

@testset "Cat 0: y-vector packing round-trips" begin
    game   = make_double_integrator_game(N=5, np=2)
    ws     = ALGAMESWorkspace(game, ALGAMES(), nothing)

    # Randomise all workspace fields
    ws.X[:, 2:end] .= randn(ws.n, ws.N)
    for i in 1:ws.np
        ws.U[i]      .= randn(size(ws.U[i]))
        ws.Λ_dyn[i]  .= randn(ws.n, ws.N)
    end

    y   = _pack_y(ws)
    ws2 = _unpack_y(ws, y)

    @test ws2.X[:, 2:end] ≈ ws.X[:, 2:end]
    @test ws2.X[:, 1]     ≈ ws.X[:, 1]           # x₀ unchanged
    for i in 1:ws.np
        @test ws2.U[i]     ≈ ws.U[i]
        @test ws2.Λ_dyn[i] ≈ ws.Λ_dyn[i]
    end
end

@testset "Cat 0: _apply_step! correctness" begin
    game   = make_double_integrator_game(N=6, np=2)
    ws     = ALGAMESWorkspace(game, ALGAMES(), nothing)
    n_y    = _y_dim(ws)
    δy     = randn(n_y)
    α      = 0.7

    X_pre  = copy(ws.X)
    U_pre  = [copy(ws.U[i]) for i in 1:ws.np]
    Λ_pre  = [copy(ws.Λ_dyn[i]) for i in 1:ws.np]

    _apply_step!(ws, δy, α)

    # x₀ untouched
    @test ws.X[:, 1] ≈ X_pre[:, 1]
    # x(k+1) updated correctly
    for k in 1:ws.N
        @test ws.X[:, k+1] ≈ X_pre[:, k+1] .+ α .* δy[_x_idx(ws, k)]
    end
    # u updated
    for i in 1:ws.np, k in 1:ws.N
        @test ws.U[i][:, k] ≈ U_pre[i][:, k] .+ α .* δy[_u_idx(ws, i, k)]
    end
    # λ_dyn updated
    for i in 1:ws.np, k in 1:ws.N
        @test ws.Λ_dyn[i][:, k] ≈ Λ_pre[i][:, k] .+ α .* δy[_λdyn_idx(ws, i, k)]
    end
end

@testset "Cat 0: G-index helpers span correct range" begin
    game = make_double_integrator_game(N=4, np=2)
    ws   = ALGAMESWorkspace(game, ALGAMES(), nothing)
    n_G  = _G_dim(ws)

    # Every row of G must be covered by exactly one index helper
    covered = zeros(Int, n_G)
    for i in 1:ws.np, k in 1:ws.N
        for idx in _Gx_idx(ws, i, k);  covered[idx] += 1;  end
        for idx in _Gu_idx(ws, i, k);  covered[idx] += 1;  end
    end
    for k in 1:ws.N
        for idx in _Gdyn_idx(ws, k);   covered[idx] += 1;  end
    end
    @test all(covered .== 1)   # each row covered exactly once
end

@testset "Cat 0: residual Jacobian = ForwardDiff Jacobian" begin
    # H = ∂G/∂y computed via _build_jacobian! must equal ForwardDiff.jacobian
    # through _G_from_y. They use the same code path so this is a consistency check,
    # but it also catches any bug where buf.H is not filled correctly.
    game   = make_double_integrator_game(N=3, np=2)
    solver = ALGAMES()
    ws     = ALGAMESWorkspace(game, solver, nothing)
    n_y    = _y_dim(ws)
    buf    = ALGAMESBuffers(Float64, n_y)

    _build_residual!(buf, game, ws)
    _build_jacobian!(buf, game, ws)

    y0   = _pack_y(ws)
    H_fd = ForwardDiff.jacobian(y -> _G_from_y(game, ws, y), y0)
    @test buf.H ≈ H_fd atol=1e-10
end

@testset "Cat 0: finite-difference Jacobian cross-check" begin
    game   = make_double_integrator_game(N=3, np=2)
    ws     = ALGAMESWorkspace(game, ALGAMES(), nothing)
    n_y    = _y_dim(ws)
    buf    = ALGAMESBuffers(Float64, n_y)
    _build_jacobian!(buf, game, ws)

    y0   = _pack_y(ws)
    ε    = 1e-6
    H_findiff = zeros(n_y, n_y)
    for j in 1:n_y
        yp = copy(y0); yp[j] += ε
        ym = copy(y0); ym[j] -= ε
        H_findiff[:, j] .= (_G_from_y(game, ws, yp) .- _G_from_y(game, ws, ym)) ./ (2ε)
    end
    @test buf.H ≈ H_findiff atol=1e-5
end

# ============================================================================
# Cat 1 — Analytical ground truth
# ============================================================================

@testset "Cat 1: 1-player LQR dynamics feasibility at solution" begin
    # Single-player double integrator.  ALGAMES must find dynamics-feasible X.
    game   = make_double_integrator_game(N=10, np=1)
    solver = ALGAMES(tol_opt=1e-8, tol_dyn=1e-8, outer_iter=100, inner_iter=20)
    sol    = solve(game, solver)

    @test sol.converged

    X   = sol.state_trajectory
    U   = get_trajectory(sol, 1).controls
    dyn = game.dynamics
    N   = n_steps(game)
    for k in 1:N
        x_next = evaluate_dynamics(dyn, X[:,k], U[:,k], nothing, k)
        @test X[:, k+1] ≈ x_next atol=1e-7
    end
    @test sol.solver_info[:opt_vio] < 1e-6
    @test sol.solver_info[:dyn_vio] < 1e-6
end

@testset "Cat 1: 1-player cost is non-negative (Q, R, Qf PSD/PD)" begin
    game = make_double_integrator_game(N=10, np=1)
    sol  = solve(game, ALGAMES(tol_opt=1e-8, tol_dyn=1e-8))
    @test get_cost(sol, 1) >= 0.0
end

@testset "Cat 1: 2-player symmetric game — costs equal at equilibrium" begin
    # Symmetric game (identical B, Q, R, Qf, x₀ shared):
    # by symmetry u¹=u² and therefore cost¹=cost².
    game   = make_double_integrator_game(N=8, np=2)
    solver = ALGAMES(tol_opt=1e-7, tol_dyn=1e-7, outer_iter=100, inner_iter=20)
    sol    = solve(game, solver)

    @test sol.converged

    c1 = get_cost(sol, 1)
    c2 = get_cost(sol, 2)
    @test abs(c1 - c2) / max(abs(c1), 1e-10) < 1e-4

    # Controls must be identical by symmetry
    U1 = get_trajectory(sol, 1).controls
    U2 = get_trajectory(sol, 2).controls
    @test U1 ≈ U2 atol=1e-5
end

@testset "Cat 1: dynamics satisfied (2-player)" begin
    game   = make_double_integrator_game(N=10, np=2)
    solver = ALGAMES(tol_dyn=1e-8)
    sol    = solve(game, solver)

    X      = sol.state_trajectory
    U1     = get_trajectory(sol, 1).controls
    U2     = get_trajectory(sol, 2).controls
    U_flat = vcat(U1, U2)
    dyn    = game.dynamics
    N      = n_steps(game)
    for k in 1:N
        x_next = evaluate_dynamics(dyn, X[:,k], U_flat[:,k], nothing, k)
        @test X[:, k+1] ≈ x_next atol=1e-6
    end
end

# ============================================================================
# Cat 2 — Optimality conditions at solution
# ============================================================================

@testset "Cat 2: stationarity and dynamics residuals below tolerance" begin
    game   = make_double_integrator_game(N=8, np=2)
    solver = ALGAMES(tol_opt=1e-6, tol_dyn=1e-6)
    sol    = solve(game, solver)

    @test sol.solver_info[:opt_vio] < solver.tol_opt * 10
    @test sol.solver_info[:dyn_vio] < solver.tol_dyn * 10
end

@testset "Cat 2: dual feasibility (λ ≥ 0) for inequality constraints" begin
    T  = Float64
    tf = T(1.0); dt = T(0.1); N = Int(round(tf/dt))
    f  = (x, u, p, t) -> [x[2]; u[1]]
    Q  = Matrix{T}(I, 2, 2)
    Qf = 5.0 * Matrix{T}(I, 2, 2)
    R  = Matrix{T}(I, 1, 1)
    o1 = PlayerObjective(1, LQStageCost(Q, R), LQTerminalCost(Qf))
    o2 = PlayerObjective(2, LQStageCost(Q, R), LQTerminalCost(Qf))

    # Collision avoidance: players must stay ≥ 0.5 apart in x[1] dimension
    prx = collision_avoidance([1,2]; i_offset=0, j_offset=2, pos_dim=1, d_min=0.5)
    p1  = Player{T}(1, 2, 1, [2.0, 0.0], f, o1, [])
    p2  = Player{T}(2, 2, 1, [-2.0, 0.0], f, o2, [])

    game   = DifferentialGame([p1, p2], tf, dt; shared_constraints=[prx])
    solver = ALGAMES(tol_con=1e-3, outer_iter=60)
    sol    = solve(game, solver)

    λ = sol.solver_info[:dual_variables][:λ]
    @test all(λ .>= -1e-8)    # numerical noise tolerance
    @test sol.solver_info[:con_vio] < solver.tol_con * 10
end

@testset "Cat 2: complementarity holds at solution" begin
    # At convergence, inactive constraints should have λ ≈ 0.
    T  = Float64
    tf = T(1.0); dt = T(0.1); N = Int(round(tf/dt))
    f  = (x, u, p, t) -> [x[2]; u[1]]
    Q  = Matrix{T}(I, 2, 2)
    Qf = 5.0 * Matrix{T}(I, 2, 2)
    R  = Matrix{T}(I, 1, 1)
    o1 = PlayerObjective(1, LQStageCost(Q, R), LQTerminalCost(Qf))
    o2 = PlayerObjective(2, LQStageCost(Q, R), LQTerminalCost(Qf))

    # Large d_min so agents start far apart → constraint inactive at solution
    prx = collision_avoidance([1,2]; i_offset=0, j_offset=2, pos_dim=1, d_min=0.1)
    p1  = Player{T}(1, 2, 1, [5.0, 0.0], f, o1, [])
    p2  = Player{T}(2, 2, 1, [-5.0, 0.0], f, o2, [])

    game   = DifferentialGame([p1, p2], tf, dt; shared_constraints=[prx])
    solver = ALGAMES(tol_opt=1e-5, tol_dyn=1e-5, tol_con=1e-4)
    sol    = solve(game, solver)

    # Agents start and stay far apart → constraint inactive → λ ≈ 0
    λ = sol.solver_info[:dual_variables][:λ]
    @test maximum(abs.(λ)) < 1e-3   # λ stays near zero when inactive
end

# ============================================================================
# Cat 3 — Cross-check: finite-difference gradient at solution
# ============================================================================

@testset "Cat 3: player 1 cost gradient ≈ 0 at solution (FD check)" begin
    # Perturb player 1's controls and verify cost changes at 2nd order only.
    game   = make_double_integrator_game(N=6, np=2)
    solver = ALGAMES(tol_opt=1e-8, tol_dyn=1e-8, inner_iter=20, outer_iter=100)
    sol    = solve(game, solver)

    @test sol.converged

    X    = sol.state_trajectory
    U1   = copy(get_trajectory(sol, 1).controls)
    U2   = copy(get_trajectory(sol, 2).controls)
    dyn  = game.dynamics
    N    = n_steps(game)

    # Finite-difference gradient ∂J¹/∂u¹ₖ at the solution
    ε   = 1e-5
    max_grad = 0.0
    for k in 1:N
        # Perturb u¹(k) in place, recompute trajectory from k onwards
        U1p = copy(U1); U1p[1, k] += ε
        U1m = copy(U1); U1m[1, k] -= ε
        U_p = vcat(U1p, U2)
        U_m = vcat(U1m, U2)
        # Roll out from current X[:,k] with the perturbed control
        Xp = copy(X); Xm = copy(X)
        for kk in k:N
            Xp[:, kk+1] .= evaluate_dynamics(dyn, Xp[:,kk], U_p[:,kk], nothing, kk)
            Xm[:, kk+1] .= evaluate_dynamics(dyn, Xm[:,kk], U_m[:,kk], nothing, kk)
        end
        cp = traj_cost(game, 1, Xp, U1p)
        cm = traj_cost(game, 1, Xm, U1m)
        max_grad = max(max_grad, abs((cp - cm) / (2ε)))
    end
    @test max_grad < 1e-4
end

# ============================================================================
# Cat 4 — Convergence
# ============================================================================

@testset "Cat 4: solver_capabilities contains expected symbols" begin
    caps = solver_capabilities(ALGAMES)
    @test :LQGame          in caps
    @test :NonlinearGame   in caps
    @test :GNEP            in caps
    @test :ConstrainedGame in caps
    @test :DiscreteTime    in caps
end

@testset "Cat 4: converged flag and iteration count" begin
    game   = make_double_integrator_game(N=8, np=2)
    solver = ALGAMES()
    sol    = solve(game, solver)

    @test sol.converged
    @test 1 <= sol.iterations <= solver.outer_iter
    @test sol.solve_time > 0.0
end

@testset "Cat 4: warm-start — same solution as cold start" begin
    game   = make_double_integrator_game(N=8, np=2)
    solver = ALGAMES(tol_opt=1e-7, tol_dyn=1e-7)

    sol_cold = solve(game, solver)
    sol_warm = solve(game, solver; warmstart=sol_cold)

    @test sol_warm.converged
    # Warm-started solve should not require more outer iterations
    @test sol_warm.iterations <= sol_cold.iterations + 2
    # Solutions should agree in cost
    @test abs(get_cost(sol_warm, 1) - get_cost(sol_cold, 1)) < 1e-5
    @test abs(get_cost(sol_warm, 2) - get_cost(sol_cold, 2)) < 1e-5
end

@testset "Cat 4: constraint violation satisfied at convergence" begin
    T  = Float64
    tf = T(1.0); dt = T(0.1); N = Int(round(tf/dt))
    f  = (x, u, p, t) -> [x[2]; u[1]]
    Q  = Matrix{T}(I, 2, 2)
    Qf = 5.0 * Matrix{T}(I, 2, 2)
    R  = Matrix{T}(I, 1, 1)
    o1 = PlayerObjective(1, LQStageCost(Q, R), LQTerminalCost(Qf))
    o2 = PlayerObjective(2, LQStageCost(Q, R), LQTerminalCost(Qf))
    prx = collision_avoidance([1,2]; i_offset=0, j_offset=2, pos_dim=1, d_min=0.3)
    p1  = Player{T}(1, 2, 1, [1.5, 0.0], f, o1, [])
    p2  = Player{T}(2, 2, 1, [-1.5, 0.0], f, o2, [])
    game   = DifferentialGame([p1, p2], tf, dt; shared_constraints=[prx])
    solver = ALGAMES(tol_con=1e-3)
    sol    = solve(game, solver)
    @test sol.solver_info[:con_vio] < solver.tol_con * 10
end

@testset "Cat 4: verbose mode runs without error" begin
    game   = make_double_integrator_game(N=3, np=1)
    solver = ALGAMES(outer_iter=2, inner_iter=2)
    @test_nowarn solve(game, solver; verbose=true)
end

@testset "Cat 4: GNEPSolution has OpenLoopStrategy" begin
    game = make_double_integrator_game(N=5, np=2)
    sol  = solve(game, ALGAMES())
    @test has_strategy(sol)
    @test get_strategy(sol) isa OpenLoopStrategy
    @test sol.equilibrium_type == :OpenLoopNash
end

println("\n✓ ALGAMES test suite complete.")