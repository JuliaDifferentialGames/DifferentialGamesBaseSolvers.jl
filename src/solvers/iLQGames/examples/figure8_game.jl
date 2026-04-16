# ============================================================================
# figure8_game.jl — Two-player shared unicycle game (figure-8 example)
#
# Reproduced from the iLQGames.jl README example:
#   https://github.com/JuliaGameTheoreticPlanning/iLQGames.jl
#
# ── Problem ──────────────────────────────────────────────────────────────────
#
# Two players share control of a single 4D unicycle. The equilibrium produces
# a figure-8 trajectory: Player 1 steers toward the origin while Player 2
# maintains a target speed, and their interaction forces a looping path.
#
# ── Dynamics ─────────────────────────────────────────────────────────────────
#
# State:    x = (px, py, φ, v) ∈ ℝ⁴
#   px, py  : Cartesian position (m)
#   φ       : Heading angle (rad)
#   v       : Speed (m/s)
#
# Controls:
#   Player 1:  u₁ = ω  (steering rate, rad/s)
#   Player 2:  u₂ = a  (longitudinal acceleration, m/s²)
#
# Continuous-time ODE (unicycle):
#   ṗx = v·cos(φ)
#   ṗy = v·sin(φ)
#   φ̇  = ω
#   v̇  = a
#
# ── Objectives ───────────────────────────────────────────────────────────────
#
# Player 1 (steering): Regulate position to origin, penalise steering effort.
#   ℓ₁(x, u₁) = w_pos·(px² + py²) + w_ω·ω²
#
# Player 2 (throttle): Track reference speed v_ref, penalise acceleration.
#   ℓ₂(x, u₂) = w_v·(v - v_ref)² + w_a·a²
#
# Both players have mild quadratic terminal costs on their respective terms.
#
# ── Parameters (matching iLQGames.jl defaults) ───────────────────────────────
#   dt = 0.1 s,  N = 200 steps,  tf = 20.0 s
#   x0 = [1.0, 1.0, 0.0, 0.5]   (1 m off-origin, v₀ = 0.5 m/s)
#   v_ref = 1.0 m/s
#
# ── Reference ────────────────────────────────────────────────────────────────
# Peters, L. & Sunberg, Z. (2020). iLQGames.jl: Rapidly Designing and Solving
#   Differential Games in Julia. arXiv:2002.10185.
# Fridovich-Keil, D. et al. (2020). Efficient Iterative Linear-Quadratic
#   Approximations for Nonlinear Multi-Player General-Sum Differential Games.
#   ICRA 2020. arXiv:1909.04694.
# ============================================================================

using DifferentialGamesBase
using DifferentialGamesBaseSolvers
using LinearAlgebra

# ── Parameters ───────────────────────────────────────────────────────────────

const F8_DT     = 0.1      # time step (s)
const F8_N      = 200      # horizon steps
const F8_TF     = F8_N * F8_DT
const F8_VREF   = 1.0      # Player 2 reference speed (m/s)
const F8_X0     = [1.0, 1.0, 0.0, 0.5]

# ── Unicycle dynamics ─────────────────────────────────────────────────────────
#
# Joint control u = [ω; a] ∈ ℝ².
# CoupledNonlinearDynamics takes f(x, u, p, t) -> ẋ.

function unicycle_ode(x, u, p, t)
    px, py, φ, v = x[1], x[2], x[3], x[4]
    ω, a         = u[1], u[2]
    return [v * cos(φ), v * sin(φ), ω, a]
end

# ── Stage costs ───────────────────────────────────────────────────────────────
#
# Exactly matching iLQGames.jl README:
#   costs = (FunctionPlayerCost((g, x, u, t) -> (x[1]^2 + x[2]^2 + u[1]^2)),
#            FunctionPlayerCost((g, x, u, t) -> ((x[4] - 1)^2 + u[2]^2)))
#
# No terminal costs — the iLQGames.jl example uses stage costs only.
# Using NonlinearTerminalCost(x -> 0.0) gives a zero Hessian, which is
# equivalent to no terminal cost and avoids changing the game structure.

function stage_cost_p1(x, u, p, t)
    return x[1]^2 + x[2]^2 + u[1]^2
end

function stage_cost_p2(x, u, p, t)
    return (x[4] - F8_VREF)^2 + u[1]^2   # u[1] = a, only element in slice
end

terminal_cost_p1(x, p) = zero(eltype(x))
terminal_cost_p2(x, p) = zero(eltype(x))

# ── Game constructor ──────────────────────────────────────────────────────────

"""
    figure8_game(; dt=F8_DT, N=F8_N, x0=F8_X0) -> GameProblem{Float64}

Construct the two-player shared-unicycle figure-8 game.

# State layout
  x = [px, py, φ, v] — shared unicycle state (dim 4)

# Player layout
  Player 1: controls u₁ = ω  (steering rate);  aims to stay near origin
  Player 2: controls u₂ = a  (acceleration);   aims to track speed F8_VREF

# Equilibrium
iLQGames finds a local open-loop Nash equilibrium. The emergent figure-8
trajectory arises because Player 1's origin-seeking creates a centripetal
steering force, while Player 2's speed-tracking creates forward acceleration —
the interaction produces a looping orbit around the origin.

# Reference
Peters & Sunberg (2020), iLQGames.jl README minimal example.
"""
function figure8_game(;
    dt::Float64 = F8_DT,
    N::Int      = F8_N,
    x0::Vector  = F8_X0,
    v_ref::Float64 = F8_VREF
)
    tf = Float64(N * dt)

    # Shared unicycle dynamics — total m=2, Player 1 gets u[1:1], Player 2 gets u[2:2]
    # Per-player breakdown is carried by the Players' .m fields and GameMetadata.
    dyn = CoupledNonlinearDynamics(unicycle_ode, 4, 2)

    sc1 = NonlinearStageCost(stage_cost_p1)
    sc2 = NonlinearStageCost(stage_cost_p2)
    tc1 = NonlinearTerminalCost(terminal_cost_p1)
    tc2 = NonlinearTerminalCost(terminal_cost_p2)

    obj1 = PlayerObjective(1, sc1, tc1)
    obj2 = PlayerObjective(2, sc2, tc2)

    # Dynamics placeholder — not used; shared CoupledNonlinearDynamics is the actual ODE.
    _noop = (x, u, p, t) -> zeros(eltype(x), length(x))

    p1 = Player{Float64}(1, 4, 1, x0, _noop, obj1, AbstractPrivateConstraint[])
    p2 = Player{Float64}(2, 4, 1, x0, _noop, obj2, AbstractPrivateConstraint[])

    return DifferentialGame(dyn, [p1, p2], tf, dt)
end

# ── Solver calls ─────────────────────────────────────────────────────────────

"""
    solve_figure8(; verbose=false, solver=nothing, kwargs...) -> GNEPSolution{Float64}

Solve the figure-8 game from zero-control initialisation, matching the C++
ilqgames default: `P⁰ᵢ = 0`, `α⁰ᵢ = 0`, straight-line rollout from x₀.

Keyword arguments other than `verbose` and `solver` are forwarded to
`figure8_game` (e.g. `dt`, `N`, `x0`).

# Solver defaults
  max_iter = 400, ε_conv = 0.05 (trajectory ∞-norm), μ_init = 0.01

# Warmstart
To initialise from a prior solution, use `solve_figure8_warmstart`.
"""
function solve_figure8(;
    verbose::Bool   = false,
    solver::Union{Nothing, iLQGames} = nothing,
    kwargs...
)
    game    = figure8_game(; kwargs...)
    slv     = solver !== nothing ? solver :
                iLQGames(max_iter=400, ε_conv=0.05, μ_init=0.0, μ_max=1e4, μ_scale=10.0, μ_decay=0.5)
    return solve(game, slv; verbose, check_compatibility=false)
end

"""
    figure8_warmstart(game) -> (X_warm, U_warm)

Compute a sinusoidal reference trajectory for the figure-8 game that
lies in the basin of attraction of the figure-8 Nash equilibrium.

The warmstart prescribes:
  ω(t) = (π/2)·(2π/T)·cos(2πt/T)   — heading traces φ = (π/2)·sin(2πt/T)
  a(t) = (v_ref - v₀)/3  for t ≤ 3s, 0 otherwise  — speed ramp

With T = 10 s over a 20-second horizon, the unicycle makes two full loops
of radius ≈ v_ref/(2π/T) ≈ 1.6 m centred near the starting point,
closely approximating the figure-8 Nash trajectory.

Returns `(X_warm, U_warm)` where `X_warm` is the rollout of `U_warm` from
`game.initial_state`.
"""
function figure8_warmstart(game::GameProblem)
    N    = n_steps(game)
    dt   = game.time_horizon.dt
    x0   = game.initial_state
    v0   = x0[4]
    T    = 10.0   # steering period (s) — two full cycles over 20 s

    U = zeros(2, N)
    for k in 1:N
        t       = (k-1) * dt
        U[1, k] = (π/2) * (2π/T) * cos(2π * t / T)
        U[2, k] = k <= round(Int, 3.0/dt) ? (F8_VREF - v0) / 3.0 : 0.0
    end

    t_vec = collect(range(0.0, game.time_horizon.tf, length=N+1))
    X     = rollout(game.dynamics, x0, U, nothing, t_vec)
    return X, U
end

"""
    solve_figure8_warmstart(; verbose=false, solver=nothing, kwargs...) -> GNEPSolution{Float64}

Solve the figure-8 game from a sinusoidal warmstart trajectory computed by
`figure8_warmstart`. Use this when zero-control initialisation converges to
the degenerate zero-velocity Nash rather than the figure-8 equilibrium.

The warmstart is constructed analytically: heading traces φ(t) = (π/2)·sin(2πt/10),
which prescribes the exact steering rate and produces figure-8-shaped loops when
combined with a brief acceleration ramp from v₀ to v_ref.

All keyword arguments are forwarded to `figure8_game`.
"""
function solve_figure8_warmstart(;
    verbose::Bool   = false,
    solver::Union{Nothing, iLQGames} = nothing,
    kwargs...
)
    game    = figure8_game(; kwargs...)
    slv     = solver !== nothing ? solver :
              iLQGames(max_iter=400, ε_conv=0.05, μ_init=0.01, μ_max=1e3)

    X_warm, U_warm = figure8_warmstart(game)

    N     = n_steps(game)
    cd    = game.metadata.control_dims
    co    = [0; cumsum(cd)[1:end-1]]
    t_vec = collect(range(0.0, game.time_horizon.tf, length=N+1))

    # Package as GNEPSolution for the warmstart interface.
    # Only the trajectory (X, U) matters — costs and convergence flag are placeholders.
    trajs = [
        Trajectory(i, X_warm, U_warm[co[i]+1:co[i]+cd[i], :], t_vec, 0.0)
        for i in 1:2
    ]
    strat_warm = OpenLoopStrategy(
        [U_warm[co[i]+1:co[i]+cd[i], :] for i in 1:2], cd, t_vec
    )
    ws = GNEPSolution(game, trajs;
                      state_trajectory = X_warm,
                      strategy         = strat_warm,
                      equilibrium_type = :OpenLoopNash,
                      converged        = false,
                      iterations       = 0,
                      solve_time       = 0.0)

    return solve(game, slv; warmstart=ws, verbose, check_compatibility=false)
end

using Plots
sol  = solve_figure8(x0=[1.0, 1.0, 0.0, 0.5], verbose=true)
anim = animate_solution(sol; fps=20, subsample=2)
save_animation(anim, "figure8.gif")
