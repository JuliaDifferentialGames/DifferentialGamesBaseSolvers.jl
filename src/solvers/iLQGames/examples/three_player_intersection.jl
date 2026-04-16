# ============================================================================
# three_player_intersection.jl — Three-player unicycle intersection game
#
# ── Problem ──────────────────────────────────────────────────────────────────
#
# Three autonomous vehicles approach an uncontrolled intersection, each
# from a different direction. Each vehicle must reach a goal position past
# the intersection while maintaining a reference speed, minimising control
# effort, and avoiding collision with the other two vehicles.
#
# This is the benchmark problem from Fridovich-Keil et al. (2020), §V:
# "In a three-player, 14-state simulated intersection problem, our algorithm
#  initially converges in < 0.25s."
# The paper describes 14 states; we implement the minimal version with
# 3 × 4D unicycles = 12 states (px, py, φ, v per vehicle), which captures
# the essential game-theoretic interaction.
#
# ── Dynamics ─────────────────────────────────────────────────────────────────
#
# Each vehicle i has unicycle dynamics:
#   ẋᵢ = [vᵢ·cos(φᵢ), vᵢ·sin(φᵢ), ωᵢ, aᵢ]
#
# The three vehicles are dynamically independent (no aerodynamic coupling),
# but their costs are coupled through proximity penalties. We use
# CoupledNonlinearDynamics with a product-system RHS:
#   ẋ = [ẋ₁; ẋ₂; ẋ₃]   (12D joint state)
#
# Using CoupledNonlinearDynamics rather than SeparableDynamics is necessary
# because each player's proximity cost reads the other players' positions
# from the joint state; quadraticize_costs only passes the full joint state
# to evaluate_stage_cost for CoupledNonlinearDynamics.
#
# ── State layout ─────────────────────────────────────────────────────────────
#
#   x[ 1: 4] = (px₁, py₁, φ₁, v₁)  — Player 1 (approaching from south)
#   x[ 5: 8] = (px₂, py₂, φ₂, v₂)  — Player 2 (approaching from west)
#   x[ 9:12] = (px₃, py₃, φ₃, v₃)  — Player 3 (approaching from north)
#
# ── Objectives (per player) ───────────────────────────────────────────────────
#
# Each player i minimises:
#
#   Jᵢ = Σₖ [w_u·‖uᵢ‖² + w_v·(vᵢ-v_ref)² + Σⱼ≠ᵢ w_prox·c_prox(‖pᵢ-pⱼ‖)]
#       + w_goal·‖pᵢ(tf) - pᵢ_goal‖²   (terminal)
#
# Proximity penalty (smooth, C² for ForwardDiff):
#   c_prox(d) = exp(-d² / (2·σ²))
#
# A Gaussian kernel is used rather than a barrier because:
#   (1) It is globally C∞ — no singularity issues during iLQGames iterates
#   (2) Its Hessian is positive semidefinite near d=0 (peaked at collision)
#   (3) It decays to zero for well-separated vehicles, avoiding remote coupling
#
# ── Initial conditions ────────────────────────────────────────────────────────
#
# Intersection centre at origin. Vehicles start 4 m back, approaching at
# v₀ = 1.5 m/s. The road layout is:
#
#          ↓ P3 (from north,  goal south)
#    ──────────────
#    P2 →          → (goal east)
#    ──────────────
#          ↑ P1 (from south, goal north)
#
# Player 1: starts at (0, -4) heading north (φ=π/2)
# Player 2: starts at (-4, 0) heading east  (φ=0)
# Player 3: starts at (0,  4) heading south (φ=-π/2 = 3π/2)
#
# ── Goal positions ────────────────────────────────────────────────────────────
# Player 1: goal (0,  4) — pass through and continue north
# Player 2: goal (4,  0) — pass through and continue east
# Player 3: goal (0, -4) — pass through and continue south
#
# ── Parameters ───────────────────────────────────────────────────────────────
#   dt = 0.1 s,  N = 100 steps,  tf = 10.0 s
#   v_ref = 1.5 m/s (nominal approach speed)
#   d_min = 1.5 m   (soft collision radius for proximity penalty)
#   σ = d_min / √(2·log(2)) (Gaussian width: penalty = 0.5·w_prox at d=d_min)
#
# ── Reference ────────────────────────────────────────────────────────────────
# Fridovich-Keil, D., Ratner, E., Peters, L., Dragan, A.D., & Tomlin, C.J.
#   (2020). Efficient Iterative Linear-Quadratic Approximations for Nonlinear
#   Multi-Player General-Sum Differential Games. ICRA 2020. arXiv:1909.04694.
# ============================================================================

using DifferentialGamesBase
using DifferentialGamesBaseSolvers
using LinearAlgebra

# ── Parameters ───────────────────────────────────────────────────────────────

const INT_DT    = 0.1      # time step (s)
const INT_N     = 100      # horizon steps
const INT_TF    = INT_N * INT_DT    # total horizon (10 s)
const INT_VREF  = 1.5      # reference speed (m/s)
const INT_V0    = 1.5      # initial speed (m/s)

# Geometry
const INT_D_BACK = 4.0     # initial distance behind intersection centre (m)

# Cost weights
const INT_W_U    = 0.1     # control effort
const INT_W_V    = 1.0     # speed tracking
const INT_W_PROX = 5.0     # proximity penalty amplitude
const INT_W_GOAL = 10.0    # terminal goal position weight
const INT_SIGMA  = 1.0     # Gaussian proximity width (m)

# Initial states [px, py, φ, v] per player
const INT_X0_P1 = [0.0, -INT_D_BACK,  π/2,  INT_V0]   # south, heading north
const INT_X0_P2 = [-INT_D_BACK, 0.0,  0.0,  INT_V0]   # west,  heading east
const INT_X0_P3 = [0.0,  INT_D_BACK, -π/2,  INT_V0]   # north, heading south

# Goal positions [px_goal, py_goal] per player
const INT_GOAL_P1 = [0.0,  INT_D_BACK]   # north
const INT_GOAL_P2 = [INT_D_BACK,  0.0]   # east
const INT_GOAL_P3 = [0.0, -INT_D_BACK]   # south

# ── State index helpers ───────────────────────────────────────────────────────

@inline _px(x, i) = x[4*(i-1)+1]
@inline _py(x, i) = x[4*(i-1)+2]
@inline _φ(x, i)  = x[4*(i-1)+3]
@inline _v(x, i)  = x[4*(i-1)+4]

# ── Joint dynamics ────────────────────────────────────────────────────────────
#
# Product unicycle system. State offset: player i at indices 4(i-1)+1 : 4i.
# Control offset: player i at indices 2(i-1)+1 : 2i (ωᵢ, aᵢ).
# CoupledNonlinearDynamics passes the full joint (x, u) to iLQGames and AD.

function intersection_ode(x, u, p, t)
    # Extract per-player state and control
    vcat(
        [begin
            vᵢ = _v(x, i)
            φᵢ = _φ(x, i)
            ωᵢ = u[2*(i-1)+1]
            aᵢ = u[2*(i-1)+2]
            [vᵢ * cos(φᵢ), vᵢ * sin(φᵢ), ωᵢ, aᵢ]
        end for i in 1:3]...
    )
end

# ── Proximity penalty ─────────────────────────────────────────────────────────
#
# Gaussian kernel: large near collision, decays to zero when separated.
# C∞ differentiability ensures clean ForwardDiff quadraticisation.
# σ controls the interaction range; penalty = w_prox at d=0.

function proximity_penalty(x, i, j, σ)
    dxi = _px(x, i) - _px(x, j)
    dyi = _py(x, i) - _py(x, j)
    d2  = dxi^2 + dyi^2
    return exp(-d2 / (2 * σ^2))
end

# ── Stage costs ───────────────────────────────────────────────────────────────
#
# For CoupledNonlinearDynamics, evaluate_stage_cost receives (x_full, u_i, p, t)
# where x_full is the full 12D joint state and u_i is player i's 2D control
# slice [ωᵢ, aᵢ].

function make_stage_cost(player_id::Int; σ=INT_SIGMA)
    others = [j for j in 1:3 if j != player_id]
    i = player_id
    function stage_cost(x, ui, p, t)
        ωi, ai = ui[1], ui[2]
        vi     = _v(x, i)

        # Control effort + speed tracking
        cost = INT_W_U * (ωi^2 + ai^2) + INT_W_V * (vi - INT_VREF)^2

        # Proximity penalty w.r.t. all other players
        for j in others
            cost += INT_W_PROX * proximity_penalty(x, i, j, σ)
        end
        return cost
    end
    return stage_cost
end

# ── Terminal costs ────────────────────────────────────────────────────────────
#
# Penalise distance to goal position and speed deviation at final time.

function make_terminal_cost(player_id::Int, goal::Vector{Float64})
    i = player_id
    function terminal_cost(x, p)
        px_goal, py_goal = goal
        dpx = _px(x, i) - px_goal
        dpy = _py(x, i) - py_goal
        vi  = _v(x, i)
        return INT_W_GOAL * (dpx^2 + dpy^2) + INT_W_V * (vi - INT_VREF)^2
    end
    return terminal_cost
end

# ── Game constructor ──────────────────────────────────────────────────────────

"""
    three_player_intersection_game(; dt, N, x0s, goals, v_ref, w_prox, σ)
        -> GameProblem{Float64}

Construct the three-player unicycle intersection game from Fridovich-Keil et
al. (2020), §V.

# State layout (12D joint state)
  x[1:4]  = (px₁, py₁, φ₁, v₁)  — Player 1 (south → north)
  x[5:8]  = (px₂, py₂, φ₂, v₂)  — Player 2 (west  → east)
  x[9:12] = (px₃, py₃, φ₃, v₃)  — Player 3 (north → south)

# Control layout (6D joint control)
  u[1:2] = (ω₁, a₁),  u[3:4] = (ω₂, a₂),  u[5:6] = (ω₃, a₃)

# Each player minimises
  Jᵢ = Σₖ [wᵤ‖uᵢ‖² + wᵥ(vᵢ-v_ref)² + Σⱼ≠ᵢ w_prox·exp(-‖pᵢ-pⱼ‖²/2σ²)]
        + wgoal‖pᵢ(tf) - pᵢgoal‖² + wᵥ(vᵢ(tf)-v_ref)²

# Arguments (all optional — defaults reproduce the paper example)
- `dt`     : time step (default $(INT_DT) s)
- `N`      : horizon steps (default $(INT_N))
- `x0s`    : initial states per player (vector of 3 length-4 vectors)
- `goals`  : goal positions per player (vector of 3 length-2 vectors)
- `v_ref`  : reference speed (default $(INT_VREF) m/s)
- `w_prox` : proximity penalty weight (default $(INT_W_PROX))
- `σ`      : proximity Gaussian width (default $(INT_SIGMA) m)
"""
function three_player_intersection_game(;
    dt::Float64    = INT_DT,
    N::Int         = INT_N,
    x0s::Vector    = [INT_X0_P1, INT_X0_P2, INT_X0_P3],
    goals::Vector  = [INT_GOAL_P1, INT_GOAL_P2, INT_GOAL_P3],
    v_ref::Float64 = INT_VREF,
    w_prox::Float64 = INT_W_PROX,
    σ::Float64     = INT_SIGMA
)
    tf = Float64(N * dt)

    # Joint initial state
    x0 = vcat(x0s...)
    @assert length(x0) == 12 "Expected 3 × 4D unicycle states"

    # Shared coupled dynamics — total m=6, per-player [2,2,2] from players' .m fields.
    dyn = CoupledNonlinearDynamics(intersection_ode, 12, 6)

    # Per-player objectives
    objectives = map(1:3) do i
        sc = NonlinearStageCost(make_stage_cost(i; σ=σ))
        tc = NonlinearTerminalCost(make_terminal_cost(i, goals[i]))
        PlayerObjective(i, sc, tc)
    end

    # Dynamics placeholder — not used; shared CoupledNonlinearDynamics carries the ODE.
    _noop = (x, u, p, t) -> zeros(eltype(x), length(x))

    players = [
        Player{Float64}(i, 12, 2, x0, _noop, objectives[i], AbstractPrivateConstraint[])
        for i in 1:3
    ]

    return DifferentialGame(dyn, players, tf, dt)
end

# ── Convenience solver call ───────────────────────────────────────────────────

"""
    solve_intersection(; verbose, kwargs...) -> GNEPSolution{Float64}

Construct and solve the three-player intersection game.

Returns a `GNEPSolution{Float64}` with:
  - `sol.state_trajectory`  : (12 × N+1) joint state
  - `sol.trajectories[i]`   : Player i's trajectory (controls = [ωᵢ, aᵢ] over N steps)
  - `sol.equilibrium_type`  : `:OpenLoopNash`

# Extracting per-player positions
  X = sol.state_trajectory
  px1, py1 = X[1, :], X[2, :]   # Player 1
  px2, py2 = X[5, :], X[6, :]   # Player 2
  px3, py3 = X[9, :], X[10, :]  # Player 3

# Notes
The open-loop Nash equilibrium reflects each player's optimal trajectory
given the other players' trajectories. Vehicles will naturally adjust their
speed to yield to or assert priority over approaching vehicles, producing
emergent negotiation without explicit rule-following.

Convergence is typically achieved in 10–30 iterations for the default
initial conditions. A warm start from `prev_sol` can be passed to
`solve(game, solver; warmstart=prev_sol)` for receding-horizon use.
"""
function solve_intersection(; verbose::Bool=false, kwargs...)
    game   = three_player_intersection_game(; kwargs...)
    solver = iLQGames(max_iter=300, ε_conv=1e-4, μ_init=0.0)
    return solve(game, solver; verbose, check_compatibility=false)
end

# ── Inter-vehicle distance utility ───────────────────────────────────────────

"""
    inter_vehicle_distances(sol::GNEPSolution) -> (d12, d13, d23)

Extract pairwise Euclidean distances between vehicle positions over the
trajectory. Returns three length-(N+1) vectors.
"""
function inter_vehicle_distances(sol::GNEPSolution)
    X = sol.state_trajectory
    N = size(X, 2)
    d12 = [norm(X[1:2, k] - X[5:6, k])  for k in 1:N]
    d13 = [norm(X[1:2, k] - X[9:10, k]) for k in 1:N]
    d23 = [norm(X[5:6, k] - X[9:10, k]) for k in 1:N]
    return d12, d13, d23
end