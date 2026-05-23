# ============================================================================
# miller_mitra_example.jl
#
# Reproduction of the multi-agent navigation example from:
#   Miller, K. & Mitra, S. (2022).
#   "Multi-agent motion planning using differential games with
#    lexicographic preferences." IEEE CDC 2022, pp. 5751–5757.
#
# ── Problem ──────────────────────────────────────────────────────────────────
#   4 agents navigate a shared 2D workspace. Each agent must reach a
#   target position while avoiding collisions with all others. Collision
#   avoidance has strict priority (lexicographic first), followed by
#   minimizing distance to goal and control effort.
#
# ── Dynamics (double integrator, 4D state) ───────────────────────────────────
#   State:    xᵢ = (pxᵢ, pyᵢ, vxᵢ, vyᵢ) ∈ ℝ⁴   position + velocity
#   Control:  uᵢ = (axᵢ, ayᵢ) ∈ ℝ²              acceleration
#   ODE:      ṗ = v,  v̇ = u
#   Discretized by Euler with dt = 0.1 s.
#
# ── Initial positions and goals ───────────────────────────────────────────────
#   Agents start on a unit circle and swap to the opposite position:
#     Agent 1: start (1,0)   → goal (-1,0)
#     Agent 2: start (-1,0)  → goal (1,0)
#     Agent 3: start (0,1)   → goal (0,-1)
#     Agent 4: start (0,-1)  → goal (0,1)
#   All start at rest (v = 0). Control horizon N = 50, dt = 0.1 s.
#
# ── Cost structure ────────────────────────────────────────────────────────────
#   Collision cost fᵢⱼ(zᵢ, zⱼ):
#     Quadratic barrier on proximity: Σₖ max(0, d²_safe − ‖pᵢ(k)−pⱼ(k)‖²)²
#     d_safe = 0.3 m  (soft collision avoidance)
#
#   Personal cost gᵢ(zᵢ):
#     Goal-tracking + control regularization:
#     Σₖ w_goal ‖pᵢ(k) − p_goal_i‖² + w_ctrl ‖uᵢ(k)‖²
#
# ── Usage ────────────────────────────────────────────────────────────────────
#   using DifferentialGamesBase, DifferentialGamesBaseSolvers
#   include("miller_mitra_example.jl")
#   sol = solve_miller_mitra(; verbose=true)
#   println("Converged:  ", sol.converged)
#   println("Iterations: ", sol.iterations)
#   println("J_col:      ", round.(sol.solver_info[:J_col]; digits=4))
#   println("J_per:      ", round.(sol.solver_info[:J_per]; digits=4))
# ============================================================================

using LinearAlgebra
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# ── Parameters ───────────────────────────────────────────────────────────────

const MM_N      = 50        # horizon steps
const MM_DT     = 0.1       # time step (s)
const MM_TF     = MM_N * MM_DT
const MM_DSAFE  = 0.3       # minimum safe separation (m)
const MM_WGOAL  = 1.0       # goal-tracking weight
const MM_WCTRL  = 0.1       # control regularization weight
const MM_NPLAYERS = 4

# ── Initial states and goals ─────────────────────────────────────────────────
#   State layout per player: [px, py, vx, vy]

const MM_X0S = [
    [1.0, 0.0, 0.0, 0.0],    # Agent 1
    [-1.0, 0.0, 0.0, 0.0],   # Agent 2
    [0.0, 1.0, 0.0, 0.0],    # Agent 3
    [0.0, -1.0, 0.0, 0.0],   # Agent 4
]

const MM_GOALS = [
    [-1.0, 0.0],   # Agent 1 goal
    [1.0, 0.0],    # Agent 2 goal
    [0.0, -1.0],   # Agent 3 goal
    [0.0, 1.0],    # Agent 4 goal
]

# ── Double integrator dynamics ────────────────────────────────────────────────
#   Private state: [px, py, vx, vy]
#   Private control: [ax, ay]

const _mm_dynamics = (x, u, p, t) -> [x[3], x[4], u[1], u[2]]

# ── Cost functions ────────────────────────────────────────────────────────────
#
# z_i = (X_i, U_i) where:
#   X_i :: Matrix (4 × N+1)  — state trajectory (rows = [px, py, vx, vy])
#   U_i :: Matrix (2 × N)    — control trajectory (rows = [ax, ay])

"""
    _mm_collision_cost_pair(d_safe) -> Function

Returns fᵢⱼ(zᵢ, zⱼ) = Σₖ max(0, d²_safe − ‖pᵢ(k) − pⱼ(k)‖²)²

This is symmetric (fᵢⱼ = fⱼᵢ) and zero when separation ≥ d_safe everywhere.
"""
function _mm_collision_cost_pair(d_safe::Float64)
    d2 = d_safe^2
    return (z_i, z_j) -> begin
        X_i, _ = z_i
        X_j, _ = z_j
        N1 = size(X_i, 2)
        val = zero(eltype(X_i))
        for k in 1:N1
            Δx = X_i[1, k] - X_j[1, k]
            Δy = X_i[2, k] - X_j[2, k]
            sep2 = Δx^2 + Δy^2
            viol = max(zero(eltype(X_i)), d2 - sep2)
            val += viol^2
        end
        return val
    end
end

"""
    _mm_personal_cost(goal, w_goal, w_ctrl) -> Function

Returns gᵢ(zᵢ) = Σₖ w_goal ‖pᵢ(k) − goal‖² + Σₖ w_ctrl ‖uᵢ(k)‖²

The velocity terms are omitted (agents are allowed to drift freely).
"""
function _mm_personal_cost(goal::Vector{Float64}, w_goal::Float64, w_ctrl::Float64)
    return (z_i,) -> begin
        X_i, U_i = z_i
        N1 = size(X_i, 2); N2 = size(U_i, 2)
        val = zero(eltype(X_i))
        for k in 1:N1
            val += w_goal * ((X_i[1, k] - goal[1])^2 + (X_i[2, k] - goal[2])^2)
        end
        for k in 1:N2
            val += w_ctrl * (U_i[1, k]^2 + U_i[2, k]^2)
        end
        return val
    end
end

# ── Problem construction ──────────────────────────────────────────────────────

"""
    build_miller_mitra(; d_safe, w_goal, w_ctrl, tf, dt, N) -> LexicographicGameProblem

Construct the 4-agent circle-crossing lexicographic game.
"""
function build_miller_mitra(;
    d_safe  ::Float64 = MM_DSAFE,
    w_goal  ::Float64 = MM_WGOAL,
    w_ctrl  ::Float64 = MM_WCTRL,
    tf      ::Float64 = MM_TF,
    dt      ::Float64 = MM_DT,
)
    T = Float64
    n = MM_NPLAYERS

    # Per-player specs: dynamics, initial state, placeholder LQ objective.
    # LIBR does NOT use the forward_game objectives — they are placeholders
    # required by the GameProblem structure.
    players = map(1:n) do i
        Q  = zeros(T, 4, 4)
        R  = Matrix{T}(I, 2, 2) .* 1e-6   # tiny placeholder R (must be PD)
        Qf = zeros(T, 4, 4)
        stage    = LQStageCost(Q, R)
        terminal = LQTerminalCost(Qf)
        obj      = PlayerObjective(i, stage, terminal)
        PlayerSpec{T}(i, 4, 2, MM_X0S[i], _mm_dynamics, obj, T[])
    end

    # Symmetric collision cost matrix (zeros on diagonal, same f_ij for all pairs)
    f_pair = _mm_collision_cost_pair(d_safe)
    col_pairs = Matrix{Function}(undef, n, n)
    for i in 1:n, j in 1:n
        col_pairs[i, j] = (i == j) ? ((zi, zj) -> 0.0) : f_pair
    end

    # Per-player personal costs (explicit Function element type for dispatch)
    pers = Function[_mm_personal_cost(MM_GOALS[i], w_goal, w_ctrl) for i in 1:n]

    return LexicographicGameProblem(players, col_pairs, pers, T(tf), T(dt))
end

# ── Main solve function ───────────────────────────────────────────────────────

"""
    solve_miller_mitra(; verbose=false, kwargs...) -> GNEPSolution

Build and solve the 4-agent Miller & Mitra circle-crossing example using L-IBR.

# Example
```julia
sol = solve_miller_mitra(verbose=true)
println("Converged:  ", sol.converged)
println("J_col each: ", round.(sol.solver_info[:J_col]; digits=4))
println("J_per each: ", round.(sol.solver_info[:J_per]; digits=4))
```
"""
function solve_miller_mitra(;
    verbose ::Bool    = false,
    max_iter::Int     = 50,
    inner_iter::Int   = 100,
    step_size::Float64 = 0.005,
    kwargs...
)
    game   = build_miller_mitra(; kwargs...)
    solver = LIBR(; max_iter, inner_iter, step_size)
    return solve(game, solver; verbose)
end
