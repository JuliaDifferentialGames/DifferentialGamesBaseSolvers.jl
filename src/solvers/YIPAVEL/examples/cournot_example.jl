# ============================================================================
# cournot_example.jl
#
# Reproduction of the networked Cournot competition example from:
#   Yi, P. & Pavel, L. (2017). "A distributed primal-dual algorithm for
#   computation of generalized Nash equilibria via operator splitting."
#   IEEE CDC 2017, Section VI (pp. 3845–3846).
#
# ── Problem ──────────────────────────────────────────────────────────────────
#   N companies (F₁,…,Fₙ) compete on m markets (M₁,…,Mₘ).
#   Each company i controls xᵢ ∈ Ωᵢ ⊂ ℝ^{nᵢ}: its production per market.
#   Total market supply: Ax = Σᵢ Aᵢ xᵢ ∈ ℝ^m (must not exceed capacity r).
#
#   Player i's objective (minimise cost minus revenue):
#     fᵢ(xᵢ, x₋ᵢ) = cᵢ(xᵢ) − P(Ax)ᵀ Aᵢ xᵢ
#   where
#     cᵢ(xᵢ) = πᵢ ‖xᵢ‖²  +  bᵢᵀ xᵢ          (production cost, strictly convex)
#     P(y)   = P̄ − D·y                         (linear inverse-demand price)
#
#   Shared coupling constraint: Ax ≤ r (market capacity)
#   Private constraint:         0 ≤ xᵢ ≤ Θᵢ   (production upper bound)
#
# ── Implementation notes ──────────────────────────────────────────────────────
#   This example uses N = 5 companies and m = 2 markets (scaled down from the
#   paper's N=20/m=7 for tractability).  All companies participate in both
#   markets with identity market-selection matrices Aᵢ = I₂.
#   Parameters are fixed (not random) for reproducibility.
#
# ── Usage ────────────────────────────────────────────────────────────────────
#   using DifferentialGamesBase, DifferentialGamesBaseSolvers
#   include("cournot_example.jl")
#   sol = solve_cournot(verbose=true)
#   println("Converged:  ", sol.converged)
#   println("Iterations: ", sol.iterations)
#   println("Δx_hist:    ", round.(sol.solver_info[:Δ_hist][end]; digits=8))
# ============================================================================

using LinearAlgebra
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# ── Game parameters ───────────────────────────────────────────────────────────

const CN_NPLAYERS = 5         # number of companies
const CN_NMARKETS = 2         # number of markets

# Production cost coefficients (πᵢ — quadratic convexity, strictly > 0)
const CN_PI = [2.0, 1.5, 2.5, 1.8, 2.2]

# Linear cost coefficients bᵢ ∈ ℝ^m  (one row per player)
const CN_B = [
    [0.5, 0.3],   # company 1
    [0.4, 0.6],   # company 2
    [0.7, 0.2],   # company 3
    [0.3, 0.5],   # company 4
    [0.6, 0.4],   # company 5
]

# Production upper bounds Θᵢ ∈ ℝ^m (per player per market)
const CN_THETA = [
    [8.0, 10.0],
    [9.0,  8.0],
    [7.0, 11.0],
    [8.5,  9.5],
    [7.5,  8.5],
]

# Market capacity r ∈ ℝ^m (total supply ≤ r per market)
const CN_R = [18.0, 20.0]

# Maximum price vector P̄ ∈ ℝ^m
const CN_PBAR = [12.0, 10.0]

# Price sensitivity matrix D ∈ ℝ^{m×m} (diagonal)
const CN_D = Diagonal([0.5, 0.4])

# ── Market-selection matrices ─────────────────────────────────────────────────
#   All companies participate in all markets: Aᵢ = Iₘ for all i.
#   Total supply: Ax = Σᵢ xᵢ.

const CN_A_I = Matrix{Float64}(I, CN_NMARKETS, CN_NMARKETS)

# ── Cost functions ────────────────────────────────────────────────────────────

"""
    _cn_cost(i, segs...) -> Real

Player i's objective fᵢ(x₁,…,xₙ):
  fᵢ = πᵢ ‖xᵢ‖² + bᵢᵀxᵢ − P(Ax)ᵀ Aᵢxᵢ

where `segs[j]` is player j's vectorised control vec(Uⱼ) = xⱼ (for the
1-step static game, the single control step equals the decision variable).
"""
function _cn_cost(i::Int, segs::AbstractVector{<:AbstractVector})
    xi  = segs[i]
    Ax  = sum(segs)               # Σⱼ xⱼ  (Aⱼ = I for all j)
    P   = CN_PBAR .- CN_D * Ax   # P(Ax) = P̄ − D·Ax
    return CN_PI[i] * dot(xi, xi) + dot(CN_B[i], xi) - dot(P, xi)
end

# ── Game construction ─────────────────────────────────────────────────────────

"""
    build_cournot(; N, m, r, kwargs...) -> (game, cost_fns, coupling_A, coupling_b, lb, ub)

Construct the N-player, m-market Cournot game as a `ConvexGameProblem{Float64}`.

# Returns
- `game`       : `ConvexGameProblem` with placeholder LQ objectives and separable
                 trivial dynamics.  Real objectives are in `cost_fns`.
- `cost_fns`   : Per-player closures `fᵢ(seg₁,…,segₙ) -> Real`.
- `coupling_A` : Full coupling matrix `[I | I | … | I] ∈ ℝ^{m × N·m}`.
- `coupling_b` : Market capacity vector `r ∈ ℝ^m`.
- `lb`         : Per-player lower bounds (zero vectors).
- `ub`         : Per-player upper bounds (Θᵢ vectors).
"""
function build_cournot(;
    N      ::Int     = CN_NPLAYERS,
    m      ::Int     = CN_NMARKETS,
    π_vals ::Vector{Float64} = CN_PI[1:N],
    b_vals ::Vector{Vector{Float64}} = CN_B[1:N],
    Θ_vals ::Vector{Vector{Float64}} = CN_THETA[1:N],
    r      ::Vector{Float64} = CN_R,
    P̄      ::Vector{Float64} = CN_PBAR,
    D      ::AbstractMatrix{Float64} = CN_D,
)
    T  = Float64
    dt = T(1.0)
    tf = T(1.2)   # N_steps = round(1.2/1.0) = 1; dt < tf ✓

    # ── Trivial dynamics: ẋᵢ = 0 (state is irrelevant) ─────────────────────
    # Decision variable = control uᵢ ∈ ℝ^m at the single time step.
    # With x₀ = 0, ẋ = 0 gives X = [0 | 0], U = reshape(xᵢ, m, 1).
    dyn = (x, u, p, t) -> zeros(T, m)

    # ── Players with placeholder LQ objectives ───────────────────────────────
    # Real objectives are provided via cost_fns; LQ placeholder has Q=0, R≈0.
    players = map(1:N) do i
        Q  = zeros(T, m, m)
        R  = Matrix{T}(I, m, m) * 1e-8   # tiny but positive definite
        Qf = zeros(T, m, m)
        stage    = LQStageCost(Q, R)
        terminal = LQTerminalCost(Qf)
        obj      = PlayerObjective(i, stage, terminal)
        PlayerSpec{T}(i, m, m, zeros(T, m), dyn, obj, T[])
    end

    game = ConvexGameProblem(players, tf, dt; assume_convex=true)

    # ── Cost functions ─────────────────────────────────────────────────────
    # Capture parameters explicitly to avoid closure issues with loop variable i.
    cost_fns = Function[
        let ii = i, πi = π_vals[i], bi = b_vals[i]
            (segs...) -> begin
                xi  = segs[ii]
                Ax  = sum(segs)          # Σⱼ xⱼ
                Py  = P̄ .- D * Ax       # price vector
                πi * dot(xi, xi) + dot(bi, xi) - dot(Py, xi)
            end
        end
        for i in 1:N
    ]

    # ── Coupling constraint: Ax ≤ r  (market capacity) ──────────────────────
    # Full coupling matrix [I | I | … | I] ∈ ℝ^{m × N·m}  (since Aᵢ = Iₘ)
    coupling_A = hcat(fill(Matrix{T}(I, m, m), N)...)
    coupling_b = r

    # ── Box constraints ────────────────────────────────────────────────────
    lb = [zeros(T, m) for _ in 1:N]         # xᵢ ≥ 0
    ub = [copy(Θ_vals[i]) for i in 1:N]     # xᵢ ≤ Θᵢ

    return game, cost_fns, coupling_A, coupling_b, lb, ub
end

# ── Solve wrapper ─────────────────────────────────────────────────────────────

"""
    solve_cournot(; verbose, kwargs...) -> GNEPSolution

Build and solve the N-company, m-market networked Cournot competition example
from Yi & Pavel (2017, Section VI) using the distributed primal-dual algorithm.

# Example
```julia
sol = solve_cournot(verbose=true)
println("Converged:  ", sol.converged)
println("Iterations: ", sol.iterations)
for i in 1:5
    xi = sol.trajectories[i].controls[:, 1]
    println("Company ", i, ": x = ", round.(xi; digits=3))
end
```
"""
function solve_cournot(;
    verbose   ::Bool    = false,
    max_iter  ::Int     = 5000,
    tol       ::Float64 = 1e-6,
    τ         ::Float64 = 0.03,
    ν         ::Float64 = 0.01,
    σ         ::Float64 = 0.03,
    N         ::Int     = CN_NPLAYERS,
    m         ::Int     = CN_NMARKETS,
    kwargs...
)
    game, cost_fns, coup_A, coup_b, lb, ub =
        build_cournot(; N=N, m=m, kwargs...)

    solver = YiPavel(; max_iter=max_iter, tol=tol, τ=τ, ν=ν, σ=σ)

    return solve(
        game, solver;
        verbose      = verbose,
        cost_fns     = cost_fns,
        coupling_A   = coup_A,
        coupling_b   = coup_b,
        coupling_leq = true,
        lb           = lb,
        ub           = ub
    )
end
