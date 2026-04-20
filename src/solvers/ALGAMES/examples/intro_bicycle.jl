# ============================================================================
# intro_bicycle.jl
#
# Three-player bicycle game — port of the Algames.jl intro example.
#
# Reference: Le Cleac'h et al. (2021)
#   github.com/RoboticExplorationLab/Algames.jl/examples/intro_example.jl
#
# ── Problem ──────────────────────────────────────────────────────────────────
#   3 players, bicycle (kinematic car) dynamics, N=20 steps, dt=0.1 s.
#   Each player drives from a staggered start to a goal while:
#     (objective)   soft proximity penalty between all pairs
#     (constraint)  hard collision avoidance, d_min = 0.08
#     (constraint)  control limits |a|, |δ| ≤ 5  (all players)
#     (constraint)  state bounds ±5  (player 1 only)
#     (constraint)  wall at y = −0.4, active for x ∈ [0, 1]
#     (constraint)  three circular keep-out zones
#
# ── Bicycle dynamics (continuous-time, integrated by ALGAMES RK4) ─────────
#   PDGNEProblem assembles the joint state by concatenating private states:
#   Joint state (12D): x = [x₁y₁v₁ψ₁ | x₂y₂v₂ψ₂ | x₃y₃v₃ψ₃]  (block layout)
#   Joint control (6D): u = [a₁δ₁ | a₂δ₂ | a₃δ₃]               (block layout)
#   Player i's position in joint state: xᵢ at (i-1)*4+1, yᵢ at (i-1)*4+2
#
#   Each player's RHS receives its PRIVATE 4D state [xᵢ,yᵢ,vᵢ,ψᵢ] and
#   2D control [aᵢ,δᵢ] (SeparableDynamics slices before calling):
#     β  = atan(lr·tan(δᵢ), lf+lr)
#     [ẋᵢ, ẏᵢ, v̇ᵢ, ψ̇ᵢ] = [vᵢcos(β+ψᵢ), vᵢsin(β+ψᵢ), aᵢ, vᵢsin(β)/lr]
#
# ── Usage ────────────────────────────────────────────────────────────────────
#   using DifferentialGamesBase, DifferentialGamesBaseSolvers
#   include("intro_bicycle.jl")
#   sol = solve_bicycle_intro(verbose=true)
#   println("Converged:  ", sol.converged)
#   println("Costs:      ", [round(get_cost(sol, i); digits=3) for i in 1:3])
# ============================================================================

using LinearAlgebra
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# ============================================================================
# Per-player bicycle RHS (used by SeparableDynamics)
# ============================================================================

"""
    _bicycle_rhs(lf, lr) -> Function

Continuous-time RHS for one bicycle player.

`SeparableDynamics.evaluate_dynamics` calls each player's function with
the player's PRIVATE state and control slices (not the full joint vectors):
  - `x` : length 4  [xᵢ, yᵢ, vᵢ, ψᵢ]  (private state)
  - `u` : length 2  [aᵢ, δᵢ]            (private control)

Returns the 4D private state derivative [ẋᵢ, ẏᵢ, v̇ᵢ, ψ̇ᵢ].
SeparableDynamics stacks all players' returns to form the 12D joint ẋ.
"""
function _bicycle_rhs(lf::Float64, lr::Float64)
    L = lf + lr
    return (x, u, params, t) -> begin
        # Private state layout: x = [xᵢ, yᵢ, vᵢ, ψᵢ]
        # Private control layout: u = [aᵢ, δᵢ]
        β = atan(lr * tan(u[2]), L)   # slip angle
        v = x[3];  ψ = x[4]
        [v * cos(β + ψ),    # ẋᵢ
         v * sin(β + ψ),    # ẏᵢ
         u[1],               # v̇ᵢ = aᵢ
         v * sin(β) / lr]   # ψ̇ᵢ
    end
end

# ============================================================================
# Problem constructor
# ============================================================================

"""
    build_bicycle_intro(; kwargs...) -> GameProblem{Float64}

Construct the 3-player bicycle intro game. All keyword arguments are
optional with defaults matching the Algames.jl reference.

| Keyword     | Default | Meaning                                       |
|:------------|:--------|:----------------------------------------------|
| `lf`        | 0.05    | Front half-wheelbase (m)                      |
| `lr`        | 0.05    | Rear half-wheelbase (m)                       |
| `N`         | 20      | Number of timesteps                           |
| `dt`        | 0.1     | Timestep (s)                                  |
| `d_min`     | 0.08    | Hard collision avoidance radius (m)           |
| `u_lim`     | 5.0     | Symmetric control bound for a and δ           |
| `x_lim`     | 5.0     | Symmetric state bound for player 1            |
| `Q_diag`    | 10.0    | Diagonal state cost weight                    |
| `R_diag`    | 0.1     | Diagonal control cost weight                  |
| `Qf_scale`  | 10.0    | Terminal cost = Qf_scale × stage Q            |
"""
function build_bicycle_intro(;
    lf      ::Float64 = 0.05,
    lr      ::Float64 = 0.05,
    N       ::Int     = 20,
    dt      ::Float64 = 0.1,
    d_min   ::Float64 = 0.08,
    u_lim   ::Float64 = 5.0,
    x_lim   ::Float64 = 5.0,
    Q_diag  ::Float64 = 10.0,
    R_diag  ::Float64 = 0.1,
    Qf_scale::Float64 = 10.0,
)
    T  = Float64
    tf = T(N * dt)
    p  = 3;  ni = 4;  mi = 2   # players, private state dim, private control dim
    n  = p * ni                  # joint state dim = 12

    # Zero-based offsets in joint vectors
    soff = (i::Int) -> (i-1) * ni    # state offset for player i
    coff = (i::Int) -> (i-1) * mi    # control offset for player i

    # Position indices in joint state (block layout from PDGNEProblem).
    # PDGNEProblem sets initial_state = vcat([p.x0 for p in players]...), so
    # player i's private state occupies indices [(i-1)*ni+1 : i*ni] in the joint
    # state — a contiguous block, NOT the interleaved layout of the reference
    # Algames.jl BicycleGame.  The private state per player is [x,y,v,ψ].
    xi_idx(i) = (i-1)*ni + 1    # x-position of player i in joint state
    yi_idx(i) = (i-1)*ni + 2    # y-position of player i in joint state

    # ── Initial and goal states (private 4D per player) ───────────────────
    # Reference x0 uses interleaved layout [x₁x₂x₃ | y₁y₂y₃ | v₁v₂v₃ | ψ₁ψ₂ψ₃]:
    #   [0.1, 0.0, 0.5,   ← x-positions
    #   -0.4, 0.0, 0.7,   ← y-positions
    #    0.0, 0.0, 0.0,   ← speeds (all zero)
    #    0.0, 0.0, 0.0]   ← headings (all zero)
    # De-interleaved into per-player private states [x, y, v, ψ]:
    x0_priv = [
        T[0.1, -0.4, 0.0, 0.0],   # player 1
        T[0.0,  0.0, 0.0, 0.0],   # player 2
        T[0.5,  0.7, 0.0, 0.0],   # player 3
    ]
    xf_priv = [
        T[2.0,  0.4, 0.0, 0.0],   # player 1 goal
        T[2.0,  0.0, 0.0, 0.0],   # player 2 goal
        T[3.0, -0.4, 0.0, 0.0],   # player 3 goal
    ]

    # ── Dynamics ──────────────────────────────────────────────────────────
    fdyn = [_bicycle_rhs(lf, lr) for _ in 1:p]
    dyn  = SeparableDynamics(fdyn, fill(ni, p), fill(mi, p))

    # ── Objectives ────────────────────────────────────────────────────────
    # Stage: quadratic tracking cost only.
    #
    # The reference adds soft proximity costs in the objective for smoother
    # convergence. We omit them here because CompositeCostTerm cost_term_gradient
    # is called with the PRIVATE state slice (4D) for SeparableDynamics games,
    # but proximity costs need cross-player position indices into the full joint
    # state (12D). The hard SharedInequality collision avoidance below handles
    # collision avoidance correctly; the solver converges without the soft term.

    objectives = map(1:p) do i
        so = soff(i);  co = coff(i)

        stage_cost = (
            QuadraticStateCost(
                Q_diag * Matrix{T}(I, ni, ni),
                xf_priv[i], so, ni
            ) +
            QuadraticControlCost(
                R_diag * Matrix{T}(I, mi, mi),
                co, mi
            )
        )

        terminal_cost = QuadraticTerminalCost(
            Qf_scale * Q_diag * Matrix{T}(I, ni, ni),
            xf_priv[i], so, ni
        )

        minimize(stage_cost; terminal=terminal_cost, player_id=i)
    end

    # ── Private constraints ────────────────────────────────────────────────

    # 1. Control bounds for all players
    ctrl_cons = [
        control_bounds(i;
            control_offset = coff(i),
            control_dim    = mi,
            lower          = fill(-u_lim, mi),
            upper          = fill( u_lim, mi),
        ) for i in 1:p
    ]

    # 2. State bounds for player 1 (all 4 private states in [-x_lim, x_lim])
    sb_p1 = state_bounds(1;
        state_offset = soff(1),
        state_dim    = ni,
        lower        = fill(-x_lim, ni),
        upper        = fill( x_lim, ni),
    )

    # 3. Wall constraint: yᵢ ≥ -0.4, active when xᵢ ∈ [0, 1]
    #    Half-plane: (pos - p1)·v ≤ 0 with p1=[0,-0.4], v=[0,-1]
    #    → -(yᵢ + 0.4) ≤ 0  when xᵢ ∈ [0,1]
    wall_cons = map(1:p) do i
        ixᵢ = xi_idx(i);  iyᵢ = yi_idx(i)
        PrivateInequality(i;
            func = (x, u, params, t) -> begin
                S = eltype(x)
                # Gate: active only when x ∈ [0, 1]
                x_val = x[ixᵢ];  y_val = x[iyᵢ]
                active = (x_val >= zero(S) && x_val <= one(S))
                c = -(y_val + S(0.4))   # -(y - (-0.4)) ≤ 0  ↔  y ≥ -0.4
                [active ? c : zero(S)]
            end,
            dim = 1
        )
    end

    # 4. Circular keep-out zones: r² - (xᵢ-cx)² - (yᵢ-cy)² ≤ 0
    xc    = T[1.0, 2.0, 3.0]
    yc    = T[1.0, 2.0, 3.0]
    r_obs = T[0.1, 0.2, 0.3]
    circle_cons = [
        PrivateInequality(i;
            func = let cx = xc[k], cy = yc[k], r = r_obs[k],
                       ixᵢ = xi_idx(i), iyᵢ = yi_idx(i)
                (x, u, params, t) -> begin
                    [r^2 - (x[ixᵢ]-cx)^2 - (x[iyᵢ]-cy)^2]
                end
            end,
            dim = 1
        )
        for k in 1:length(xc) for i in 1:p
    ]

    # ── Shared constraints ────────────────────────────────────────────────

    # Hard collision avoidance for all player pairs.
    # Cannot use DGB's collision_avoidance() because positions are non-contiguous
    # in the interleaved layout. Use SharedInequality with explicit distance.
    ε_dist = T(1e-6)
    collision_cons = [
        SharedInequality([i, j];
            func = let ixᵢ = xi_idx(i), iyᵢ = yi_idx(i),
                       ixⱼ = xi_idx(j), iyⱼ = yi_idx(j),
                       dmin = T(d_min),  ε = ε_dist
                (x, u, params, t) -> [dmin - sqrt((x[ixᵢ]-x[ixⱼ])^2 +
                                                   (x[iyᵢ]-x[iyⱼ])^2 + ε^2)]
            end,
            dim = 1
        )
        for i in 1:p for j in i+1:p
    ]

    # ── Assemble PlayerSpecs ──────────────────────────────────────────────
    players = [
        PlayerSpec(
            i, ni, mi, x0_priv[i], fdyn[i], objectives[i],
            # Private constraints for player i:
            AbstractConstraint[
                ctrl_cons[i],                                  # control bounds
                (i == 1 ? [sb_p1] : AbstractConstraint[])..., # state bounds (player 1 only)
                wall_cons[i],                                  # wall segment
                [c for c in circle_cons if c.player == i]..., # circle obstacles
            ]
        )
        for i in 1:p
    ]

    # Shared constraints: only the inter-player collision avoidance.
    # Control bounds and state bounds are private (per-player) and live in
    # PlayerSpec.constraints above. ALGAMES evaluates all constraints via
    # Iterators.flatten((game.private_constraints, game.shared_constraints)).
    game = DifferentialGame(players, tf, T(dt);
                            shared_constraints=AbstractConstraint[collision_cons...])
    return game
end

# ============================================================================
# Solver entry point
# ============================================================================

"""
    solve_bicycle_intro(; verbose=false, kwargs...) -> GNEPSolution{Float64}

Build and solve the 3-player bicycle intro problem with ALGAMES.

Problem keyword arguments (passed to `build_bicycle_intro`) and solver
keyword arguments are both accepted.

## Example
```julia
using DifferentialGamesBase, DifferentialGamesBaseSolvers
include("intro_bicycle.jl")

sol = solve_bicycle_intro(verbose=true)
println("Converged:  ", sol.converged)
println("Iterations: ", sol.iterations)
println("Solve time: ", round(sol.solve_time; digits=2), "s")
for i in 1:3
    println("  Player \$i cost: ", round(get_cost(sol, i); digits=4))
end
println("Con violation: ", sol.solver_info[:con_vio])
```
"""
function solve_bicycle_intro(;
    verbose     ::Bool    = false,
    outer_iter  ::Int     = 200,
    inner_iter  ::Int     = 20,
    tol_opt     ::Float64 = 1e-4,
    tol_dyn     ::Float64 = 1e-4,
    tol_con     ::Float64 = 1e-3,
    ρ_init      ::Float64 = 1.0,
    ρ_increase  ::Float64 = 10.0,
    ρ_max       ::Float64 = 1e8,
    # Problem parameters (forwarded to build_bicycle_intro)
    lf       ::Float64 = 0.05,
    lr       ::Float64 = 0.05,
    N        ::Int     = 20,
    dt       ::Float64 = 0.1,
    d_min    ::Float64 = 0.08,
    u_lim    ::Float64 = 5.0,
    x_lim    ::Float64 = 5.0,
    Q_diag   ::Float64 = 10.0,
    R_diag   ::Float64 = 0.1,
    Qf_scale ::Float64 = 10.0,
)
    game = build_bicycle_intro(;
        lf, lr, N, dt, d_min, u_lim, x_lim, Q_diag, R_diag, Qf_scale
    )
    solver = ALGAMES(;
        outer_iter, inner_iter,
        tol_opt, tol_dyn, tol_con,
        ρ_init, ρ_increase, ρ_max,
    )

    # ── Initial guess ─────────────────────────────────────────────────────────
    # All players start with v=0. Zero-control rollout leaves them stationary,
    # making the initial Jacobian near-singular (∂ẋ/∂u ≈ 0 via β≈0 at v=0).
    # Seed with a constant gentle acceleration so the initial trajectory has
    # nonzero speed and the Newton system is well-conditioned from iteration 1.
    #
    # a_seed (m/s²) gives each player a speed of a_seed * N * dt by step N.
    # With N=20, dt=0.1, a_seed=0.5: final speed ≈ 1 m/s — well within u_lim=5.
    a_seed = 0.5
    warmstart = _bicycle_zero_speed_warmstart(game, Float64(a_seed))
    return solve(game, solver; warmstart, verbose)
end

"""
    _bicycle_zero_speed_warmstart(game, a_seed) -> GNEPSolution

Construct a warmstart by rolling out each player with constant acceleration
`a_seed` and zero steering. This breaks the v=0 degeneracy in the bicycle
dynamics and gives the Newton solver a well-conditioned starting point.
"""
function _bicycle_zero_speed_warmstart(game::GameProblem{T}, a_seed::T) where {T}
    np   = game.n_players
    N    = n_steps(game)
    dt   = T(game.time_horizon.dt)
    cdim = game.metadata.control_dims
    n    = total_state_dim(game.dynamics)

    # Constant acceleration, zero steering for all players
    U_warm = [fill(T(0), cdim[i], N) for i in 1:np]
    for i in 1:np
        U_warm[i][1, :] .= a_seed   # aᵢ = a_seed, δᵢ = 0
    end

    # Forward rollout using SeparableDynamics.evaluate_dynamics + RK4.
    # Inlined here (not calling _rollout_flat) so this function works when the
    # example file is run directly from Main via include().
    U_flat   = vcat(U_warm...)
    X_warm   = Matrix{T}(undef, n, N+1)
    X_warm[:, 1] .= game.initial_state
    for k in 1:N
        x_k = X_warm[:, k]
        u_k = U_flat[:, k]
        # Inlined RK4 (mirrors _discrete_step for SeparableDynamics)
        f  = (x, u) -> evaluate_dynamics(game.dynamics, x, u, nothing, T((k-1)*dt))
        k1 = f(x_k,             u_k)
        k2 = f(x_k .+ dt/2 .* k1, u_k)
        k3 = f(x_k .+ dt/2 .* k2, u_k)
        k4 = f(x_k .+ dt    .* k3, u_k)
        X_warm[:, k+1] .= x_k .+ (dt/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    end

    times    = collect(range(T(0), game.time_horizon.tf, length=N+1))
    trajs    = [Trajectory{T}(i, X_warm, U_warm[i], times, T(0)) for i in 1:np]
    strategy = OpenLoopStrategy(U_warm, copy(cdim), times)
    return GNEPSolution(
        game, trajs;
        state_trajectory = X_warm,
        strategy         = strategy,
        equilibrium_type = :OpenLoopNash,
        converged        = false,
        iterations       = 0,
        solve_time       = 0.0,
        solver_info      = Dict{Symbol,Any}(),
    )
end

sol = solve_bicycle_intro(verbose=true)
plot_trajectories(sol)
