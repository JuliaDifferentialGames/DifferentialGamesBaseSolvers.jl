using LinearAlgebra
using Printf
using DifferentialGamesBase
using DifferentialGamesBaseSolvers
using Statistics

# ============================================================================
# Trajectory builder (warmstart utility)
# ============================================================================

function _build_trajectory(
    z_flat::Vector{T}, player_id::Int, n_x::Int, n_u::Int,
    N::Int, times::Vector{T}, cost::T
) where {T}
    stride   = n_x + n_u
    states   = Matrix{T}(undef, n_x, N + 1)
    controls = Matrix{T}(undef, n_u, N)
    for k in 0:(N-1)
        xs = k * stride + 1
        states[:, k+1]   = z_flat[xs : xs+n_x-1]
        controls[:, k+1] = z_flat[xs+n_x : xs+n_x+n_u-1]
    end
    states[:, N+1] = z_flat[N*stride+1 : N*stride+n_x]
    return Trajectory{T}(player_id, states, controls, times, cost)
end

# ============================================================================
# Warmstart
# ============================================================================

"""
    hallway_warmstart(game, goals, lateral_offsets; dt)

Sinusoidal lateral-bump warmstart with unicycle inverse kinematics.
Interpolates position linearly, adds lateral bump, recovers (v, ω).
"""
function hallway_warmstart(game, goals, lateral_offsets; dt=0.1)
    N   = game.time_horizon.N
    trajs = Trajectory{Float64}[]
    for i in 1:2
        n_x    = game.metadata.state_dims[i]
        n_u    = game.metadata.control_dims[i]
        stride = n_x + n_u
        off    = game.metadata.state_offsets[i]
        x0_i   = game.initial_state[(off+1):(off+n_x)]
        g_i    = goals[i]
        lat    = lateral_offsets[i]

        z = zeros((n_x + n_u) * N + n_x)
        for k in 0:N
            α   = k / N
            x_k = (1 - α) .* x0_i .+ α .* g_i
            x_k[2] += lat * sin(π * α)
            z[k*stride + 1 : k*stride + n_x] = x_k

            if k < N
                α_next = (k + 1) / N
                x_next = (1 - α_next) .* x0_i .+ α_next .* g_i
                x_next[2] += lat * sin(π * α_next)
                dx = x_next .- x_k
                θ  = x_k[3]
                v  = (dx[1]*cos(θ) + dx[2]*sin(θ)) / dt
                ω  = dx[3] / dt
                z[k*stride + n_x + 1 : (k+1)*stride] = [v, ω]
            end
        end

        times = collect(range(0.0, game.time_horizon.tf, length=N+1))
        push!(trajs, _build_trajectory(z, i, n_x, n_u, N,
                                        convert(Vector{Float64}, times), 0.0))
    end
    return WarmstartData{Float64}(trajs, nothing, Dict{Symbol,Any}())
end

# ============================================================================
# Problem definition
# ============================================================================

"""
    build_hallway_problem(; kwargs...)

2-player unicycle hallway-passing GNEP.

# Cost structure (iLQGames Eq. 9–12 + ALGAMES Eq. 15)
- Stage cost: control effort + running goal-tracking at every timestep
- Terminal cost: stronger goal-tracking at final state
- Hard constraints: AL-enforced collision avoidance + hallway bounds

Running goal cost distributes the gradient signal across the horizon,
preventing terminal cost from being overwhelmed by AL penalty terms.

# Constraints (unnormalized, ALGAMES Eq. 13–14)
- Collision: d_safe² - ‖p₁ - p₂‖² ≤ 0
- Hallway:  ±py - hw_bound ≤ 0  for each player
"""
function build_hallway_problem(;
    N::Int            = 40,
    dt::Float64       = 0.1,
    ρ_ctrl::Float64   = 0.01,
    ρ_run::Float64    = 1.0,
    ρ_goal::Float64   = 100.0,
    d_safe::Float64   = 0.3,
    hw::Float64       = 0.5,
    goal_offset::Float64 = 0.15
)
    T  = Float64
    tf = N * dt

    x0_1 = T[-1.0, -goal_offset, 0.0]
    x0_2 = T[ 1.0,  goal_offset, 0.0]
    g1   = T[ 1.0,  goal_offset, 0.0]
    g2   = T[-1.0, -goal_offset, 0.0]

    hw_bound = hw + d_safe / 1.25   # wall constraint threshold

    # Unicycle Euler dynamics: x = [px, py, θ],  u = [v, ω]
    dyn = (x, u, p, t) -> x .+ dt .* [u[1]*cos(x[3]), u[1]*sin(x[3]), u[2]]

    # ── Stage cost: control effort + time-weighted running goal cost ──────────
    # Uniform running cost opposes collision avoidance at the crossing: the goal
    # gradient at k≈N/2 pulls agents toward each other's territory. A ramp
    # α=t/N is near-zero early (agents cross freely) and near-one late
    # (strong pull to goal after crossing). This breaks the frozen-robot
    # equilibrium without fighting collision avoidance at the crossing.
    stage_1 = NonlinearStageCost(
        (x, u, p, t) -> begin
            α = T(t) / T(N)
            ρ_ctrl * (u[1]^2 + u[2]^2) + ρ_run * α * sum(abs2, x .- g1)
        end;
        is_separable = true
    )
    stage_2 = NonlinearStageCost(
        (x, u, p, t) -> begin
            α = T(t) / T(N)
            ρ_ctrl * (u[1]^2 + u[2]^2) + ρ_run * α * sum(abs2, x .- g2)
        end;
        is_separable = true
    )

    # ── Terminal cost ─────────────────────────────────────────────────────────
    term_1 = NonlinearTerminalCost((x, p) -> ρ_goal * sum(abs2, x .- g1))
    term_2 = NonlinearTerminalCost((x, p) -> ρ_goal * sum(abs2, x .- g2))

    # ── Terminal goal constraint: ‖x_N - g‖² ≤ r_goal² ──────────────────────
    # Encoded as a stage inequality active only at k=N so it enters through
    # Q[N] in the backward pass (not Qf), avoiding Riccati blow-up from
    # large ρ in Qf. Returns -1 (inactive) at all other timesteps.
    r_goal = T(0.05)
    goal_c1 = NonlinearConstraint(
        (x, u, p, t) -> [t == N ? sum(abs2, x .- g1) - r_goal^2 : -one(T)],
        1; constraint_type = :inequality
    )
    goal_c2 = NonlinearConstraint(
        (x, u, p, t) -> [t == N ? sum(abs2, x .- g2) - r_goal^2 : -one(T)],
        1; constraint_type = :inequality
    )
    priv1 = [PrivateConstraint(goal_c1, 1)]
    priv2 = [PrivateConstraint(goal_c2, 2)]

    p1 = PlayerSpec(1, 3, 2, x0_1, dyn, PlayerObjective(1, stage_1, term_1), priv1)
    p2 = PlayerSpec(2, 3, 2, x0_2, dyn, PlayerObjective(2, stage_2, term_2), priv2)

    # ── Shared constraints (unnormalized, ALGAMES Eq. 13–14) ─────────────────
    coupling_c = NonlinearConstraint(
        (x_joint, u, p, t) -> begin
            p1_xy = x_joint[1:2]
            p2_xy = x_joint[4:5]
            py1   = x_joint[2]
            py2   = x_joint[5]
            [
                d_safe^2 - sum(abs2, p1_xy .- p2_xy),
                 py1 - hw_bound,
                -py1 - hw_bound,
                 py2 - hw_bound,
                -py2 - hw_bound
            ]
        end,
        5; constraint_type = :inequality
    )
    shared = SharedConstraint(coupling_c, [1, 2])

    return PDGNEProblem([p1, p2], [shared], T(tf), T(dt))
end

# ============================================================================
# Problem + solver instantiation
# ============================================================================

game = build_hallway_problem(;
    N           = 60,
    dt          = 0.1,
    ρ_ctrl      = 0.01,
    ρ_run       = 5.0,    # ramps from 0 to 5 over horizon — light near crossing, strong near goal
    ρ_goal      = 1000.0, # very strong terminal pull
    d_safe      = 0.25,   # smaller safe distance gives geometric slack at crossing
    hw          = 0.5,
    goal_offset = 0.35    # separation at crossing = 2*0.35=0.7 > 2*d_safe=0.5 → 0.2m slack
)

solver = PANGOLIN(;
    max_iter             = 500,
    convergence_tol      = 1e-3,
    max_elwise_diff_step = 1.0,
    max_line_search_iter = 20,
    α_init               = 0.5,
    α_step               = 0.5,
    armijo_param         = 1e-4,
    state_regularization = 5.0,
    constraint_opts      = ALOptions(;
        ρ_init        = 0.1,    # start very small so early iters converge unconstrained
        ρ_max         = 200.0,
        φ             = 1.2,    # slow growth: reaches ρ_max in ~50 iters from ρ_init=0.1
        violation_tol = 1e-3
    )
)

g1_ws = [1.0,  0.35, 0.0]
g2_ws = [-1.0, -0.35, 0.0]
# Large lateral offset (0.6) pushes crossing geometry wide; the sinusoidal
# bump peaks at the midpoint (k=N/2=30), but agents are already well-separated
# at the actual crossing point (k≈10-15) because goal_offset=0.35 gives
# 0.7m separation. This warmstart keeps agents far apart through the crossing.
ws = hallway_warmstart(game, [g1_ws, g2_ws], [0.6, -0.6]; dt=0.1)

sol = solve(game, solver; warmstart=ws, verbose=true)

# ============================================================================
# Diagnostics
# ============================================================================

function check_dynamic_consistency(sol, game)
    println("\n=== Dynamic Consistency Check ===")
    dt  = game.time_horizon.dt
    N   = game.time_horizon.N
    for i in 1:game.n_players
        traj = sol.trajectories[i]
        max_defect = 0.0; worst_k = 0
        for k in 1:N
            x_k = traj.states[:, k]; u_k = traj.controls[:, k]
            x_true = x_k .+ dt .* [u_k[1]*cos(x_k[3]), u_k[1]*sin(x_k[3]), u_k[2]]
            d = maximum(abs.(x_true .- traj.states[:, k+1]))
            if d > max_defect; max_defect = d; worst_k = k; end
        end
        println("Player ", i, ": max defect = ", round(max_defect, sigdigits=6),
                "  (k=", worst_k, ")")
    end
end

function check_costs(sol, game)
    # Uses the actual problem parameters — no hardcoded values
    println("\n=== Manual Cost Check ===")
    N        = game.time_horizon.N
    dt       = game.time_horizon.dt
    offsets  = game.metadata.state_offsets
    dims     = game.metadata.state_dims
    for i in 1:game.n_players
        traj = sol.trajectories[i]
        sc   = game.objectives[i].stage_cost
        tc   = game.objectives[i].terminal_cost
        stage_sum = 0.0
        for k in 1:N
            xi = traj.states[:, k]
            ui = traj.controls[:, k]
            stage_sum += evaluate_stage_cost(sc, xi, ui, nothing, k)
        end
        xf  = traj.states[:, end]
        term = evaluate_terminal_cost(tc, xf, nothing)
        expected = stage_sum + term
        stored   = traj.cost
        println("Player ", i,
                ":  stage=",   round(stage_sum, sigdigits=5),
                "  terminal=", round(term,      sigdigits=5),
                "  expected=", round(expected,  sigdigits=5),
                "  stored=",   round(stored,    sigdigits=5),
                "  MISMATCH=", abs(expected - stored) > 0.5)
    end
end

function check_convergence_trigger(sol, solver)
    println("\n=== Convergence Trigger ===")
    iter_costs = sol.solver_info[:iter_costs]
    viols      = sol.solver_info[:violation_history]
    ρs         = sol.solver_info[:ρ_history]
    println("Iterations: ", length(iter_costs), "  converged: ", sol.converged)
    println("tol=", solver.convergence_tol, "  max_dev=", solver.max_elwise_diff_step,
            "  β=", solver.armijo_param)
    # Print first 10 and last 5
    show_iters = vcat(1:min(10, length(iter_costs)),
                      max(11, length(iter_costs)-4):length(iter_costs))
    for k in unique(show_iters)
        k > length(iter_costs) && continue
        println("  iter ", k, ":  costs=", round.(iter_costs[k], sigdigits=5),
                "  viol=", round(viols[k], sigdigits=4),
                "  ρ=",   round(ρs[k],    sigdigits=4))
    end
    if length(iter_costs) >= 2
        Δ = maximum(abs.(iter_costs[end] .- iter_costs[end-1]))
        println("Cost change (last 2 iters): ", round(Δ, sigdigits=4))
    end
end

function check_warmstart(ws, game; d_safe=0.25, hw=0.5)
    println("\n=== Warmstart Check ===")
    N     = game.time_horizon.N
    bound = hw + d_safe / 1.25
    for i in 1:2
        traj = ws.trajectories[i]
        max_py = maximum(abs.(traj.states[2, :]))
        println("Player ", i,
                ":  x0=",  round.(traj.states[:, 1],     sigdigits=4),
                "  xmid=", round.(traj.states[:, N÷2+1], sigdigits=4),
                "  xf=",   round.(traj.states[:, end],   sigdigits=4))
        println("   max|py|=", round(max_py, sigdigits=4),
                "  bound=", round(bound, sigdigits=4), "  ",
                max_py <= bound ? "✓ feasible" : "✗ INFEASIBLE",
                "  max|u|=", round(maximum(abs.(traj.controls)), sigdigits=4))
    end
    println("\nCollision check (d_safe=", d_safe, "):")
    any_coll = false
    for k in 1:N+1
        d = norm(ws.trajectories[1].states[1:2, k] .- ws.trajectories[2].states[1:2, k])
        if d < d_safe
            println("  COLLISION k=", k, " dist=", round(d, sigdigits=4))
            any_coll = true
        end
    end
    any_coll || println("  No collisions ✓")
end

function check_constraints(sol, game; d_safe=0.25, hw=0.5)
    println("\n=== Constraint Violations Along Solution ===")
    N     = game.time_horizon.N
    bound = hw + d_safe / 1.25
    coll  = Float64[]; h1 = Float64[]; h2 = Float64[]
    for k in 1:N+1
        p1  = sol.trajectories[1].states[1:2, k]
        p2  = sol.trajectories[2].states[1:2, k]
        py1 = sol.trajectories[1].states[2, k]
        py2 = sol.trajectories[2].states[2, k]
        push!(coll, max(0.0, d_safe^2 - sum(abs2, p1 .- p2)))
        push!(h1,   max(0.0, abs(py1) - bound))
        push!(h2,   max(0.0, abs(py2) - bound))
    end
    println("Collision:  max=", round(maximum(coll), sigdigits=4),
            "  mean=", round(mean(coll), sigdigits=4))
    println("Hallway P1: max=", round(maximum(h1),   sigdigits=4),
            "  mean=", round(mean(h1),   sigdigits=4))
    println("Hallway P2: max=", round(maximum(h2),   sigdigits=4),
            "  mean=", round(mean(h2),   sigdigits=4))

    # Goal inference: agents swap sides including y-offset.
    # Player 1: x0=[-1,-off,0] → goal=[1,off,0] = [-x0[1], -x0[2], 0]
    # Player 2: x0=[1,off,0]   → goal=[-1,-off,0]
    r_goal = 0.05
    xf1  = sol.trajectories[1].states[:, end]
    xf2  = sol.trajectories[2].states[:, end]
    x0_1 = sol.trajectories[1].states[:, 1]
    x0_2 = sol.trajectories[2].states[:, 1]
    goal1 = [-x0_1[1], -x0_1[2], 0.0]
    goal2 = [-x0_2[1], -x0_2[2], 0.0]
    dist1_sq = sum(abs2, xf1 .- goal1)
    dist2_sq = sum(abs2, xf2 .- goal2)
    g1_viol = max(0.0, dist1_sq - r_goal^2)
    g2_viol = max(0.0, dist2_sq - r_goal^2)
    println("Goal P1:    ‖xf-g‖=", round(sqrt(dist1_sq), sigdigits=4),
            "  viol=", round(g1_viol, sigdigits=4),
            (g1_viol < 1e-6 ? "  ✓" : "  ✗ NOT REACHED"))
    println("Goal P2:    ‖xf-g‖=", round(sqrt(dist2_sq), sigdigits=4),
            "  viol=", round(g2_viol, sigdigits=4),
            (g2_viol < 1e-6 ? "  ✓" : "  ✗ NOT REACHED"))

    if maximum(coll) > 1e-4 || maximum(h1) > 1e-4
        println("Per-timestep violations > 1e-4:")
        for k in 1:N+1
            (coll[k] > 1e-4 || h1[k] > 1e-4 || h2[k] > 1e-4) || continue
            println("  k=", k, "  coll=", round(coll[k], sigdigits=3),
                    "  hall1=", round(h1[k], sigdigits=3),
                    "  hall2=", round(h2[k], sigdigits=3))
        end
    end
end

function check_trajectories(sol)
    println("\n=== Solution Trajectory Summary ===")
    for i in 1:length(sol.trajectories)
        traj = sol.trajectories[i]
        println("Player ", i, ":")
        println("  x0 = ", round.(traj.states[:, 1],   sigdigits=4))
        println("  xf = ", round.(traj.states[:, end], sigdigits=4))
        println("  max|x|=", round(maximum(abs.(traj.states)),   sigdigits=4),
                "  max|u|=", round(maximum(abs.(traj.controls)), sigdigits=4))
        println("  cost (stored) = ", round(traj.cost, sigdigits=5))
        (any(isnan, traj.states) || any(isnan, traj.controls)) &&
            println("  ⚠ NaN detected!")
    end
end

println("\nRunning diagnostics...")
check_warmstart(ws, game)
check_dynamic_consistency(sol, game)
check_costs(sol, game)
check_convergence_trigger(sol, solver)
check_constraints(sol, game)
check_trajectories(sol)