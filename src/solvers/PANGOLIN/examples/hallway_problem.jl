using LinearAlgebra
using Printf
using DifferentialGamesBase
using DifferentialGamesBaseSolvers
using Statistics

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


"""
Generate a straight-line interpolated warm start for the hallway problem.

Matches the toy solver's `init_interpolate` with unicycle inverse kinematics:
  - Interpolate px linearly from x0 to goal
  - Add sinusoidal lateral bump (lateral_offset * sin(π·α)) to py
  - Recover (v, ω) from the position difference via unicycle inverse kinematics
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

"""
Build the 2-player hallway passing problem as a PDGNEProblem.

Unicycle dynamics with Euler integration. Coupling constraints encoded as
NonlinearConstraint objects on the joint state at each timestep, but here
we use the shared constraint path so FALCON handles them via μ dual updates.

The coupling function evaluates:
  - Collision avoidance: d_safe² - ‖p₁-p₂‖² ≤ 0
  - Hallway upper bound:  py - (hw + d_safe/1.25) ≤ 0
  - Hallway lower bound: -py - (hw + d_safe/1.25) ≤ 0
"""
function build_hallway_problem(;
    N::Int     = 20,
    dt::Float64 = 0.1,
    ρ_ctrl::Float64 = 0.1,
    d_safe::Float64 = 0.2,
    hw::Float64     = 0.5,   # hallway half-width
    ρ_goal::Float64 = 50.0,
    goal_offset::Float64 = 0.0
)
    T = Float64
    tf = N * dt

    x0_1 = [-1.0, -goal_offset, 0.0]
    x0_2 = [ 1.0,  goal_offset, 0.0]
    g1   = [ 1.0,  goal_offset, 0.0]
    g2   = [-1.0, -goal_offset, 0.0]

    # Unicycle: x = [px, py, θ],  u = [v, ω]
    dyn = (x, u, p, t) -> x .+ dt .* [u[1]*cos(x[3]), u[1]*sin(x[3]), u[2]]

    # Stage cost: control effort only (goal encoded in terminal cost)
    stage_i(i) = NonlinearStageCost(
        (x, u, p, t) -> ρ_ctrl * (u[1]^2 + u[2]^2);   # scalar return
        is_separable = true
    )

    # Terminal cost: soft goal reaching
    term_1 = NonlinearTerminalCost(
        (x, p) -> ρ_goal * sum(abs2, x .- g1)
    )
    term_2 = NonlinearTerminalCost(
        (x, p) -> ρ_goal * sum(abs2, x .- g2)
    )

    # Initial trajectories: straight-line interpolation with lateral offset
    # (zero controls — FALCON linearizes and solves from there)

    p1 = PlayerSpec(1, 3, 2, x0_1, dyn,
                    PlayerObjective(1, stage_i(1), term_1))
    p2 = PlayerSpec(2, 3, 2, x0_2, dyn,
                    PlayerObjective(2, stage_i(2), term_2))

    # Shared coupling constraint: evaluated on joint state [x1; x2] at each t.
    # dim = 3 per timestep: (collision, hallway_upper, hallway_lower)
    # NOTE: this is a single NonlinearConstraint whose .func receives the full
    # stacked joint state [x1[1:3]; x2[1:3]] and the player's own control.
    # Normalize collision constraint to have unit gradient magnitude
    coupling_c = NonlinearConstraint(
        (x_joint, u, p, t) -> begin
            p1_xy  = x_joint[1:2]
            p2_xy  = x_joint[4:5]
            py1    = x_joint[2]
            py2    = x_joint[5]
            dist²  = sum(abs2, p1_xy .- p2_xy)
            d_safe² = d_safe^2
            [
                (d_safe² - dist²) / d_safe²,        # normalized: O(1) when barely violated
                py1 / (hw + d_safe/1.25) - 1.0,    # normalized upper bound
                -py1 / (hw + d_safe/1.25) - 1.0,    # normalized lower bound
                py2 / (hw + d_safe/1.25) - 1.0,
                -py2 / (hw + d_safe/1.25) - 1.0
            ]
        end,
        5; constraint_type = :inequality
    )
    shared = SharedConstraint(coupling_c, [1, 2])

    return PDGNEProblem([p1, p2], [shared], T(tf), T(dt))
end



# Reduced horizon for test speed: N=10, dt=0.2 → tf=2s
game = build_hallway_problem(;
    N        = 30,      # longer horizon: 3s at dt=0.1
    dt       = 0.1,
    ρ_ctrl   = 0.01,    # less control penalty → agents move faster
    ρ_goal   = 100.0,   # stronger goal attraction
    d_safe   = 0.3,     # slightly larger safe distance
    hw       = 0.5
)

game = build_hallway_problem(;
    N=30, dt=0.1, ρ_ctrl=0.05, ρ_goal=200.0, d_safe=0.3, hw=0.5
)

solver = PANGOLIN(;
    max_iter             = 500,
    convergence_tol      = 1e-3,
    max_elwise_diff_step = 5.0,
    max_line_search_iter = 30,
    state_regularization = 1e-3,
    constraint_opts      = ALOptions(;
        ρ_init        = 1.0,
        ρ_max         = 20.0,   # ρ_goal=200, so ρ_max = ρ_goal/10
        φ             = 1.5,
        violation_tol = 0.1
    )
)
# Warm start: straight-line interpolation
g1 = [1.0, 0.0, 0.0]; g2 = [-1.0, 0.0, 0.0]
ws = hallway_warmstart(game, [g1, g2], [0.4, -0.4]; dt=0.1)


sol = solve(game, solver; warmstart=ws, verbose=true)

# ── Diagnostic 1: verify operating point is dynamically consistent ────────────
function check_dynamic_consistency(sol, game)
    println("\n=== Dynamic Consistency Check ===")
    N   = game.time_horizon.N
    dt  = game.time_horizon.dt
    n_p = game.n_players
 
    for i in 1:n_p
        traj = sol.trajectories[i]
        max_defect = 0.0
        worst_k    = 0
        for k in 1:N
            x_k = traj.states[:, k]
            u_k = traj.controls[:, k]
            # Euler unicycle step
            x_next_true = x_k .+ dt .* [u_k[1]*cos(x_k[3]),
                                          u_k[1]*sin(x_k[3]),
                                          u_k[2]]
            x_next_traj = traj.states[:, k+1]
            defect = maximum(abs.(x_next_true .- x_next_traj))
            if defect > max_defect
                max_defect = defect
                worst_k    = k
            end
        end
        println("Player ", i, ": max dynamic defect = ", round(max_defect, sigdigits=6),
                "  (worst at k=", worst_k, ")")
    end
end
 
# ── Diagnostic 2: recompute costs from scratch ────────────────────────────────
function check_costs(sol, game; ρ_ctrl=0.05, ρ_goal=200.0,
                     goals=[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    println("\n=== Manual Cost Check ===")
    N   = game.time_horizon.N
    n_p = game.n_players
 
    for i in 1:n_p
        traj = sol.trajectories[i]
        stage_sum = 0.0
        for k in 1:N
            u = traj.controls[:, k]
            stage_sum += ρ_ctrl * (u[1]^2 + u[2]^2)
        end
        xf   = traj.states[:, end]
        term = ρ_goal * sum(abs2, xf .- goals[i])
        expected = stage_sum + term
        stored   = traj.cost
        println("Player ", i,
                ":  stage=",    round(stage_sum, sigdigits=5),
                "  terminal=",  round(term,      sigdigits=5),
                "  expected=",  round(expected,  sigdigits=5),
                "  stored=",    round(stored,    sigdigits=5),
                "  MISMATCH=",  abs(expected - stored) > 0.01)
    end
end
 
# ── Diagnostic 3: inspect convergence trigger ─────────────────────────────────
function check_convergence_trigger(sol, solver)
    println("\n=== Convergence Trigger ===")
    iter_costs  = sol.solver_info[:iter_costs]
    viols       = sol.solver_info[:violation_history]
    ρs          = sol.solver_info[:ρ_history]
 
    println("Iterations recorded:  ", length(iter_costs))
    println("convergence_tol:      ", solver.convergence_tol)
    println("max_elwise_diff_step: ", solver.max_elwise_diff_step)
    println("Converged:            ", sol.converged)
 
    for (k, c) in enumerate(iter_costs)
        println("  iter ", k, ":  costs=", round.(c, sigdigits=5),
                "  violation=", round(viols[k], sigdigits=4),
                "  rho=", ρs[k])
    end
 
    # Check if consecutive costs are nearly identical (premature convergence)
    if length(iter_costs) >= 2
        delta = maximum(abs.(iter_costs[end] .- iter_costs[end-1]))
        println("\nCost change last→second-last iter: ", round(delta, sigdigits=4))
    end
end
 
# ── Diagnostic 4: inspect warmstart ──────────────────────────────────────────
function check_warmstart(ws, game; d_safe=0.2, hw=0.5)
    println("\n=== Warmstart Check ===")
    N     = game.time_horizon.N
    bound = hw + d_safe/1.25
    n_p   = game.n_players
 
    for i in 1:n_p
        traj   = ws.trajectories[i]
        max_py = maximum(abs.(traj.states[2, :]))
        max_u  = maximum(abs.(traj.controls))
        x0     = traj.states[:, 1]
        xmid   = traj.states[:, N÷2+1]
        xf     = traj.states[:, end]
        println("Player ", i,
                ":  x0=",   round.(x0,   sigdigits=4),
                "  xmid=",  round.(xmid, sigdigits=4),
                "  xf=",    round.(xf,   sigdigits=4))
        println("         max|py|=", round(max_py, sigdigits=4),
                "  bound=", round(bound, sigdigits=4),
                "  hallway ", max_py <= bound ? "✓ feasible" : "✗ INFEASIBLE",
                "  max|u|=", round(max_u, sigdigits=4))
    end
 
    # Collision check
    println("\nCollision check (d_safe=", d_safe, "):")
    any_collision = false
    for k in 1:N+1
        p1   = ws.trajectories[1].states[1:2, k]
        p2   = ws.trajectories[2].states[1:2, k]
        dist = norm(p1 .- p2)
        if dist < d_safe
            println("  COLLISION at k=", k, ": dist=", round(dist, sigdigits=4))
            any_collision = true
        end
    end
    any_collision || println("  No collisions")
end
 
# ── Diagnostic 5: inspect constraint violations along solution trajectory ──────
function check_constraints(sol, game; d_safe=0.2, hw=0.5)
    println("\n=== Constraint Violations Along Solution ===")
    N     = game.time_horizon.N
    bound = hw + d_safe/1.25
 
    coll_viols   = Float64[]
    hall_viols_1 = Float64[]
    hall_viols_2 = Float64[]
 
    for k in 1:N+1
        p1 = sol.trajectories[1].states[1:2, k]
        p2 = sol.trajectories[2].states[1:2, k]
        py1 = sol.trajectories[1].states[2, k]
        py2 = sol.trajectories[2].states[2, k]
 
        coll_g = d_safe^2 - sum(abs2, p1 .- p2)   # ≤ 0 satisfied
        h1_g   = abs(py1) - bound                   # ≤ 0 satisfied
        h2_g   = abs(py2) - bound
 
        push!(coll_viols,   max(0.0, coll_g))
        push!(hall_viols_1, max(0.0, h1_g))
        push!(hall_viols_2, max(0.0, h2_g))
    end
 
    println("Collision:  max_viol=", round(maximum(coll_viols),   sigdigits=4),
            "  mean=", round(mean(coll_viols),   sigdigits=4))
    println("Hallway P1: max_viol=", round(maximum(hall_viols_1), sigdigits=4),
            "  mean=", round(mean(hall_viols_1), sigdigits=4))
    println("Hallway P2: max_viol=", round(maximum(hall_viols_2), sigdigits=4),
            "  mean=", round(mean(hall_viols_2), sigdigits=4))
 
    # Per-timestep breakdown if violations present
    if maximum(coll_viols) > 1e-4 || maximum(hall_viols_1) > 1e-4
        println("\nPer-timestep violations > 0:")
        for k in 1:N+1
            cv  = coll_viols[k]
            h1v = hall_viols_1[k]
            h2v = hall_viols_2[k]
            if cv > 1e-4 || h1v > 1e-4 || h2v > 1e-4
                println("  k=", k,
                        "  coll=",  round(cv,  sigdigits=3),
                        "  hall1=", round(h1v, sigdigits=3),
                        "  hall2=", round(h2v, sigdigits=3))
            end
        end
    end
end
 
# ── Diagnostic 6: operating point trajectory summary ──────────────────────────
function check_trajectories(sol)
    println("\n=== Solution Trajectory Summary ===")
    for i in 1:length(sol.trajectories)
        traj = sol.trajectories[i]
        println("Player ", i, ":")
        println("  x0  = ", round.(traj.states[:, 1],   sigdigits=4))
        println("  xf  = ", round.(traj.states[:, end], sigdigits=4))
        println("  max |x| = ", round(maximum(abs.(traj.states)),   sigdigits=4))
        println("  max |u| = ", round(maximum(abs.(traj.controls)), sigdigits=4))
        println("  cost (stored) = ", round(traj.cost, sigdigits=5))
        any_nan = any(isnan, traj.states) || any(isnan, traj.controls)
        any_nan && println("  WARNING: NaN detected in trajectory!")
    end
end
 
# ── Run all diagnostics ───────────────────────────────────────────────────────
println("Running diagnostics on hallway problem solution...")
check_warmstart(ws, game)
check_dynamic_consistency(sol, game)
check_costs(sol, game)
check_convergence_trigger(sol, solver)
check_constraints(sol, game)
check_trajectories(sol)