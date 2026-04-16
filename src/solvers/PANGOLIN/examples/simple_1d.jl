using LinearAlgebra
using DifferentialGamesBase
using DifferentialGamesBaseSolvers
using Printf

# ============================================================================
# Simple 1D Constrained Game — AL Diagnostic
#
# Two players, each with 1D linear dynamics:
#   x_i(k+1) = x_i(k) + u_i(k) * dt
#
# Player 1: start at -1, goal at +1
# Player 2: start at +1, goal at -1
#
# Cost: J_i = Σ ρ_ctrl * u_i(k)² + ρ_goal * (x_i(N) - g_i)²
#
# Constraint: separation |x1 - x2| ≥ d_safe
#
# WHY SYMMETRY BREAKING IS REQUIRED:
#   This is a 1D crossing problem. The only feasible Nash equilibria are
#   asymmetric: one player yields to let the other pass. The symmetric
#   trajectory (both agents at equal speed) is infeasible — they collide
#   at x=0. If the solver initializes symmetrically, the AL gradient is also
#   symmetric and compresses both agents against the constraint equally,
#   rather than guiding one to yield. Fix: warmstart with P1 moving faster
#   than P2 so the solver finds the "P1 passes first" equilibrium.
# ============================================================================

println("=" ^ 60)
println("Simple 1D Constrained Game — AL Diagnostic")
println("=" ^ 60)

T    = Float64
N    = 20
dt   = T(0.1)
tf   = N * dt

ρ_ctrl = T(0.01)
ρ_goal = T(100.0)
d_safe = T(0.3)

x0_1 = T[-1.0]
x0_2 = T[ 1.0]
g1   = T[ 1.0]
g2   = T[-1.0]

dyn1 = (x, u, p, t) -> x .+ dt .* u
dyn2 = (x, u, p, t) -> x .+ dt .* u

stage_1 = NonlinearStageCost((x, u, p, t) -> ρ_ctrl * u[1]^2; is_separable=true)
stage_2 = NonlinearStageCost((x, u, p, t) -> ρ_ctrl * u[1]^2; is_separable=true)
term_1  = NonlinearTerminalCost((x, p) -> ρ_goal * (x[1] - g1[1])^2)
term_2  = NonlinearTerminalCost((x, p) -> ρ_goal * (x[1] - g2[1])^2)

p1 = PlayerSpec(1, 1, 1, x0_1, dyn1, PlayerObjective(1, stage_1, term_1))
p2 = PlayerSpec(2, 1, 1, x0_2, dyn2, PlayerObjective(2, stage_2, term_2))

# Smooth separation constraint: d_safe - sqrt((x1-x2)² + ε) ≤ 0
# sqrt avoids the zero ForwardDiff subgradient that abs() gives at contact.
smooth_abs_eps = T(1e-6)
sep_c = NonlinearConstraint(
    (x_joint, u, p, t) -> begin
        d = x_joint[1] - x_joint[2]
        [d_safe - sqrt(d^2 + smooth_abs_eps)]
    end,
    1; constraint_type = :inequality
)
shared = SharedConstraint(sep_c, [1, 2])
game   = PDGNEProblem([p1, p2], [shared], T(tf), T(dt))

println("\nProblem:")
println("  N=$N, dt=$dt, tf=$tf")
println("  Player 1: x0=$(x0_1[1]) → goal=$(g1[1])")
println("  Player 2: x0=$(x0_2[1]) → goal=$(g2[1])")
println("  d_safe=$d_safe, ρ_ctrl=$ρ_ctrl, ρ_goal=$ρ_goal")

# ── Test 1: Unconstrained ────────────────────────────────────────────────────
println("\n--- Test 1: Unconstrained (no constraint) ---")
game_unc  = PDGNEProblem([p1, p2], SharedConstraint[], T(tf), T(dt))
solver_unc = PANGOLIN(;
    max_iter             = 100,
    convergence_tol      = 1e-4,
    state_regularization = 1.0,
    constraint_opts      = ALOptions(ρ_init=1.0, ρ_max=1.0, φ=1.01,
                                     violation_tol=0.1, φ_threshold=0.75,
                                     max_consecutive_increases=3)
)
sol_unc = solve(game_unc, solver_unc; verbose=false)

xf1_unc = sol_unc.trajectories[1].states[1, end]
xf2_unc = sol_unc.trajectories[2].states[1, end]
min_sep = minimum(abs(sol_unc.trajectories[1].states[1,k] -
                      sol_unc.trajectories[2].states[1,k]) for k in 1:N+1)
println("  Converged: $(sol_unc.converged)  Iterations: $(sol_unc.iterations)")
println("  P1: xf=$(round(xf1_unc,sigdigits=4))  error=$(round(abs(xf1_unc-g1[1]),sigdigits=4))")
println("  P2: xf=$(round(xf2_unc,sigdigits=4))  error=$(round(abs(xf2_unc-g2[1]),sigdigits=4))")
println("  Min sep: $(round(min_sep,sigdigits=4))  ($(min_sep>=d_safe ? "✓" : "✗ COLLISION"))")

# ── Test 2: Constrained with asymmetric warmstart ────────────────────────────
println("\n--- Test 2: Constrained ---")

# Build asymmetric warmstart by solving an unconstrained game with asymmetric
# goal weights: P1 penalizes goal error 3× more heavily than P2, so P1 moves
# faster and crosses first. This produces a dynamically consistent trajectory
# (from a real solve) that lies in the basin of the "P1 passes first" Nash,
# without requiring us to manually construct Trajectory objects.
#
# WarmstartData only accepts a GameSolution — it cannot be built from raw
# trajectories — so we obtain it via a cheap unconstrained pre-solve.
println("  Computing asymmetric warmstart via pre-solve...")

stage_1_ws = NonlinearStageCost((x, u, p, t) -> ρ_ctrl * u[1]^2; is_separable=true)
stage_2_ws = NonlinearStageCost((x, u, p, t) -> ρ_ctrl * u[1]^2; is_separable=true)
term_1_ws  = NonlinearTerminalCost((x, p) -> 3.0 * ρ_goal * (x[1] - g1[1])^2)  # 3× weight → P1 faster
term_2_ws  = NonlinearTerminalCost((x, p) -> 1.0 * ρ_goal * (x[1] - g2[1])^2)  # 1× weight → P2 slower

p1_ws = PlayerSpec(1, 1, 1, x0_1, dyn1, PlayerObjective(1, stage_1_ws, term_1_ws))
p2_ws = PlayerSpec(2, 1, 1, x0_2, dyn2, PlayerObjective(2, stage_2_ws, term_2_ws))
game_ws = PDGNEProblem([p1_ws, p2_ws], SharedConstraint[], T(tf), T(dt))

solver_ws = PANGOLIN(;
    max_iter             = 50,
    convergence_tol      = 1e-3,
    state_regularization = 1.0,
    constraint_opts      = ALOptions(ρ_init=1.0, ρ_max=1.0, φ=1.01,
                                     violation_tol=0.1, φ_threshold=0.75,
                                     max_consecutive_increases=3)
)
sol_ws = solve(game_ws, solver_ws; verbose=false)
warmstart_data = WarmstartData(sol_ws)

ws_min_sep = minimum(
    abs(warmstart_data.trajectories[1].states[1,k] -
        warmstart_data.trajectories[2].states[1,k]) for k in 1:N+1)
println("  Warmstart min sep: $(round(ws_min_sep, sigdigits=4)) ($(ws_min_sep>=d_safe ? "feasible ✓" : "infeasible — solver will fix ✗"))")

solver_con = PANGOLIN(;
    max_iter             = 300,
    convergence_tol      = 1e-3,
    max_elwise_diff_step = 0.15,
    state_regularization = 1.0,
    control_limit        = 5.0,
    init_mode            = :warmstart,
    constraint_opts      = ALOptions(;
        ρ_init                    = 10.0,
        ρ_max                     = 1e4,
        φ                         = 1.5,
        violation_tol             = 1e-3,
        φ_threshold               = 0.75,
        max_consecutive_increases = 3
    )
)
sol_con = solve(game, solver_con; warmstart=warmstart_data, verbose=false)

xf1_con = sol_con.trajectories[1].states[1, end]
xf2_con = sol_con.trajectories[2].states[1, end]
min_sep_con = minimum(abs(sol_con.trajectories[1].states[1,k] -
                          sol_con.trajectories[2].states[1,k]) for k in 1:N+1)
println("  Converged: $(sol_con.converged)  Iterations: $(sol_con.iterations)")
println("  P1: xf=$(round(xf1_con,sigdigits=4))  error=$(round(abs(xf1_con-g1[1]),sigdigits=4))")
println("  P2: xf=$(round(xf2_con,sigdigits=4))  error=$(round(abs(xf2_con-g2[1]),sigdigits=4))")
println("  Min sep: $(round(min_sep_con,sigdigits=4))  ($(min_sep_con>=d_safe-1e-3 ? "✓" : "✗ VIOLATION"))")

println("\n  Per-timestep separation:")
for k in 1:N+1
    x1k  = sol_con.trajectories[1].states[1, k]
    x2k  = sol_con.trajectories[2].states[1, k]
    sep  = abs(x1k - x2k)
    viol = max(0.0, d_safe - sep)
    mark = viol > 1e-3 ? " ← VIOLATION" : ""
    @printf("    k=%2d: x1=%6.3f  x2=%6.3f  sep=%5.3f  viol=%6.4f%s\n",
            k, x1k, x2k, sep, viol, mark)
end

viols = sol_con.solver_info[:violation_history]
ρs    = sol_con.solver_info[:ρ_history]
costs = sol_con.solver_info[:iter_costs]
println("\n  AL history (every 10 iters):")
for k in 1:10:length(viols)
    @printf("    iter %3d: costs=[%.3f, %.3f]  viol=%.5f  ρ=%.2f\n",
            k, costs[k][1], costs[k][2], viols[k], ρs[k])
end

println("\n--- Test 3: Final AL multiplier state ---")
al_state = sol_con.solver_info[:al_state]
for j in 1:length(al_state.λ_shared_traj)
    λ     = al_state.λ_shared_traj[j]
    max_λ = maximum(λ)
    max_k = argmax(vec(λ))
    println("  Constraint $j: max(λ)=$(round(max_λ,sigdigits=4)) at k=$max_k")
    println("    λ profile: ", round.(vec(λ), sigdigits=3))
end
println("  Final ρ: $(al_state.ρ)")

println("\n" * "=" ^ 60)
println("Summary:")
println("  Unconstrained: min_sep=$(round(min_sep,sigdigits=4)) (expect collision)")
println("  Constrained:   min_sep=$(round(min_sep_con,sigdigits=4)) (expect ≥ $d_safe)")
println("  Goal error P1: $(round(abs(xf1_con-g1[1]),sigdigits=4))")
println("  Goal error P2: $(round(abs(xf2_con-g2[1]),sigdigits=4))")
println("=" ^ 60)