# ============================================================================
# Export of the FALCON side of the FALCON-vs-DGSQP comparison (run in the FALCON environment)
#
# Runs FALCON's own benchmark protocol (`f1_austin_benchmark.jl` of the FALCON repository) on the 50 initial conditions
# of the paper and writes, per initial condition, FALCON's single-agent warm-start controls (the warm start both methods use in
# the paper) and FALCON's result (converged, iterations, time, χ, ΔQ) to `falcon_ws_and_results.csv`. The DGSQP side then
# reads that file: `austin_benchmark.jl`.
#
#   FALCON_EX=<FALCON checkout>/examples julia --project=<FALCON environment> falcon_export.jl 50 <output directory>
#
# NOTE: the paper numbers are from FALCON commit 6db00e7 ("paper results"). Later working trees need FALCON_CONFIG=tuned (below).
# Columns of the CSV: index, converged, iterations, time [s], χ, ΔQ, U¹ (3×25, column-major), U² (3×25).
# ============================================================================
# Runs FALCON's own benchmark protocol on the 50 ICs and exports (a) its single-agent warm-start controls and (b) its results.
ENV["GKSwstype"] = "nul"
const FALCON_EX = get(ENV, "FALCON_EX", "/home/bennet/Dropbox/Documents/Research/methods/FALCON/examples")
include(joinpath(FALCON_EX, "f1_austin_racing.jl"))
using Statistics, Printf, DelimitedFiles
data = load_npz_float64(joinpath(FALCON_EX, "dgsqp_benchmark.npz"), ["q0_1", "q0_2", "dgsqp_times", "dgsqp_iters", "dgsqp_converged"])
n = parse(Int, ARGS[1]); out = ARGS[2]
track = load_smooth_track(joinpath(FALCON_EX, "f1_austin_tenth_scale.npz"), 60.0, 80.0)
# FALCON_CONFIG selects the solver parameters:
#   "paper"  (default) the settings of the paper's benchmark script (valid at the paper commit 6db00e7);
#   "tuned"  ω_max = 1e4 and r_min = 1e-5 (the penalty cap of 1e3 and the trust-region floor of 1e-3 stop a few initial
#            conditions short of feasibility), plus, on working trees that have the mandatory ζ certificate, the loosened
#            tolerances ζ_cert_abs_tol = 2e-3, ζ_cert_rel_tol = 1e-2, ζ_cert_polish_ΔJ_tol = 1. 50 of 50 converge.
const _CONFIG = get(ENV, "FALCON_CONFIG", "paper")
const _TUNED = _CONFIG == "tuned" ?
    merge((ω_max = 1e4, r_min = 1e-5),
          hasfield(FALCONSolver, :ζ_cert_abs_tol) ? (ζ_cert_abs_tol = 2e-3, ζ_cert_rel_tol = 1e-2, ζ_cert_polish_ΔJ_tol = 1.0) : NamedTuple()) :
    NamedTuple()
solver = FALCONSolver(; r0=0.5, r_max=1.0, α_2=1.5, η=100.0, w0=10.0, β=2.0, γ_δ=0.9, ω_max=1e3, ϵ_opt=2e-2, ϵ_feas=1e-4,
    max_iterations=400, linear_solver=:qdldl, r_min=1e-3, socp_tol_loose=1e-5, socp_tol_tight=1e-5, socp_tol_switch=1e-2,
    infeasibility_preprocessing=true, preprocess_iters=5, _TUNED...)
let q1 = data["q0_1"][1, :], q2 = data["q0_2"][1, :]
    g, w = build_racing_game(track, q1, q2, 25, 0.1, 10; ws_verbose=false); solve(g, solver; warmstart=w, verbose=false, check_compatibility=false)
end
rows = Vector{Float64}[]
for i in 1:n
    q1 = data["q0_1"][i, :]; q2 = data["q0_2"][i, :]
    game, ws = build_racing_game(track, q1, q2, 25, 0.1, 10; ws_verbose=false)
    t0 = time(); sol = solve(game, solver; warmstart=ws, verbose=false, check_compatibility=false); ts = time() - t0
    ni = get(sol.solver_info, :scp_iterations, 0); χ = get(sol.solver_info, :infeasibility, NaN); ΔQ = get(sol.solver_info, :ΔQ_max, NaN)
    U1 = ws.trajectories[1].controls; U2 = ws.trajectories[2].controls   # 3×25 each
    push!(rows, vcat(Float64(i), sol.converged ? 1.0 : 0.0, Float64(ni), ts, χ, ΔQ, vec(U1), vec(U2)))
    @printf("[%2d/%d] %s iters=%3d t=%.3f s  chi=%.1e dQ=%.1e\n", i, n, sol.converged ? "CONV" : "FAIL", ni, ts, χ, ΔQ); flush(stdout)
end
writedlm(joinpath(out, "falcon_ws_and_results.csv"), reduce(hcat, rows)', ',')
