# ============================================================================
# DG-SQP on the F1 Austin racing benchmark of the FALCON paper
#
# Reproduces the setting of the paper's Table "Performance comparison of DG-SQP and FALCON for a 2-player
# F1 Austin racing scenario (N = 25, Δt = 0.1 s, M = 10, χ = 1e-3, ΔQ = 2e-2)": the dynamic-bicycle racing game on
# the sharp corner s ∈ [60, 80] m of the Austin track, solved from the 50 initial conditions of the FALCON
# repository (`dgsqp_benchmark.npz`, stored in `data/f1_austin_benchmark_ics.csv` together with the results of the
# original Python DG-SQP on them). DG-SQP runs with the settings of the reference benchmark script
# (`benchmark_dgsqp.py`): racing-paper step acceptance, ε₀ = 100 decaying by 0.95, checkpoints every 10 iterations,
# merit memory 10, line search with 20 steps, feasibility tolerance 1e-4 and stationarity/complementarity tolerance 2e-2.
#
#     include("austin_benchmark.jl")
#     rows = austin_dgsqp_benchmark(; samples = 50)
#     summarize_benchmark(rows)
# ============================================================================

include("austin_hairpin.jl")
using DelimitedFiles
import Statistics
const mean = Statistics.mean
const std = Statistics.std

"""
    benchmark_ics() -> NamedTuple

The 50 initial conditions (`q0_1`, `q0_2`, 7 states each: `[vₓ, v_y, ψ̇, x, y, ψ, s]`) and the Python DG-SQP results
(`py_time`, `py_iters`, `py_converged`) of the FALCON benchmark.
"""
function benchmark_ics()
    A = readdlm(joinpath(@__DIR__, "data", "f1_austin_benchmark_ics.csv"), ',', Float64; comments = true)
    return (q0_1 = A[:, 1:7], q0_2 = A[:, 8:14], py_time = A[:, 15], py_iters = A[:, 16], py_converged = A[:, 17] .> 0.5)
end

"""
    dgsqp_benchmark_options(; kwargs...) -> DGSQP

Settings of `benchmark_dgsqp.py` / `globals.py` of the reference implementation.
"""
dgsqp_benchmark_options(; kwargs...) =
    DGSQP(; step_acceptance = :nms, nms_frequency = 10, nms_memory_size = 10, line_search_iters = 20, sqp_iters = 500,
            p_tol = 1e-4, d_tol = 2e-2, reg = 1e2, reg_decay = 0.95, delta_decay = 0.99, merit_decrease = 0.01,
            tau = 0.5, merit_function = :stat_l1, divergence_tol = 1e10, kwargs...)

"""
    austin_dgsqp_benchmark(; samples = 50, N = 25, q_l = 20.0, solver = dgsqp_benchmark_options(), verbose = true)

Run DG-SQP on the first `samples` initial conditions. `warm_starts[i] = (U1, U2)` (`3×N` control matrices per car) supplies
the warm start of sample `i`, e.g. the single-agent plans computed by FALCON's benchmark script, as in the paper;
without it a closed-loop controller plan is used. Returns one `NamedTuple` per sample with the solve time (warm start
excluded), the number of SQP iterations, convergence and the final KKT residuals.
"""
function austin_dgsqp_benchmark(; samples::Int = 50, N::Int = 25, q_l::Float64 = 20.0,
                                  solver::DGSQP = dgsqp_benchmark_options(), verbose::Bool = true, warm_starts = nothing)
    ics = benchmark_ics()
    track = load_race_track(joinpath(pkgdir(DifferentialGamesBenchmarks), "data", "f1_austin_segment.csv"))
    rows = NamedTuple[]
    for i in 1:samples
        q0 = [ics.q0_1[i, :], ics.q0_2[i, :]]
        game, w0, info = warm_starts === nothing ?
            austin_hairpin_game(; N = N, q_l = q_l, track = track, q0 = q0,
                                  v_ref = [max(q0[1][1], 1.0), max(q0[2][1], 1.0)], collision_aware = true) :
            austin_hairpin_game(; N = N, q_l = q_l, track = track, q0 = q0, warm_start = warm_starts[i])
        res = dgsqp_solve(game, w0; solver = solver)
        m = hairpin_metrics(game, res.w, info)
        push!(rows, (sample = i, converged = res.converged, status = res.status, iterations = res.iterations,
                     time = res.time, cond = res.cond, min_separation = m.min_separation,
                     track_clearance = m.track_clearance,
                     py_converged = ics.py_converged[i], py_iters = ics.py_iters[i], py_time = ics.py_time[i]))
        verbose && @printf("[%2d/%d] %-8s iters %4d  time %7.2f s  (python: %s, %4d iters, %6.1f s)  KKT %.1e/%.1e/%.1e\n",
                           i, samples, res.converged ? "CONVERGED" : string(res.status), res.iterations, res.time,
                           ics.py_converged[i] ? "conv" : "fail", ics.py_iters[i], ics.py_time[i],
                           res.cond.p_feas, res.cond.comp, res.cond.stat)
        flush(stdout)
    end
    return rows
end

_ms(x) = isempty(x) ? "—" : @sprintf("%.3f ± %.3f (%.3f–%.3f)", mean(x), std(x), minimum(x), maximum(x))

"""
    summarize_benchmark(rows)

Print the table of the FALCON paper for DG-SQP: convergence rate, time and iterations of the converged runs
(mean ± std (min–max)), and the mean time over all runs; next to the same statistics of the Python DG-SQP on the same
initial conditions.
"""
function summarize_benchmark(rows)
    conv = filter(r -> r.converged, rows)
    pyc = filter(r -> r.py_converged, rows)
    println("                       conv. rate [%]   time converged [s]               iterations (converged)         time all [s]")
    @printf("Julia DG-SQP            %6.1f          %-32s %-30s %.3f\n", 100 * length(conv) / length(rows), _ms([r.time for r in conv]),
            _ms([r.iterations for r in conv]), mean(r.time for r in rows))
    @printf("Python DG-SQP (stored)  %6.1f          %-32s %-30s %.3f\n", 100 * length(pyc) / length(rows), _ms([r.py_time for r in pyc]),
            _ms([r.py_iters for r in pyc]), mean(r.py_time for r in rows))
    both = count(r -> r.converged && r.py_converged, rows)
    println("converged in both: $both, Julia only: $(count(r -> r.converged && !r.py_converged, rows)), ",
            "Python only: $(count(r -> !r.converged && r.py_converged, rows)), neither: $(count(r -> !r.converged && !r.py_converged, rows))")
    return nothing
end

"""
    austin_dgsqp_benchmark_tuned(; kwargs...)

[`austin_dgsqp_benchmark`](@ref) with [`dgsqp_tuned_options`](@ref).
"""
austin_dgsqp_benchmark_tuned(; kwargs...) = austin_dgsqp_benchmark(; solver = dgsqp_tuned_options(), kwargs...)
