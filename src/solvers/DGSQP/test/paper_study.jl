# ============================================================================
# Replication of the curved-track Monte Carlo study (Table I of Zhu & Borrelli, ICRA 2023)
#
# Runs DGSQP with the paper's settings on random collision-free initial conditions for
# θ ∈ {45°, 75°, 90°} and N ∈ {10, 15, 20, 25} and compares the number of converged trials with
# the DG-SQP row of Table I (out of 200 trials). Used by the scheduled `paper-study` job of
# .github/workflows/DGSQP.yml and runnable by hand:
#
#     DGSQP_STUDY_TRIALS=200 julia --project=. src/solvers/DGSQP/test/paper_study.jl
#
# The initial states are drawn from a different random stream than the paper's, so counts agree
# statistically, not trial by trial. The check allows a 5-percentage-point shortfall plus twice the
# sampling error of a binomial proportion.
# ============================================================================

using Test
using Printf
include(joinpath(@__DIR__, "..", "examples", "curved_track.jl"))

# DG-SQP "# Conv" row of Table I (out of 200 trials). Keys: (θ, N).
const PAPER_CONVERGED = Dict(
    (45, 10) => 200, (45, 15) => 199, (45, 20) => 199, (45, 25) => 185,
    (75, 10) => 200, (75, 15) => 199, (75, 20) => 195, (75, 25) => 169,
    (90, 10) => 200, (90, 15) => 197, (90, 20) => 193, (90, 25) => 172,
)

trials = parse(Int, get(ENV, "DGSQP_STUDY_TRIALS", "50"))
rows = run_curved_track_study(; trials = trials, seed = 1)

println()
@printf("%4s %4s | %9s %9s | %6s %6s %9s | %8s\n", "θ", "N", "paper", "ours", "fail", "maxit", "QPs/trial", "time [s]")
@testset "paper Table I (DG-SQP row), $trials trials per setting" begin
    for r in rows
        paper = PAPER_CONVERGED[(r.theta, r.N)] / 200
        ours = r.converged / r.trials
        @printf("%4d %4d | %8.1f%% %8.1f%% | %6d %6d %9.1f | %8.3f\n", r.theta, r.N, 100paper, 100ours,
                r.fail, r.max_iters, r.mean_qp_solves, r.mean_time)
        σ = sqrt(max(paper * (1 - paper), 0.0025) / r.trials)      # binomial sampling error
        @testset "θ = $(r.theta)°, N = $(r.N)" begin
            @test ours ≥ paper - 0.05 - 2σ
        end
    end
end
