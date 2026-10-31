# Timing and scaling report for the stochastic solvers.
#
#   julia --project=. benchmark/stochastic_scaling.jl [out.csv]
#
# Reports wall time of SFNELQ (perfect information and shared observation) and one
# CCiLQGames solve (default options, no Monte Carlo) versus horizon N and number of
# players P (each player: a 2D point with a scalar control). CCiLQGames time is a single
# post-compile run and its final status is recorded: a `max_iterations` row means the
# default iteration limits were not enough for that size, not that the time is comparable.
# SFNELQ is the minimum of three post-compile runs. Numbers are for regression tracking
# on one machine, not absolute claims; a per-phase breakdown is not implemented.
using LinearAlgebra, Printf
using DifferentialGamesBase, DifferentialGamesBaseSolvers

function lq_game(P, N; dt = 0.1)
    n = 2P
    A = Matrix(1.0I, n, n)
    B = [begin b = zeros(n, 1); b[2i-1] = dt; b end for i in 1:P]
    Q = [begin q = zeros(n, n); q[2i-1, 2i-1] = 4.0; q[2i, 2i] = 1.0; q end for i in 1:P]
    R = [fill(0.1, 1, 1) for _ in 1:P]
    qlin = [begin v = zeros(n); v[2i-1] = -4.0; v end for i in 1:P]
    LQGameProblem(A, B, Q, R, Q, zeros(n), N * dt; dt = dt, q = qlin)
end

function best_time(f; reps = 3)
    f()                                  # compile
    minimum(@elapsed(f()) for _ in 1:reps)
end

function run(out)
    rows = ["players,horizon,sfnelq_perfect_s,sfnelq_shared_s,ccilqgames_s,ccilqgames_status"]
    for P in (2, 4), N in (10, 20, 40)
        g = lq_game(P, N); n = 2P
        W = Matrix(0.01I, n, n)
        perfect = StochasticGameProblem(g, ProcessNoiseModel(W); initial_covariance = Matrix(0.01I, n, n))
        shared = StochasticGameProblem(g, ProcessNoiseModel(W); initial_covariance = Matrix(0.01I, n, n),
                                       information = linear_observation(Matrix(1.0I, n, n), Matrix(0.05I, n, n)))
        cc = ChanceConstraint(StateBounds(1; state_offset = 0, state_dim = 1, lower = [-5.0], upper = [0.6]), 0.1)
        chance = StochasticGameProblem(g, ProcessNoiseModel(W); chance_constraints = [cc])
        t1 = best_time(() -> solve(perfect, SFNELQ()))
        t2 = best_time(() -> solve(shared, SFNELQ()))
        status = Ref(:none)
        t3 = best_time(() -> (status[] = solve(chance, CCiLQGames(mc_samples = 0)).solver_info[:status]); reps = 1)
        push!(rows, @sprintf("%d,%d,%.5f,%.5f,%.5f,%s", P, N, t1, t2, t3, status[]))
        println(rows[end])
    end
    out === nothing || write(out, join(rows, "\n") * "\n")
end

run(isempty(ARGS) ? nothing : ARGS[1])
