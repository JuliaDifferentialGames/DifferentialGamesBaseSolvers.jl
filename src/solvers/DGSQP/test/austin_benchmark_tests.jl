# Tests for the Austin racing benchmark of the FALCON paper (examples/austin_hairpin.jl, examples/austin_benchmark.jl).
# They need DifferentialGamesBenchmarks (track data and vehicle model); they are skipped when it is not installed.
using Test
using LinearAlgebra
using DelimitedFiles

const _HAVE_BENCHMARKS = Base.find_package("DifferentialGamesBenchmarks") !== nothing

if _HAVE_BENCHMARKS
@testset "Austin racing benchmark (FALCON paper, Table: DG-SQP vs FALCON)" begin
    include(joinpath(@__DIR__, "..", "examples", "austin_benchmark.jl"))
    import DifferentialGamesBaseSolvers: _evaluate

    @testset "stored data reproduces the paper's DG-SQP row" begin
        ics = benchmark_ics()
        @test size(ics.q0_1) == (50, 7) && size(ics.q0_2) == (50, 7)
        # Table of the FALCON paper, DG-SQP: 72 % converged, 50.021 ± 19.269 s (31.696–124.621), 422.8 ± 138.0 iterations (276–947),
        # 81.471 s over all runs. These are the Python results stored with the initial conditions.
        c = ics.py_converged
        @test count(c) / 50 ≈ 0.72
        @test mean(ics.py_time[c]) ≈ 50.021 atol = 1e-3
        @test std(ics.py_time[c]) ≈ 19.269 atol = 1e-3
        @test minimum(ics.py_time[c]) ≈ 31.696 atol = 1e-3
        @test maximum(ics.py_time[c]) ≈ 124.621 atol = 1e-3
        @test mean(ics.py_iters[c]) ≈ 422.8 atol = 0.05
        @test (minimum(ics.py_iters[c]), maximum(ics.py_iters[c])) == (276, 947)
        @test mean(ics.py_time) ≈ 81.471 atol = 1e-3
        # FALCON (paper commit) on this machine, stored: 100 %, 0.873 ± 0.407 s, 85.9 ± 43.8 iterations
        F = readdlm(joinpath(@__DIR__, "..", "examples", "data", "falcon_results_paper_commit.csv"), ',', Float64; comments = true)
        @test all(F[:, 2] .== 1)
        @test mean(F[:, 4]) ≈ 0.873 atol = 0.05
        @test mean(F[:, 3]) ≈ 85.9 atol = 0.1
    end

    @testset "game: structured derivatives and feasible-looking initial states" begin
        ics = benchmark_ics()
        track = load_race_track(joinpath(pkgdir(DifferentialGamesBenchmarks), "data", "f1_austin_segment.csv"))
        q0 = [ics.q0_1[7, :], ics.q0_2[7, :]]
        game, w0, info = austin_hairpin_game(; N = 3, track = track, q0 = q0, v_ref = [1.5, 1.5], collision_aware = true)
        @test length(w0) == 18 && game.state_dims == [7, 7]
        og = OpenLoopGame(game)
        w = w0 .+ 0.01 .* sin.(1:length(w0))
        λ = abs.(cos.(1:length(og.constraints(w))))
        Qs, hs, Gs, gs = _evaluate(game, w, λ)
        Qo, ho, Go, go = _evaluate(og, w, λ)
        @test hs ≈ ho atol = 1e-10
        @test Gs ≈ Go atol = 1e-10
        @test gs ≈ go atol = 1e-10
        @test Qs ≈ Qo atol = 1e-8
        # the metrics are consistent with the plan
        m = hairpin_metrics(game, w0, info)
        @test m.min_separation > 0 && isfinite(m.track_clearance)
    end

    @testset "DG-SQP solves a hairpin instance (N = 12) and the solution satisfies the KKT conditions" begin
        track = hairpin_track()
        game, w0, info = austin_hairpin_game(; N = 12, track = track)
        res = dgsqp_solve(game, w0; solver = dgsqp_racing_options())
        @test res.converged
        @test res.cond.p_feas < 1e-4
        @test res.cond.stat < 1e-4
        m = hairpin_metrics(game, res.w, info)
        @test m.min_separation ≥ 0.5946 - 1e-4 && m.track_clearance ≥ -1e-4
    end
end
else
    @info "DifferentialGamesBenchmarks not found: skipping the Austin racing benchmark tests"
end
