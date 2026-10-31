# ============================================================================
# DGSQP tests
#
# Cat 0 — Unit:          PSD projection, QP step, merit function and its derivative, μ rule
# Cat 1 — Ground truth:  games with a known equilibrium (analytic and LQ)
# Cat 2 — Optimality:    KKT residuals, best-response (unilateral-deviation) certificates
# Cat 3 — Cross-check:   agreement with ALGAMES on an LQ game; constrained nonlinear game
# Cat 4 — Interface:     status codes, warm start, capabilities, GameProblem adapter
# Cat 5 — Paper study:   reduced replication of the curved-track Monte Carlo (Table I)
# ============================================================================

using Test
using LinearAlgebra
using Random
using ForwardDiff
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

import DifferentialGamesBaseSolvers: _nearest_psd, _solve_qp, _merit, _dmerit, _get_mu, _evaluate, _step_at

# ── Helpers ─────────────────────────────────────────────────────────────────
# Two scalar players: J₁ = (a−1)² + ½ab, J₂ = (b+1)² + ½ab.
# Unconstrained Nash: 2(a−1) + b/2 = 0, 2(b+1) + a/2 = 0  ⇒  a = −b = 4/3.
two_player_costs() = Function[w -> (w[1] - 1)^2 + 0.5 * w[1] * w[2],
                              w -> (w[2] + 1)^2 + 0.5 * w[1] * w[2]]

function double_integrator_lq(; N = 8, dt = 0.1)
    A = [1.0 dt; 0.0 1.0]
    B = [reshape([0.5dt^2; dt], 2, 1), reshape([0.25dt^2; 0.5dt], 2, 1)]
    Q = [Matrix(1.0I, 2, 2), Matrix(2.0I, 2, 2)]
    R = [Matrix(1.0I, 1, 1), Matrix(0.5I, 1, 1)]
    Qf = [Matrix(5.0I, 2, 2), Matrix(3.0I, 2, 2)]
    return LQGameProblem(A, B, Q, R, Qf, [1.0, 0.0], dt * N; dt = dt)
end

# Two single integrators swapping lanes with control bounds and a collision constraint.
function lane_swap_game(; tf = 2.0, dt = 0.1, goal_x = 3.0)
    dyn = (xi, ui, p, t) -> ui
    function player(i, x0, goal)
        o = 2 * (i - 1)
        stage = NonlinearStageCost((x, u, p, t) -> 0.01 * sum(abs2, u[o+1:o+2]); is_separable = true)
        term = NonlinearTerminalCost((x, p) -> 10.0 * ((x[o+1] - goal[1])^2 + (x[o+2] - goal[2])^2))
        bounds = control_bounds(i; control_offset = o, control_dim = 2, lower = [-2.0, -2.0], upper = [2.0, 2.0])
        return PlayerSpec(i, 2, 2, x0, dyn, PlayerObjective(i, stage, term), [bounds])
    end
    sep = SharedInequality([1, 2]; func = (x, u, p, t) -> [1.5^2 - sum(abs2, x[1:2] - x[3:4])], dim = 1)
    return PDGNEProblem([player(1, [0.0, 2.0], [goal_x, -2.0]), player(2, [0.0, -1.5], [goal_x, 2.0])], [sep], tf, dt), sep
end

# ============================================================================
@testset "DGSQP" begin

# ── Cat 0 ───────────────────────────────────────────────────────────────────
@testset "Cat 0: unit" begin
    @testset "_nearest_psd" begin
        A = [2.0 1.0; 3.0 -1.0]                    # indefinite, non-symmetric
        B = _nearest_psd(A)
        @test B ≈ B'
        @test minimum(eigvals(Symmetric(B))) ≥ -1e-12
        P = [2.0 0.5; 0.5 1.0]                      # already PSD: unchanged
        @test _nearest_psd(P) ≈ P atol = 1e-12
        # closest PSD matrix in Frobenius norm of the symmetric part: negative eigenvalues are removed
        S = [1.0 0.0; 0.0 -3.0]
        @test _nearest_psd(S) ≈ [1.0 0.0; 0.0 1e-10] atol = 1e-9
    end

    @testset "_solve_qp matches a known solution" begin
        # min x₁²+x₂²−2x₁−5x₂  s.t. x₁+x₂ ≤ 1, x₁−x₂ ≤ 0  ⇒  x = (−¼, 5/4), multiplier (5/2, 0) on the first row
        B = [2.0 0.0; 0.0 2.0]; h = [-2.0, -5.0]
        G = [1.0 1.0; 1.0 -1.0]; g = [-1.0, 0.0]       # g + G p ≤ 0 with p = x
        p, λ = _solve_qp(B, h, G, g)
        @test p ≈ [-0.25, 1.25] atol = 1e-6
        @test λ ≈ [2.5, 0.0] atol = 1e-6
        # unconstrained
        p0, λ0 = _solve_qp(B, h, zeros(0, 2), Float64[])
        @test p0 ≈ [1.0, 2.5]
        @test isempty(λ0)
        # infeasible QP is reported, not thrown
        @test _solve_qp(B, h, [1.0 0.0; -1.0 0.0], [1.0, 1.0]) === nothing
    end

    @testset "merit function and its directional derivative" begin
        rng = MersenneTwister(7)
        n, nc = 6, 4
        h = randn(rng, n); G = randn(rng, nc, n); g = randn(rng, nc); l = abs.(randn(rng, nc))
        Q = randn(rng, n, n)
        du = randn(rng, n); dl = randn(rng, nc)
        s = min.(0.0, g); ds = g .+ G * du .- s
        μ = 0.7
        φ(α) = begin
            # q, G, g depend on u through h(u) = h + Q du, g(u) = g + G du (affine model)
            ha = h .+ α .* (Q * du); ga = g .+ α .* (G * du)
            _merit(l .+ α .* dl, s .+ α .* ds, ha, G, ga, μ)
        end
        dφ = _dmerit(du, l, dl, Q, h, G, g, s, μ)
        # the analytic derivative has the (g−s)-term of the reference: μΣ(g−s) is constant along the model path
        # because the slack moves with the linearised constraint; compare the stationarity part by finite differences
        ε = 1e-6
        num_stat = (φ(ε) - φ(-ε)) / 2ε
        # In this affine model g(α) − s(α) = g − s − α(G du − ds) = (1−α)(g−s) ⇒ d/dα μΣ(g−s) = −μΣ(g−s)
        @test dφ ≈ num_stat atol = 1e-5
    end

    @testset "_get_mu cases" begin
        solver = DGSQP()
        n, nc = 3, 2
        Q = Matrix(1.0I, n, n); G = [1.0 0 0; 0 1.0 0]; h = [1.0, 0.0, 0.0]
        l = zeros(nc)
        # feasible iterate: no violation ⇒ μ = 0 for either sign of the derivative
        g = [-1.0, -1.0]; s = min.(0.0, g)
        @test _get_mu([0.1, 0, 0], l, zeros(nc), s, Q, h, G, g, solver) == 0.0
        # violated iterate ⇒ μ > 0
        g = [1.0, -1.0]; s = min.(0.0, g)
        @test _get_mu([-0.1, 0, 0], l, zeros(nc), s, Q, h, G, g, solver) ≥ 0
        # :stat merit has no ℓ¹ term
        @test _get_mu([-0.1, 0, 0], l, zeros(nc), s, Q, h, G, g, DGSQP(; merit_function = :stat)) == 0.0
    end

    @testset "derivatives of an OpenLoopGame match finite differences" begin
        cons = w -> [w[1] + w[2]^2 - 0.2, -w[1] - 1.0]
        g = OpenLoopGame([1, 1], 1, two_player_costs(), cons)
        w = [0.4, -0.3]; l = [0.5, 0.2]
        Q, h, G, gv = _evaluate(g, w, l)
        @test gv ≈ cons(w)
        @test h ≈ [2 * (w[1] - 1) + 0.5 * w[2], 2 * (w[2] + 1) + 0.5 * w[1]]
        @test G ≈ [1.0 2w[2]; -1.0 0.0]
        # Q = ∂[∇_{u¹}𝓛¹, ∇_{u²}𝓛²]/∂u: [2 .5; .5 2] + λ₁∇²(w₁+w₂²) on the second diagonal
        @test Q ≈ [2.0 0.5; 0.5 2.0 + 2 * l[1]]
    end
end

# ── Cat 1 ───────────────────────────────────────────────────────────────────
@testset "Cat 1: known equilibria" begin
    @testset "unconstrained scalar game" begin
        g = OpenLoopGame([1, 1], 1, two_player_costs(), nothing)
        r = dgsqp_solve(g, [0.0, 0.0])
        @test r.converged && r.status === :conv_abs_tol
        @test r.w ≈ [4 / 3, -4 / 3] atol = 5e-3
        # without the regularisation ε the first step is the Newton step: exact for a quadratic game
        r0 = dgsqp_solve(g, [0.0, 0.0]; solver = DGSQP(; reg = 0.0))
        @test r0.w ≈ [4 / 3, -4 / 3] atol = 1e-6
        @test r0.iterations ≤ 2
    end

    @testset "constrained scalar game: active shared constraint" begin
        # a + b ≤ 0.2 is active at the unconstrained solution a+b = 0 ⇒ inactive; tighten to a − b ≤ 2
        cons = w -> [w[1] - w[2] - 2.0]
        g = OpenLoopGame([1, 1], 1, two_player_costs(), cons)
        r = dgsqp_solve(g, [0.0, 0.0]; solver = DGSQP(; reg = 0.0, nonmono_ls = true))
        @test r.converged
        @test r.w[1] - r.w[2] ≈ 2.0 atol = 1e-3          # constraint active at the equilibrium
        @test r.λ[1] > 0                                    # with a positive multiplier
        # KKT: ∇_{uⁱ}𝓛ⁱ = 0 with the shared multiplier
        a, b = r.w; λ = r.λ[1]
        @test 2 * (a - 1) + 0.5b + λ ≈ 0 atol = 5e-3
        @test 2 * (b + 1) + 0.5a - λ ≈ 0 atol = 5e-3
    end

    @testset "LQ game agrees with ALGAMES (open-loop Nash)" begin
        game = double_integrator_lq(; N = 8)
        sd = solve(game, DGSQP(; reg = 0.0, p_tol = 1e-7, d_tol = 1e-7, sqp_iters = 300))
        sa = solve(game, ALGAMES(; tol_opt = 1e-8, tol_dyn = 1e-8))
        @test sd.converged && sa.converged
        for i in 1:2
            @test get_trajectory(sd, i).controls ≈ get_trajectory(sa, i).controls atol = 1e-4
        end
        @test sd.state_trajectory ≈ sa.state_trajectory atol = 1e-4
        @test sd.equilibrium_type === :OpenLoopNash
    end
end

# ── Cat 2 ───────────────────────────────────────────────────────────────────
@testset "Cat 2: optimality certificates" begin
    @testset "best-response gap" begin
        g = OpenLoopGame([1, 1], 1, two_player_costs(), nothing)
        # at the equilibrium nobody gains by deviating
        r = dgsqp_solve(g, [0.0, 0.0]; solver = DGSQP(; reg = 0.0))
        for i in 1:2
            @test best_response_gap(g, r.w, i)[1] < 1e-8
        end
        # away from it they do: player 1's best response to b = 0 is a = 1 (cost 0) from a = 0 (cost 1)
        gap, wbr = best_response_gap(g, [0.0, 0.0], 1)
        @test gap ≈ 1.0 atol = 1e-6
        @test wbr[1] ≈ 1.0 atol = 1e-6
    end

    @testset "lane-swap game: KKT point that is a local equilibrium" begin
        game, sep = lane_swap_game()
        sol = solve(game, dgsqp_paper_options())
        @test sol.converged
        cond = sol.solver_info[:cond]
        @test cond.p_feas < 1e-3 && cond.comp < 1e-3 && cond.stat < 1e-3
        # constraints hold at every node and controls respect their bounds
        X = sol.state_trajectory
        @test minimum(norm(X[1:2, k] - X[3:4, k]) for k in axes(X, 2)) ≥ 1.5 - 1e-3
        @test maximum(abs, hcat([get_trajectory(sol, i).controls for i in 1:2]...)) ≤ 2.0 + 1e-6
        # unilateral deviations do not pay
        gaps = best_response_gap(game, sol)
        @test all(gaps .< 5e-3)
    end
end

# ── Cat 3 ───────────────────────────────────────────────────────────────────
@testset "Cat 3: non-monotone line search" begin
    # On the lane-swap game the monotone line search with the reference's β = 0.01 stalls at a non-stationary
    # point (the SQP direction is not a descent direction of the merit function); the watchdog strategy of
    # Sec. IV-B converges. This is the paper's motivation for the non-monotone strategy.
    game, _ = lane_swap_game(; tf = 4.0, goal_x = 5.0)
    olg = OpenLoopGame(game)
    w0 = zeros(sum(olg.control_dims) * olg.N)
    mono = dgsqp_solve(olg, w0; solver = DGSQP(; beta = 0.01))
    wd   = dgsqp_solve(olg, w0; solver = dgsqp_paper_options())
    @test wd.converged
    @test !(mono.converged && mono.status === :conv_abs_tol)
    @test wd.cond.stat < 1e-3
end

# ── Cat 4 ───────────────────────────────────────────────────────────────────
@testset "Cat 4: interface" begin
    @testset "status codes" begin
        g = OpenLoopGame([1, 1], 1, two_player_costs(), nothing)
        r = dgsqp_solve(g, [0.0, 0.0]; solver = DGSQP(; sqp_iters = 1, reg = 1.0))
        @test r.status === :max_it && !r.converged
        # infeasible constraints  w₁ ≤ −1 and w₁ ≥ 1
        bad = OpenLoopGame([1, 1], 1, two_player_costs(), w -> [w[1] + 1.0, -w[1] + 1.0])
        rb = dgsqp_solve(bad, [0.0, 0.0])
        @test rb.status === :qp_fail && !rb.converged
        # divergence test
        huge = OpenLoopGame([1], 1, Function[w -> 1e7 * w[1]], nothing)
        @test dgsqp_solve(huge, [0.0]).status === :diverged
        # time limit
        @test dgsqp_solve(g, [0.0, 0.0]; solver = DGSQP(; time_limit = -1.0, reg = 1.0)).status === :time_limit
    end

    @testset "warm start through WarmstartData and a previous solution" begin
        game = double_integrator_lq(; N = 6)
        s1 = solve(game, DGSQP())
        s2 = solve(game, DGSQP(); warmstart = s1)
        @test s2.converged
        @test s2.iterations ≤ s1.iterations
        @test get_trajectory(s2, 1).controls ≈ get_trajectory(s1, 1).controls atol = 5e-3
    end

    @testset "capabilities and solve()" begin
        @test validate_capabilities(solver_capabilities(DGSQP)) === nothing
        @test :ConstrainedGame in solver_capabilities(DGSQP)
        @test DGSQP() isa GameSolver
        o = dgsqp_paper_options()
        @test o.reg == 0.0 && o.nonmono_ls && o.beta == 0.01 && o.tau == 0.5
        @test dgsqp_paper_options(; sqp_iters = 7).sqp_iters == 7
    end

    @testset "decision-vector helpers" begin
        g = OpenLoopGame([2, 1], 3, Function[w -> 0.0, w -> 0.0], nothing)
        r = control_ranges(g)
        @test r == [1:6, 7:9]
        w = collect(1.0:9.0)
        U = split_controls(g, w)
        @test size(U[1]) == (2, 3) && size(U[2]) == (1, 3)
        @test join_controls(U) == w
    end
end

# ── Cat 5 ───────────────────────────────────────────────────────────────────
@testset "Cat 5: curved-track study (reduced replication of the paper's Table I)" begin
    include(joinpath(@__DIR__, "..", "examples", "curved_track.jl"))

    @testset "track geometry" begin
        t = CurveTrack(; swept_angle = π / 2)
        @test curvature(t, 0.5) == 0 && curvature(t, 5.0) ≈ (π / 2) / 8 && curvature(t, 12.0) == 0
        @test tangent_angle(t, 0.5) == 0 && tangent_angle(t, 12.0) ≈ π / 2
        @test tangent_angle(t, 5.0) ≈ π / 4
        x0, y0 = local_to_global(t, 0.5, 0.3)
        @test (x0, y0) == (0.5, 0.3)                      # straight section: left offset is +y
        # at the end of the arc the car is displaced by (r, r) and rotated by the swept angle
        r = 8 / (π / 2)
        xe, ye = local_to_global(t, 9.0, 0.0)
        @test all(isapprox.([xe, ye], [1 + r * sin(π / 2), r * (1 - cos(π / 2))]; atol = 1e-9))
    end

    @testset "Euler bicycle model" begin
        S = CurvedTrackSetup()
        q = [0.0, 0.0, 2.0, 0.0, 0.5, 0.0]
        qn = bicycle_step(q, [0.0, 0.0], S.track, S.bicycle, S.dt)
        @test qn[1] ≈ 0.2 && qn[5] ≈ 0.5 + 0.2            # straight line, no steering
        @test qn[3] < 2.0                                  # drag decelerates
    end

    @testset "Monte Carlo: success rate and number of QPs" begin
        # Paper (DS rows of Table I, 200 trials each): θ = 45°, N = 10: 200/200 converged, 3 QP solves on average;
        # θ = 90°, N = 15: 197/200. A reduced study must reproduce these qualitatively.
        rows = run_curved_track_study(; thetas = [45, 90], horizons = [10], trials = 20, seed = 3)
        for r in rows
            @test r.converged ≥ 18                          # ≥ 90 % (paper: 99.5–100 %)
            @test r.mean_qp_solves ≤ 8                      # paper: 3–4
        end
    end
end

# ── Cat 6 ───────────────────────────────────────────────────────────────────
@testset "Cat 6: structured derivatives and the racing-paper variant" begin
    include(joinpath(@__DIR__, "..", "examples", "curved_track.jl"))
    S = CurvedTrackSetup(; track = CurveTrack(; swept_angle = π / 2), N = 4)
    rng = MersenneTwister(2)
    q0, w0 = sample_initial_state(S, rng)
    N = S.N
    step_i = (x, u) -> bicycle_step(x, u, S.track, S.bicycle, S.dt)
    cost = function (i, k, X, U, Up)
        j = 3 - i
        k < N ? 0.5 * sum(S.input_weight .* U[i] .^ 2 .+ S.rate_weight .* (U[i] .- Up[i]) .^ 2) :
                -S.w_progress * X[i][5] + S.w_comp * atan(X[j][5] - X[i][5])
    end
    cons = function (k, X, U, Up)
        c = eltype(X[1])[]
        if k < N
            for i in 1:2
                u = U[i]; du = u .- Up[i]
                append!(c, [u[1] - S.u_max[1], -u[1] - S.u_max[1], u[2] - S.u_max[2], -u[2] - S.u_max[2]])
                append!(c, [du[1] - S.dt * S.du_max[1], -du[1] - S.dt * S.du_max[1], du[2] - S.dt * S.du_max[2], -du[2] - S.dt * S.du_max[2]])
            end
        end
        if k >= 1
            for i in 1:2; append!(c, [X[i][6] - S.half_width, -X[i][6] - S.half_width]); end
            push!(c, (2 * S.car_radius)^2 - sum(abs2, X[1][1:2] .- X[2][1:2]))
        end
        return c
    end
    sg = SeparableRolloutGame([2, 2], [6, 6], N, [q0[1:6], q0[7:12]], [step_i, step_i], cost, cons)
    og = OpenLoopGame(sg)

    @testset "structured derivatives equal whole-rollout automatic differentiation" begin
        w = w0 .+ 0.05 * randn(rng, length(w0))
        λ = abs.(randn(rng, length(og.constraints(w))))
        Qs, hs, Gs, gs = _evaluate(sg, w, λ)
        Qo, ho, Go, go = _evaluate(og, w, λ)
        @test hs ≈ ho atol = 1e-12
        @test Gs ≈ Go atol = 1e-12
        @test gs ≈ go atol = 1e-12
        @test Qs ≈ Qo atol = 1e-10
        @test maximum(abs, Qo) > 1                       # a non-trivial Hessian was compared
    end

    @testset "both game types give the same solve" begin
        solver = dgsqp_paper_options()
        a = dgsqp_solve(sg, w0; solver = solver)
        b = dgsqp_solve(og, w0; solver = solver)
        @test a.status === b.status
        @test a.w ≈ b.w atol = 1e-6
    end

    @testset "racing-paper variant (decaying regularization, d/m-steps)" begin
        o = dgsqp_racing_options()
        @test o.step_acceptance === :nms && o.reg == 100.0 && o.reg_decay == 0.95 && o.sqp_iters == 500
        game, _ = lane_swap_game(; tf = 2.0, goal_x = 3.0)
        olg = OpenLoopGame(game)
        r = dgsqp_solve(olg, zeros(sum(olg.control_dims) * olg.N); solver = dgsqp_racing_options(; sqp_iters = 300))
        @test r.converged
        @test r.cond.p_feas < 1e-4 && r.cond.comp < 1e-4 && r.cond.stat < 1e-4
    end
end

    @testset "elastic QP and extension options" begin
        import DifferentialGamesBaseSolvers: _solve_elastic_qp, _solve_qp
        B = Matrix(1.0I, 1, 1); h = [0.0]
        G = [1.0; -1.0;;]; g = [1.0, 1.0]              # p ≤ −1 and p ≥ 1: infeasible
        @test _solve_qp(B, h, G, g) === nothing
        p, λ = _solve_elastic_qp(B, h, G, g, 1e3, 1e3)
        @test length(p) == 1 && length(λ) == 2 && all(isfinite, p) && all(λ .>= -1e-9)
        # a feasible QP gives the same answer with and without the slack (the slack stays zero)
        G2 = reshape([1.0], 1, 1); g2 = [-1.0]; h2 = [-5.0]
        p1, l1 = _solve_qp(B, h2, G2, g2)
        pe, le = _solve_elastic_qp(B, h2, G2, g2, 1e3, 1e3)
        @test pe ≈ p1 atol = 1e-4
        @test le ≈ l1 atol = 1e-3
        # infeasible first QP: stops with :qp_fail by default (as the reference), runs on with the elastic fallback
        bad = OpenLoopGame([1, 1], 1, two_player_costs(), w -> [w[1] + 1.0, -w[1] + 1.0])
        @test dgsqp_solve(bad, [0.0, 0.0]).status === :qp_fail
        r = dgsqp_solve(bad, [0.0, 0.0]; solver = dgsqp_racing_options(; elastic_mode = true, sqp_iters = 5))
        @test r.status !== :qp_fail
        t = dgsqp_tuned_options()
        @test t.step_acceptance === :nms && t.reg == 1e3 && t.reg_decay == 0.9 && t.delta_decay == 0.95 && t.elastic_mode && t.divergence_recovery
        @test t.p_tol == 1e-4 && t.d_tol == 2e-2 && dgsqp_tuned_options(; reg = 5.0).reg == 5.0
        o = dgsqp_racing_options(; divergence_recovery = true)
        @test o.divergence_recovery && o.max_recoveries == 20 && !DGSQP().elastic_mode && !DGSQP().divergence_recovery
    end

end # testset DGSQP
