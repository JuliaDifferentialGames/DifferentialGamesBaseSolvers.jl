# ============================================================================
# dtmpc_tube_demo.jl
#
# Generic (non-jamming-game-specific) demonstration of DTMPC's Theory-Spine-v0.2
# faithful tube mechanism, producing publication figures on the SAME minimal
# 2-agent 1-D double-integrator PD-GNEP fixture used by the DTRHN theory
# validation tests (`test/theory_validation.jl`). Companion to (and generated
# the same way as) the on-orbit-jamming-game figures, but package-generic.
#
# Produces:
#   fig_tube_invariance.{pdf,png}  — closed-loop Monte-Carlo tube corridor
#                                    (analogous to Lopez–Slotine–How Fig. 3)
#   fig_tightening.{pdf,png}      — decision-dependent collision tightening,
#                                    WITH vs WITHOUT (§3, Eq. 5)
#   fig_theory_validation.{pdf,png} — 3-panel Props 2/3/4 Monte-Carlo summary
#                                      (constraint satisfaction, recursive
#                                      feasibility, cost-spread scaling)
#
# Requires FALCON (external, unregistered — this package has no FALCON
# dependency). Run via an environment with both packages, e.g. the shared
# ~/.julia/environments/v1.12 dev environment:
#
#   julia --project=@v1.12 -e '
#       using DifferentialGamesBase, DifferentialGamesBaseSolvers, FALCON, Plots
#       include("src/solvers/DTRHN/examples/dtmpc_tube_demo.jl")'
# ============================================================================

using LinearAlgebra
using Random
using Statistics
using Printf
using DifferentialGamesBase
using DifferentialGamesBaseSolvers
using DifferentialGamesBaseSolvers:
    create_independent_uncertainty_models, build_augmented_game,
    dtmpc_free_initial_mask, lsh_boundary_layer_control
using FALCON
using Plots

const FIG_DIR = joinpath(@__DIR__, "figures")
mkpath(FIG_DIR)
_save(fig, name) = (savefig(fig, joinpath(FIG_DIR, name * ".pdf")); savefig(fig, joinpath(FIG_DIR, name * ".png")))

gr()
default(fontfamily="Helvetica", framestyle=:box, grid=true, gridalpha=0.18, gridlinewidth=0.6,
        linewidth=2.0, legend=:best, dpi=300, guidefontsize=10, tickfontsize=8,
        legendfontsize=7, titlefontsize=11, foreground_color_legend=nothing,
        background_color_legend=RGBA(1, 1, 1, 0.7))
# Colorblind-safe (Okabe–Ito) palette, used consistently across every figure.
const C_A     = RGB(0/255, 114/255, 178/255)    # blue    — primary / "with DTMPC"
const C_B     = RGB(230/255, 159/255, 0/255)    # orange  — secondary
const C_TUBE  = RGB(0/255, 158/255, 115/255)    # green   — tube bound / guarantee
const C_TRIAL = RGB(0.55, 0.55, 0.55)           # gray    — Monte-Carlo ensembles
const C_ON    = RGB(0/255, 114/255, 178/255)    # blue    — protected / tightening ON
const C_OFF   = RGB(213/255, 94/255, 0/255)     # vermillion — unprotected / OFF

# ============================================================================
# Fixture: discrete 1-D double integrator, 2 agents, collision-avoidance coupling
# (same shape as DTRHN's test/theory_validation.jl, so this script exercises
# exactly what the test suite already validates — now rendered as figures).
# ============================================================================

const DT = 0.05
_di_discrete(xi, ui, p, t) = [xi[1] + DT * xi[2], xi[2] + DT * ui[1]]

function di_game(; n_players, N_win = 12, x0s = nothing,
                 q_pos = 1.0, q_vel = 0.1, r = 0.1, qf_pos = 20.0, qf_vel = 2.0)
    T = Float64
    stage(x, u, p, t) = q_pos * x[1]^2 + q_vel * x[2]^2 + r * u[1]^2
    term(x, p)        = qf_pos * x[1]^2 + qf_vel * x[2]^2
    players = PlayerSpec{T}[]
    for i in 1:n_players
        obj = PlayerObjective(i, NonlinearStageCost(stage), NonlinearTerminalCost(term))
        x0_i = x0s === nothing ? T[Float64(i), 0.0] : x0s[i]
        push!(players, PlayerSpec(i, 2, 1, x0_i, _di_discrete, obj))
    end
    return PDGNEProblem(players, T(N_win * DT), DT)
end

make_falcon(; kwargs...) = FALCONSolver(; r0=0.5, r_min=1e-4, r_max=5.0,
    ρ_0=0.0, ρ_1=0.25, ρ_2=0.70, w0=1.0, β=2.0, ω_max=1e4, η=1.0,
    ϵ_opt=1e-3, ϵ_feas=1e-4, max_iterations=200, kwargs...)

# ============================================================================
# Figure 1 — closed-loop tube invariance (Monte-Carlo, LSH Fig. 3 analogue)
# ============================================================================

function fig_tube_invariance(; N_win = 12, N_sim = 40, D_bound = 0.02, M = 30, α_fixed = 5.0)
    println("Solving closed-loop DTMPC (single agent, tube invariance demo)…")
    game = di_game(; n_players = 1, N_win = N_win)
    p = DTMPCParameters(Float64, 1; s₀=[0.3], α_fixed=[α_fixed], D=[D_bound], η=[0.02],
                       pos_dim=1, ancillary=:lsh_boundary_layer, λ=[2.0])
    um = create_independent_uncertainty_models(x -> D_bound, 1)
    x0 = [1.0, 0.0]
    dtmpc = DynamicTubeMPC(make_falcon(); warm_start=true, discrete_dynamics=true)
    sol = solve(DynamicTubeMPCProblem(game, x0, N_sim, p, um), dtmpc; verbose=false)
    @printf("  converged=%s\n", sol.converged)

    X̄ = sol.solver_info[:nominal_X]
    ū = sol.U[1]
    s = sol.s_trajectories[1]
    α = p.α_fixed[1]; λ = p.λ[1]

    rng = MersenneTwister(3)
    function trial()
        x = copy(x0)
        dev = Float64[norm(x .- X̄[:, 1])]
        for k in 1:N_sim
            u = lsh_boundary_layer_control(ū[:, k], x, X̄[:, k], s[k], α, λ, 1)
            xn = _di_discrete(x, u, nothing, (k - 1) * DT)
            d = [0.0, (2rand(rng) - 1) * D_bound * DT]
            x = xn .+ d
            push!(dev, norm(x .- X̄[:, k + 1]))
        end
        return dev
    end
    trials = [trial() for _ in 1:M]
    contained = all(all(tr[k] <= s[k] + 1e-9 for k in eachindex(tr)) for tr in trials)
    @printf("  contained=%s over %d trials\n", contained, M)

    tmin = (0:N_sim) .* DT
    p1 = plot(tmin, s; color=C_TUBE, lw=3, label="tube width s(t)",
              xlabel="time [s]", ylabel="|x - x_nominal|", title="Closed-loop tube invariance (Prop. 2)")
    for (i, tr) in enumerate(trials)
        plot!(p1, tmin, tr; color=C_TRIAL, lw=0.6, alpha=0.4, label = i == 1 ? "disturbed trials" : "")
    end
    _save(p1, "fig_tube_invariance")
    println("Wrote fig_tube_invariance")
    return contained
end

# ============================================================================
# Figure 2 — decision-dependent collision tightening (§3, Eq. 5), on vs off
# ============================================================================

function fig_tightening(; N_win = 20, D_bound = 0.05, d_min = 0.4)
    println("Solving single-window augmented games (tightening on/off)…")
    game = di_game(; n_players = 2, N_win = N_win, x0s = [[-1.0, 0.0], [1.0, 0.0]], qf_pos = 80.0)
    um = create_independent_uncertainty_models(x -> D_bound, 2)
    x0 = [-1.0, 0.0, 1.0, 0.0]

    function run(tighten)
        p = DTMPCParameters(Float64, 2; s₀=[0.2, 0.2], α_fixed=[1.0, 1.0], D=[D_bound, D_bound],
                            coupled_pairs=[(1, 2)], collision_min_distance=d_min, pos_dim=1,
                            constraint_tightening=tighten)
        aug = build_augmented_game(game, p, um, x0, p.s₀)
        sol = solve(aug, make_falcon())
        r1 = get_trajectory(sol, 1).states[1, :]
        r2 = get_trajectory(sol, 2).states[1, :]
        return abs.(r1 .- r2)
    end
    sep_on = run(true); sep_off = run(false)
    tmin = (0:N_win) .* DT
    @printf("  min separation: tightened=%.3f  untightened=%.3f  (d_min=%.2f)\n",
            minimum(sep_on), minimum(sep_off), d_min)

    p2 = plot(tmin, sep_on; color=C_ON, lw=2.4, label="tightening ON",
              xlabel="time [s]", ylabel="separation |r1-r2|",
              title="Decision-dependent tightening (§3, Eq. 5)")
    plot!(p2, tmin, sep_off; color=C_OFF, lw=2.4, ls=:dash, label="tightening OFF")
    hline!(p2, [d_min]; color=:black, ls=:dot, lw=1.4, label="d_min")
    _save(p2, "fig_tightening")
    println("Wrote fig_tightening")
end

# ============================================================================
# Figure 3 — 3-panel Theory Spine validation (Props 2/3/4), Monte-Carlo
# ============================================================================

function fig_theory_validation(; N_win = 12, N_sim = 20, M = 25)
    println("Running 3-way Monte-Carlo (Props 2/3/4)…")
    game = di_game(; n_players = 1, N_win = N_win)
    x0 = [1.0, 0.0]

    function run_trial(D_bound, seed)
        p = DTMPCParameters(Float64, 1; s₀=[0.2], α_fixed=[5.0], D=[D_bound], η=[0.01],
                            pos_dim=1, tighten_actuator=true, actuator_max=[3.0],
                            ancillary=:lsh_boundary_layer, λ=[2.0])
        um = create_independent_uncertainty_models(x -> D_bound, 1)
        rng = MersenneTwister(seed)
        dist = t -> [0.0, (2rand(rng) - 1) * D_bound * DT]
        dtmpc = DynamicTubeMPC(make_falcon(); warm_start=true, discrete_dynamics=true)
        sol = solve(DynamicTubeMPCProblem(game, x0, N_sim, p, um; disturbance=dist), dtmpc; verbose=false)
        actuator_ok = all(abs(sol.U[1][1, k]) <= 3.0 + 1e-6 for k in 1:N_sim)
        window_frac = sol.solver_info[:n_inner_converged] / N_sim
        return (; actuator_ok, window_frac, cost = sol.costs[1])
    end

    D_REF = 0.03
    results = [run_trial(D_REF, s) for s in 1:M]
    frac_constraint = mean(r.actuator_ok for r in results)
    frac_windows = mean(r.window_frac for r in results)
    @printf("  (A) Prop 2 constraint satisfaction: %.1f%%\n", 100 * frac_constraint)
    @printf("  (B) Prop 3 recursive feasibility:   %.1f%% windows converged\n", 100 * frac_windows)

    D_levels = [0.01, 0.03, 0.06]
    costs = Dict(D => [run_trial(D, 100 + s).cost for s in 1:M] for D in D_levels)
    spreads = [maximum(costs[D]) - minimum(costs[D]) for D in D_levels]
    @printf("  (C) Prop 4 cost spreads: %s\n", round.(spreads, digits=4))

    p3a = bar(["executed |a| <= a_max"], [100 * frac_constraint]; color=C_A, legend=false,
        ylims=(0, 105), ylabel="trials satisfied [%]", title="(A) Prop. 2")
    p3b = histogram([100 * r.window_frac for r in results]; color=C_B, bins=8, legend=false,
        xlabel="windows converged [%]", ylabel="# trials", title="(B) Prop. 3")
    p3c = plot(D_levels, spreads; color=C_A, marker=:circle, ms=6, lw=2.2, legend=false,
        xlabel="disturbance bound D", ylabel="cost spread", title="(C) Prop. 4: ε=O(diam D)")
    fig3 = plot(p3a, p3b, p3c; layout=(1, 3), size=(1200, 380), left_margin=5Plots.mm, bottom_margin=6Plots.mm)
    _save(fig3, "fig_theory_validation")
    println("Wrote fig_theory_validation")
end

# ============================================================================
# Figure 4 — combined 4-panel summary, designed to convey DTMPC's PRACTICAL
# EFFECT (not only validate the theory in the abstract):
#   (a) tube-bounded (DTMPC) vs. open-loop-nominal (no ancillary) trajectories
#       under the SAME disturbance realizations — the direct before/after
#       comparison of what the tube+ancillary buys (Prop. 2)
#   (b) decision-dependent collision tightening on/off (§3, Eq. 5) — the
#       tightened trajectory is the ONLY one that keeps the hard safety margin
#   (c) tube margin scaling with the disturbance bound D (Prop. 2 mechanism)
#   (d) recursive-feasibility containment margin (s−residual ≥ 0) at every
#       receding-horizon window transition, free-initial convention (§4/Prop. 3)
#
# Journal-figure conventions: colorblind-safe (Okabe–Ito) palette shared across
# all four panels; ASCII-only axis/legend text (no combining diacritics or
# subscript code points — GR's default font silently drops those glyphs);
# consistent panel titles/typography; vector PDF + high-resolution PNG.
# ============================================================================

function fig_four_panel(; N_win = 12, N_sim = 150, D_bound = 0.05, M = 30, α_fixed = 5.0)
    println("Building combined 4-panel DTMPC summary figure…")

    # This figure uses a much larger canvas (2000x1450) than the other three, so
    # its font sizes are scoped LOCALLY (not the file's global default) and set
    # noticeably larger: Plots.jl/GR's PDF page size equals the `size=(w,h)`
    # canvas in points (72/inch), so a wide working canvas gets shrunk by
    # roughly 0.25-0.4x when included at a typical single/double-column
    # linewidth in a journal template — fonts are sized up accordingly so the
    # ON-PAGE result lands comfortably in the 8-10 pt range, not 3-4 pt.
    default(fontfamily="Helvetica", framestyle=:box, grid=true, gridalpha=0.18, gridlinewidth=0.8,
            linewidth=3.0, legend=:best, dpi=300, guidefontsize=30, tickfontsize=24,
            legendfontsize=22, titlefontsize=34, foreground_color_legend=nothing,
            background_color_legend=RGBA(1, 1, 1, 0.85))

    # ── (a) WITH DTMPC ancillary vs. WITHOUT (open-loop nominal), same disturbances ──
    game1 = di_game(; n_players = 1, N_win = N_win)
    p1p = DTMPCParameters(Float64, 1; s₀=[0.3], α_fixed=[α_fixed], D=[D_bound], η=[0.02],
                         pos_dim=1, ancillary=:lsh_boundary_layer, λ=[2.0])
    um1 = create_independent_uncertainty_models(x -> D_bound, 1)
    x0_1 = [1.0, 0.0]
    dtmpc1 = DynamicTubeMPC(make_falcon(); warm_start=true, discrete_dynamics=true)
    sol_a = solve(DynamicTubeMPCProblem(game1, x0_1, N_sim, p1p, um1), dtmpc1; verbose=false)
    @printf("  (a) converged=%s\n", sol_a.converged)
    X_nom = sol_a.solver_info[:nominal_X]; u_nom = sol_a.U[1]; s_a = sol_a.s_trajectories[1]
    λa = p1p.λ[1]

    rng = MersenneTwister(3)
    function trial_a(with_ancillary::Bool)
        x = copy(x0_1); dev = Float64[norm(x .- X_nom[:, 1])]
        for k in 1:N_sim
            u = with_ancillary ?
                lsh_boundary_layer_control(u_nom[:, k], x, X_nom[:, k], s_a[k], α_fixed, λa, 1) :
                u_nom[:, k]
            xn = _di_discrete(x, u, nothing, (k - 1) * DT)
            d = [0.0, (2rand(rng) - 1) * D_bound * DT]
            x = xn .+ d
            push!(dev, norm(x .- X_nom[:, k + 1]))
        end
        return dev
    end
    # Same RNG stream reused for both arms (paired comparison, common random numbers)
    Random.seed!(rng, 3); trials_on  = [trial_a(true)  for _ in 1:M]
    Random.seed!(rng, 3); trials_off = [trial_a(false) for _ in 1:M]
    contained_on = all(all(tr[k] <= s_a[k] + 1e-9 for k in eachindex(tr)) for tr in trials_on)
    frac_off_violate = mean(any(tr[k] > s_a[k] for k in eachindex(tr)) for tr in trials_off)
    @printf("  (a) with-ancillary contained=%s | without-ancillary violates tube in %.0f%% of trials\n",
            contained_on, 100 * frac_off_violate)
    tmin_a = (0:N_sim) .* DT
    pa = plot(tmin_a, s_a; color=C_TUBE, lw=5.5, ls=:dash, label="tube bound s(t)",
              xlabel="Time [s]", ylabel="Deviation from nominal",
              title="(a)  Effect of the DTMPC ancillary", legend=:topleft)
    for (i, tr) in enumerate(trials_off)
        plot!(pa, tmin_a, tr; color=C_OFF, lw=1.2, alpha=0.35, label = i == 1 ? "without ancillary" : "")
    end
    for (i, tr) in enumerate(trials_on)
        plot!(pa, tmin_a, tr; color=C_ON, lw=1.2, alpha=0.55, label = i == 1 ? "with ancillary (DTMPC)" : "")
    end

    # ── (b) decision-dependent tightening on/off (§3, Eq. 5) ────────────────────
    N_win_b = 20; D_b = 0.05; d_min = 0.4
    game_b = di_game(; n_players = 2, N_win = N_win_b, x0s = [[-1.0, 0.0], [1.0, 0.0]], qf_pos = 80.0)
    um_b = create_independent_uncertainty_models(x -> D_b, 2)
    x0_b = [-1.0, 0.0, 1.0, 0.0]
    function run_b(tighten)
        p = DTMPCParameters(Float64, 2; s₀=[0.2, 0.2], α_fixed=[1.0, 1.0], D=[D_b, D_b],
                            coupled_pairs=[(1, 2)], collision_min_distance=d_min, pos_dim=1,
                            constraint_tightening=tighten)
        aug = build_augmented_game(game_b, p, um_b, x0_b, p.s₀)
        sol = solve(aug, make_falcon())
        r1 = get_trajectory(sol, 1).states[1, :]; r2 = get_trajectory(sol, 2).states[1, :]
        return abs.(r1 .- r2)
    end
    sep_on = run_b(true); sep_off = run_b(false)
    tmin_b = (0:N_win_b) .* DT
    @printf("  (b) min separation: tightened=%.3f  untightened=%.3f  (d_min=%.2f)\n",
            minimum(sep_on), minimum(sep_off), d_min)
    pb = plot(tmin_b, sep_on; color=C_ON, lw=4.3, label="tightening ON (DTMPC)",
              xlabel="Time [s]", ylabel="Separation |r1 - r2|",
              title="(b)  Decision-dependent tightening", legend=:topright)
    plot!(pb, tmin_b, sep_off; color=C_OFF, lw=4.3, ls=:dash, label="tightening OFF")
    hline!(pb, [d_min]; color=:black, ls=:dot, lw=2.7, label="safety limit d_min")

    # ── (c) margin scaling with disturbance bound D (Prop. 2 mechanism) ─────────
    game_c = di_game(; n_players = 1, N_win = N_win)
    x0_c = [1.0, 0.0]
    function final_s(Dval)
        p = DTMPCParameters(Float64, 1; s₀=[0.001], α_fixed=[1.0], D=[Dval], η=[0.0], pos_dim=1)
        um = create_independent_uncertainty_models(x -> 0.0, 1)
        aug = build_augmented_game(game_c, p, um, x0_c, p.s₀)
        sol = solve(aug, make_falcon())
        return get_trajectory(sol, 1).states[3, end]
    end
    D_levels = collect(range(0.01, 0.10, length = 7))
    s_finals = [final_s(D) for D in D_levels]
    @printf("  (c) steady-state tube widths: %s\n", round.(s_finals, digits=4))
    pc = plot(D_levels, s_finals; color=C_A, marker=:circle, markerstrokewidth=0, ms=11, lw=4.0,
        legend=false, xlabel="Disturbance bound D", ylabel="Steady-state tube width s",
        title="(c)  Tube margin vs. disturbance level")

    # ── (d) recursive-feasibility containment margin across RH windows (§4) ────
    game_d = di_game(; n_players = 1, N_win = N_win)
    p_d = DTMPCParameters(Float64, 1; s₀=[0.3], α_fixed=[α_fixed], D=[D_bound], η=[0.02],
                         pos_dim=1, ancillary=:lsh_boundary_layer, λ=[2.0])
    um_d = create_independent_uncertainty_models(x -> D_bound, 1)
    mask_d = dtmpc_free_initial_mask(p_d, [state_dim(game_d, 1)])
    rng_d = MersenneTwister(11)
    dist_d = t -> [0.0, (2rand(rng_d) - 1) * D_bound * DT]
    dtmpc_d = DynamicTubeMPC(make_falcon(; free_initial = mask_d); warm_start=true, discrete_dynamics=true)
    sol_d = solve(DynamicTubeMPCProblem(game_d, x0_1, N_sim, p_d, um_d; disturbance=dist_d), dtmpc_d; verbose=false)
    X_d = sol_d.X; Xn_d = sol_d.solver_info[:nominal_X]; s_d = sol_d.s_trajectories[1]
    residual_d = [norm(X_d[:, k] .- Xn_d[:, k]) for k in 1:N_sim+1]
    margin_d = s_d .- residual_d
    margin_ok = all(margin_d .>= -1e-6)
    @printf("  (d) recursive-feasibility margin stays non-negative at every window: %s\n", margin_ok)
    tmin_d = (0:N_sim) .* DT
    pd = plot(tmin_d, margin_d; color=C_TUBE, lw=0, fillrange=0, fillalpha=0.35, fillcolor=C_TUBE,
              label="", xlabel="Time [s]", ylabel="Containment margin  (s - residual)",
              title="(d)  Recursive feasibility (free-initial, Sec. 4)", legend=:topright)
    plot!(pd, tmin_d, margin_d; color=C_TUBE, lw=3.7, label="margin at each window")
    hline!(pd, [0.0]; color=:black, ls=:dash, lw=2.3, label="feasibility boundary")

    fig4 = plot(pa, pb, pc, pd; layout=(2, 2), size=(2000, 1450), dpi=300,
                background_color=:white, background_color_subplot=:white,
                left_margin=10Plots.mm, right_margin=6Plots.mm,
                bottom_margin=12Plots.mm, top_margin=8Plots.mm)
    _save(fig4, "fig_four_panel")
    println("Wrote fig_four_panel")

    # Restore the file's default (smaller-canvas) plot theme, in case this
    # script is used interactively rather than run once end-to-end.
    default(fontfamily="Helvetica", framestyle=:box, grid=true, gridalpha=0.18, gridlinewidth=0.6,
            linewidth=2.0, legend=:best, dpi=300, guidefontsize=10, tickfontsize=8,
            legendfontsize=7, titlefontsize=11, foreground_color_legend=nothing,
            background_color_legend=RGBA(1, 1, 1, 0.7))
end

# ============================================================================
fig_tube_invariance()
fig_tightening()
fig_theory_validation()
fig_four_panel()
println("All DTRHN figures written to $(FIG_DIR)")
