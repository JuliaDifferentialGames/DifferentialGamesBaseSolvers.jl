# ============================================================================
# dtmpc_tests.jl — Test suite for the Dynamic Tube MPC (DTMPC) solver
#
# Validates the augmented-PD-GNEP implementation. Tube STRUCTURE follows the
# Theory Spine v0.2; tube BANDWIDTH follows the later, authoritative unified
# manuscript (Problem 2 / Eq. 9-10): each agent's tube channel carries a
# width state `s` AND a co-optimized bandwidth state `α` with its own
# artificial control `v` (`α̇=v`), advanced by the EXACT tight recursion
# `s_{k+1}=ρ(α)·s_k+δ̂(x̄,ū,α)`. The width still carries no objective
# curvature (Theorem 1) — only `s` and `α` are unpenalized; the tightening
# constraints that read `α·s` are consequently NONCONVEX (bilinear), unlike
# the fixed-α predecessor of this file.
#
# Inner solver: FALCON (the theory-faithful path; discrete SCP/PD-GNEP). Run
# via the shared v1.12 environment, which has both DifferentialGamesBaseSolvers
# and FALCON available (this package itself has no FALCON dependency):
#   julia --project=~/.julia/environments/v1.12 -e '
#       using DifferentialGamesBase, DifferentialGamesBaseSolvers, FALCON
#       include("src/solvers/DTRHN/test/dtmpc_tests.jl")'
# ============================================================================

using Test
using LinearAlgebra
using Random
using DifferentialGamesBase
using DifferentialGamesBaseSolvers
using FALCON
using DifferentialGamesBaseSolvers:
    build_augmented_player_dynamics_discrete, build_augmented_cost,
    build_tightened_collision_constraint, build_tightened_actuator_constraint,
    build_initial_containment_constraint, build_augmented_game,
    agent_physical_block, block_offsets, augmented_state_dims,
    augmented_control_dims, extract_agent_augmented, tube_channels,
    n_tube_channels_agent, tighten_coupled_distance, check_coupled_constraint_satisfaction,
    create_independent_uncertainty_models, create_coupled_uncertainty_models,
    dtmpc_free_initial_mask, ρ_of, _physical_zfull, dtmpc_bandwidth_trajectory

Random.seed!(42)

# ============================================================================
# Fixtures — DISCRETE PD-GNEP games with NONLINEAR (FALCON) cost form.
# ============================================================================

"""
    make_pdgnep_di(; n_players, N_win, dt) -> GameProblem

`n_players` agents, each a DISCRETE 2D double integrator (nᵢ=2, mᵢ=1):
`x_{k+1} = [x[1]+dt·x[2]; x[2]+dt·u[1]]`, regulated to the origin (nonlinear
cost form, required by the augmented-cost builder's FALCON path).
"""
function make_pdgnep_di(; n_players::Int = 2, N_win::Int = 8, dt::Float64 = 0.1,
                        q_pos = 1.0, q_vel = 0.1, r = 0.5, qf_pos = 5.0, qf_vel = 1.0,
                        shared_constraints::AbstractVector = AbstractSharedConstraint[])
    T = Float64
    di_discrete(xi, ui, p, t) = [xi[1] + dt * xi[2], xi[2] + dt * ui[1]]
    stage(x, u, p, t) = q_pos * x[1]^2 + q_vel * x[2]^2 + r * u[1]^2
    term(x, p) = qf_pos * x[1]^2 + qf_vel * x[2]^2
    players = PlayerSpec{T}[]
    for i in 1:n_players
        obj = PlayerObjective(i, NonlinearStageCost(stage), NonlinearTerminalCost(term))
        x0_i = T[Float64(i), 0.0]              # spread agents along the line
        push!(players, PlayerSpec(i, 2, 1, x0_i, di_discrete, obj))
    end
    return PDGNEProblem(players, shared_constraints, T(N_win * dt), dt)
end

"""
make a default single-channel tube params object for n agents. `pos_dim=1`:
the fixture is a 1-D double integrator (state = [position; velocity], nᵢ=2),
so the translational channel's state slice `2·pos_dim` must equal 2.
"""
function make_params(n; kwargs...)
    T = Float64
    DTMPCParameters(T, n;
        s₀ = fill(T(0.1), n), α_fixed = fill(T(1.0), n),
        D = fill(T(0.05), n), η = fill(T(0.01), n), pos_dim = 1, kwargs...)
end

"default small FALCON solver for the tiny double-integrator regulation problems"
make_falcon(; kwargs...) = FALCONSolver(; r0=0.5, r_min=1e-4, r_max=5.0,
    ρ_0=0.0, ρ_1=0.25, ρ_2=0.70, w0=1.0, β=2.0, ω_max=1e4, η=1.0,
    ϵ_opt=1e-3, ϵ_feas=1e-4, max_iterations=150, kwargs...)

"build a full DTMPC problem from a discrete PD-GNEP DI game"
function make_dtmpc_problem(; n_players = 2, N_win = 8, N_sim = 6, kwargs...)
    T = Float64
    game = make_pdgnep_di(; n_players = n_players, N_win = N_win)
    um = create_independent_uncertainty_models(x -> T(0.05), n_players)
    params = make_params(n_players; kwargs...)
    x0 = vcat([T[Float64(i), 0.0] for i in 1:n_players]...)
    prob = DynamicTubeMPCProblem(game, x0, N_sim, params, um)
    return prob, params, game
end

# ============================================================================
@testset "DTMPC" begin

    # ── A1: augmented layout + accessors ────────────────────────────────────
    @testset "A1 augmented layout/accessors" begin
        # x̃ᵢ = [xᵢ(nᵢ); s₁…s_nch; α₁…α_nch] — TWO scalars per channel (width + bandwidth).
        @test augmented_state_dims([2, 3], [1, 1]) == [4, 5]     # +2 (single channel: s,α)
        @test augmented_state_dims([2, 3], 1) == [4, 5]
        # ũᵢ = [uᵢ(mᵢ); v₁…v_nch] — one artificial bandwidth-rate control per channel.
        @test augmented_control_dims([1, 2], [1, 1]) == [2, 3]
        @test block_offsets([4, 5]) == [0, 4]

        x_joint = [1.0, 2.0, 3.0, 4.0]      # two agents, n=2 each
        @test agent_physical_block(x_joint, [2, 2], 1) == [1.0, 2.0]
        @test agent_physical_block(x_joint, [2, 2], 2) == [3.0, 4.0]
        @test agent_physical_block([1.0,2.0,3.0, 4.0,5.0,6.0], [3, 3], 2) == [4.0, 5.0, 6.0]

        # Per-agent augmented block: [pos, vel, s, α] (α appended AFTER the width block).
        aug_offsets = [0, 4]
        joint_aug = [10.0, 11.0, 0.3, 1.5,  20.0, 21.0, 0.5, 2.5]
        xi, s, α = extract_agent_augmented(joint_aug, aug_offsets, [2, 2], [1, 1], 1)
        @test xi == [10.0, 11.0] && s == [0.3] && α == [1.5]
        xi2, s2, α2 = extract_agent_augmented(joint_aug, aug_offsets, [2, 2], [1, 1], 2)
        @test xi2 == [20.0, 21.0] && s2 == [0.5] && α2 == [2.5]
    end

    # ── A2: params ──────────────────────────────────────────────────────────
    @testset "A2 DTMPCParameters" begin
        T = Float64
        p = DTMPCParameters(T, 3)
        @test length(p.s₀) == 3 && length(p.α_fixed) == 3
        @test all(p.s₀ .== 0.1) && all(p.α_fixed .== 1.0)
        @test p.constraint_tightening
        @test n_tube_channels_agent(p, 1) == 1
        @test length(tube_channels(p, 1)) == 1
        @test tube_channels(p, 1)[1].kind === :translational

        p2 = DTMPCParameters(T, 2; s₀ = [0.2, 0.3],
                             coupled_pairs = [(1, 2)], collision_min_distance = 1.5, pos_dim = 1)
        @test p2.s₀ == [0.2, 0.3]
        @test p2.coupled_pairs == [(1, 2)] && p2.pos_dim == 1

        p3 = DTMPCParameters(T, 2; rotational_tube = [true, false], detection_tube = [true, false])
        @test n_tube_channels_agent(p3, 1) == 3
        @test n_tube_channels_agent(p3, 2) == 1
        @test [c.kind for c in tube_channels(p3, 1)] == [:translational, :rotational, :detection]

        @test_throws AssertionError DTMPCParameters(T, 2; s₀ = [0.1])
        @test_throws AssertionError DTMPCParameters(T, 2; α_fixed = [0.0, 1.0])
        @test_throws AssertionError DTMPCParameters(T, 2; coupled_pairs = [(1, 1)])

        # Co-optimized bandwidth box (manuscript Eq. 10, α̇=v).
        @test p.α_lo ≈ 1e-2 && p.α_max ≈ 5.0 && p.v_max ≈ 1.0   # defaults
        @test_throws AssertionError DTMPCParameters(T, 2; α_lo = 0.0)          # unstable tube ODE
        @test_throws AssertionError DTMPCParameters(T, 2; α_lo = -0.1)
        @test_throws AssertionError DTMPCParameters(T, 2; α_max = 0.005)       # α_max < α_lo
        @test_throws AssertionError DTMPCParameters(T, 2; v_max = 0.0)
        @test_throws AssertionError DTMPCParameters(T, 2; α_fixed = [10.0, 1.0])  # α₀ outside [α_lo,α_max]
    end

    # ── A3: augmented dynamics (single channel, CO-OPTIMIZED α, EXACT recursion) ───
    @testset "A3 augmented dynamics" begin
        T = Float64
        dt = 0.1
        di_discrete(xi, ui, p, t) = [xi[1] + dt * xi[2], xi[2] + dt * ui[1]]
        p = make_params(1)
        chans = tube_channels(p, 1)
        α_lo = p.α_lo
        f_aug = build_augmented_player_dynamics_discrete(
            di_discrete, 2, 1, chans, T[], dt, [2], [1], 1, [1]; α_lo = α_lo)

        # x̃ = [pos, vel, s, α]; ũ = [u, v] (bandwidth-rate control).
        x̃ = T[1.0, 2.0, 0.3, 1.2]; ũ = T[0.5, 0.4]
        x_next = f_aug(x̃, ũ, (player=1, timestep=0), 0.0)
        @test length(x_next) == 4
        @test x_next[1:2] ≈ di_discrete(x̃[1:2], ũ[1:1], nothing, 0.0)

        α = x̃[4]
        ρ = ρ_of(α, dt)
        β = chans[1].D + chans[1].η                     # continuous-time rate bound
        δ̂ = (β / α) * (1 - ρ)                            # EXACT discrete-step integral at THIS α
        # atol loosened from the exact-arithmetic 1e-10: the implementation applies the smooth
        # positive guard √(α²+(α_lo/10)²) before dividing by α (see `_α_guard`), which at
        # α=1.2 (well away from zero) differs from the raw α used here by ~1e-7 relative —
        # negligible physically, but enough to fail a bit-exact comparison.
        @test x_next[3] ≈ ρ * 0.3 + δ̂  atol=1e-6
        # α̇ = v (manuscript Eq. 10), exact under ZOH v.
        @test x_next[4] ≈ α + dt * ũ[2]  atol=1e-12

        # x_other_frozen has no effect on a channel with κ_speed=0 / Δ_funcs=nothing.
        f_aug2 = build_augmented_player_dynamics_discrete(
            di_discrete, 2, 1, chans, T[99.0, 99.0], dt, [2], [1], 1, [1]; α_lo = α_lo)
        @test f_aug2(x̃, ũ, (player=1, timestep=0), 0.0) ≈ x_next

        # Guard against α ≤ 0 trial iterates: division by α must stay finite even when the
        # decision variable transiently sits outside [α_lo, α_max] mid-SCP (the box is
        # enforced by StateBounds, not by this dynamics function).
        x̃_bad = T[1.0, 2.0, 0.3, 0.0]
        @test all(isfinite, f_aug(x̃_bad, ũ, (player=1, timestep=0), 0.0))
        x̃_neg = T[1.0, 2.0, 0.3, -0.5]
        @test all(isfinite, f_aug(x̃_neg, ũ, (player=1, timestep=0), 0.0))
    end

    # ── A4: augmented cost carries NO tube penalty (Theorem 1) ──────────────
    @testset "A4 augmented cost (no width curvature)" begin
        T = Float64
        stage(x,u,p,t) = dot(x,x) + 0.1*dot(u,u)
        term(x,p) = 5.0*dot(x,x)
        obj = PlayerObjective(1, NonlinearStageCost(stage), NonlinearTerminalCost(term))
        aug = build_augmented_cost(obj, 2, 1; sdims=[2], cdims=[1], N=1, nt_vec=[1])
        x̃ = T[1.0, 2.0, 0.3]; ũ = T[0.5]
        base = stage(T[1.0,2.0], T[0.5], nothing, 1)
        got  = evaluate_stage_cost(aug.stage_cost, x̃, ũ, nothing, 1)
        @test got ≈ base
        baseT = term(T[1.0,2.0], nothing)
        gotT  = evaluate_terminal_cost(aug.terminal_cost, x̃, nothing)
        @test gotT ≈ baseT
    end

    # ── A5: decision-dependent tightened collision constraint ───────────────
    @testset "A5 tightened collision constraint" begin
        T = Float64
        aug_offsets = [0, 3]
        c = build_tightened_collision_constraint(1, 2, aug_offsets, 2, 2, 1, T(0.5))
        @test get_players(c) == [1, 2]
        x = T[2.0, 0.0, 0.1,  1.0, 0.0, 0.2]
        val = evaluate_constraint(c, x, T[], nothing, 1)[1]
        @test val ≈ 0.5 - 1.0 + 0.1 + 0.2 atol = 1e-6
        x2 = copy(x); x2[3] = 0.5; x2[6] = 0.5
        @test evaluate_constraint(c, x2, T[], nothing, 1)[1] > val
    end

    # ── A6: tightened actuator constraint is NONCONVEX now that α is a state ──
    @testset "A6 tightened actuator constraint (bilinear α·s, nonconvex)" begin
        T = Float64
        # (i, n_i, nch, chan_idx, pos_dim, a_max); reads α from x[n_i+nch+chan_idx].
        c = build_tightened_actuator_constraint(1, 2, 1, 1, 1, T(1.0))
        @test !is_convex(c)   # bilinear α·s ⇒ no longer on FALCON's hard h_cvx path
        # x = [pos, vel, s, α] = [0.0, 0.0, 0.1, 2.0]
        val = evaluate_constraint(c, T[0.0,0.0,0.1,2.0], T[0.3], (player=1,timestep=0), 0.0)[1]
        @test val ≈ 0.3 + 2.0*0.1 - 1.0 atol=1e-6
        # Same s, different α ⇒ different reserve (this IS the bilinearity: the constraint's
        # sensitivity to s now depends on the current value of α, not a fixed constant).
        val2 = evaluate_constraint(c, T[0.0,0.0,0.1,4.0], T[0.3], (player=1,timestep=0), 0.0)[1]
        @test val2 ≈ 0.3 + 4.0*0.1 - 1.0 atol=1e-6
        @test val2 > val
    end

    # ── A7: initial-containment constraint (§4) ─────────────────────────────
    @testset "A7 initial containment" begin
        T = Float64
        c = build_initial_containment_constraint(1, 2, 1, 0, 2, T[1.0, 0.0])
        @test is_convex(c)
        v0 = evaluate_constraint(c, T[1.0,0.0,0.0], T[], (player=1,timestep=0), 0.0)[1]
        @test v0 ≤ 1e-3
        v1 = evaluate_constraint(c, T[1.5,0.0,0.1], T[], (player=1,timestep=0), 0.0)[1]
        @test v1 > 0
        v2 = evaluate_constraint(c, T[100.0,0.0,0.0], T[], (player=1,timestep=1), 0.0)[1]
        @test v2 < 0
    end

    # ── B1: augmented game assembly ─────────────────────────────────────────
    @testset "B1 build_augmented_game" begin
        T = Float64
        game = make_pdgnep_di(n_players = 2)
        params = make_params(2; coupled_pairs = [(1, 2)], collision_min_distance = 0.5, pos_dim = 1)
        um = create_independent_uncertainty_models(x -> T(0.05), 2)
        x0 = T[1.0, 0.0, 2.0, 0.0]
        aug = build_augmented_game(game, params, um, x0, fill(T(0.1), 2))
        @test is_pd_gnep(aug)
        @test has_separable_dynamics(aug)
        @test state_dim(aug) == 8                       # 2*(2+2) — width + bandwidth per channel
        @test [control_dim(aug, i) for i in 1:2] == [2, 2]  # physical(1) + bandwidth-rate v(1)
        @test length(aug.shared_constraints) == 1
    end

    # ── B1b: user shared constraints survive augmentation ──────────────────
    # Regression test: `build_augmented_game` used to silently DROP every user-supplied
    # shared constraint (`shared = Any[]`, repopulated only with the tightened collision
    # constraint). A coupled constraint on the ORIGINAL game — e.g. the jamming game's
    # detection cap — would then exist in an open-loop FALCON solve but vanish the moment
    # the same game was wrapped in DTMPC.
    @testset "B1b build_augmented_game carries shared constraints" begin
        T = Float64

        # A trivial custom shared constraint coupling both agents' positions, on a game
        # with constraint_tightening OFF so exactly one wrapped constraint survives.
        mycon = SharedInequality([1, 2]; dim = 1, func = (x, u, p, t) -> [x[1] - x[3] - 10.0])
        game = make_pdgnep_di(n_players = 2, shared_constraints = [mycon])
        params = make_params(2; constraint_tightening = false)
        um = create_independent_uncertainty_models(x -> T(0.05), 2)
        x0 = T[1.0, 0.0, 2.0, 0.0]
        aug = build_augmented_game(game, params, um, x0, fill(T(0.1), 2))
        @test length(aug.shared_constraints) == 1

        # Augmented joint state: agent block = [pos, vel, s, α] (nch=1), aug_offsets=[0,4].
        # Augmented joint control: agent block = [u, v] (nch=1), cu_offsets=[0,2].
        x_aug = T[1.0, 0.0, 0.3, 1.0,  2.0, 0.0, 0.5, 1.0]
        u_aug = T[0.1, 0.05,  0.2, 0.07]   # [u1,v1, u2,v2]
        p_arg = (player = 1, timestep = 0)
        val = evaluate_constraint(aug.shared_constraints[1], x_aug, u_aug, p_arg, 0.0)
        x_phys = T[x_aug[1], x_aug[2], x_aug[5], x_aug[6]]   # physical slice, same order
        u_phys = T[u_aug[1], u_aug[3]]                       # physical control, v stripped
        expected = evaluate_constraint(mycon, x_phys, u_phys, p_arg, 0.0)
        @test val ≈ expected

        # A user ProximityConstraint on the SAME pair the tightened collision constraint
        # already covers must be SKIPPED (not double-enforced).
        prox_dup = ProximityConstraint([1, 2], 0, 2, 1, 0.5)
        game2 = make_pdgnep_di(n_players = 2, shared_constraints = [prox_dup])
        params2 = make_params(2; coupled_pairs = [(1, 2)], collision_min_distance = 0.5, pos_dim = 1)
        aug2 = build_augmented_game(game2, params2, um, x0, fill(T(0.1), 2))
        @test length(aug2.shared_constraints) == 1   # tightened only, dup skipped

        # A user ProximityConstraint on a DIFFERENT pair than the tightened one must
        # survive alongside it.
        game3 = make_pdgnep_di(n_players = 3, shared_constraints = [ProximityConstraint([1, 3], 0, 4, 1, 0.5)])
        params3 = make_params(3; coupled_pairs = [(1, 2)], collision_min_distance = 0.5, pos_dim = 1)
        um3 = create_independent_uncertainty_models(x -> T(0.05), 3)
        x03 = T[1.0, 0.0, 2.0, 0.0, 3.0, 0.0]
        aug3 = build_augmented_game(game3, params3, um3, x03, fill(T(0.1), 3))
        @test length(aug3.shared_constraints) == 2   # user's (1,3) + tightened (1,2)
    end

    # ── B/integration: full solve with FALCON ───────────────────────────────
    @testset "Integration: 2-agent solve (FALCON)" begin
        T = Float64
        prob, params, game = make_dtmpc_problem(n_players = 2, N_win = 8, N_sim = 6)
        sol = solve(prob, DynamicTubeMPC(make_falcon()))

        @test sol isa DynamicTubeMPCSolution{T}
        @test sol.converged
        @test size(sol.X) == (state_dim(game), prob.n_sim_steps + 1)
        @test length(sol.U) == 2
        @test sol.X[:, 1] ≈ prob.x0

        for i in 1:2
            @test length(sol.s_trajectories[i]) == prob.n_sim_steps + 1
            @test sol.s_trajectories[i][1] ≈ params.s₀[i]
            @test all(sol.s_trajectories[i] .≥ 0)
            @test isfinite(sol.costs[i])
        end
        @test any(abs.(diff(sol.s_trajectories[1])) .> 1e-8)

        # α is CO-OPTIMIZED (manuscript Eq. 10): the solver has a real artificial control v
        # to exercise, so unlike the removed fixed-α design, α(t) is not required to be
        # constant. Assert it is genuinely a decision variable — bounds respected, and (the
        # actual content check, mirroring what the original co-optimized-α implementation's
        # own test suite asserted) not pinned at its initial value throughout.
        α_traj = sol.solver_info[:α_trajectories]
        @test length(α_traj) == 2
        for i in 1:2
            @test length(α_traj[i]) == prob.n_sim_steps + 1
            @test α_traj[i][1] ≈ params.α_fixed[i]              # seeded from α₀
            @test all(params.α_lo - 1e-9 .≤ α_traj[i] .≤ params.α_max + 1e-9)
            @test dtmpc_bandwidth_trajectory(sol, i) === α_traj[i]
        end
        @test any(abs.(diff(α_traj[1])) .> 1e-10) || any(abs.(diff(α_traj[2])) .> 1e-10)
    end

    # ── _physical_zfull round-trip: the edit most likely to fail SILENTLY ───
    # `_physical_zfull` feeds the LIVE cross-player z_full coupling both games read
    # opponent state/control out of. A stride mistake does not error — it hands the stage
    # costs garbage. Tested directly rather than left to downstream symptoms.
    @testset "_physical_zfull round-trip" begin
        T = Float64
        sdims = [2, 3]; cdims = [1, 2]; nch = [1, 2]; N = 2
        # Agent 1: n_aug=2+2*1=4, m_aug=1+1=2 per node; Agent 2: n_aug=3+2*2=7, m_aug=2+2=4.
        # Build a fully-populated augmented z_full for each agent: N+1 state blocks, N control
        # blocks, distinct values everywhere so any stride error shows up as a mismatch.
        function fake_zfull(n_aug, m_aug, N)
            z = T[]
            for k in 0:N-1
                append!(z, T.(1:n_aug) .+ 100k)         # state block k
                append!(z, T.(1:m_aug) .+ 100k .+ 50)   # control block k
            end
            append!(z, T.(1:n_aug) .+ 100N)             # final state block
            return z
        end
        n_aug1, m_aug1 = sdims[1] + 2nch[1], cdims[1] + nch[1]
        n_aug2, m_aug2 = sdims[2] + 2nch[2], cdims[2] + nch[2]
        aug_zfull = [fake_zfull(n_aug1, m_aug1, N), fake_zfull(n_aug2, m_aug2, N)]

        phys = _physical_zfull(aug_zfull, sdims, cdims, N, nch)
        @test length(phys) == 2
        for j in 1:2
            n, m = sdims[j], cdims[j]
            pstr = n + m
            @test length(phys[j]) == pstr * N + n
            for k in 0:N-1
                pb = k * pstr
                # Physical state block must equal the FIRST n entries of the augmented state
                # block at this node (widths/bandwidths stripped).
                @test phys[j][pb+1:pb+n] == T.(1:n) .+ 100k
                # Physical control block must equal the FIRST m entries of the augmented
                # control block (bandwidth-rate v stripped) — not the raw offset, which is
                # exactly the bug this test guards: before the fix the control slice was
                # never re-offset for the +nch control shift at all.
                @test phys[j][pb+n+1:pb+n+m] == T.(1:m) .+ 100k .+ 50
            end
            pb = N * pstr
            @test phys[j][pb+1:pb+n] == T.(1:n) .+ 100N
        end
    end

    @testset "Integration: 3-agent solve" begin
        prob, _, game = make_dtmpc_problem(n_players = 3, N_win = 6, N_sim = 5)
        sol = solve(prob, DynamicTubeMPC(make_falcon()))
        @test sol.converged
        @test dtmpc_n_players(sol) == 3
        @test size(sol.X, 2) == prob.n_sim_steps + 1
    end

    @testset "Integration: coupled collision constraint" begin
        T = Float64
        game = make_pdgnep_di(n_players = 2, N_win = 8)
        um = create_independent_uncertainty_models(x -> T(0.05), 2)
        params = make_params(2; coupled_pairs = [(1, 2)], collision_min_distance = T(0.5), pos_dim = 1)
        x0 = T[-1.0, 0.0, 1.0, 0.0]
        prob = DynamicTubeMPCProblem(game, x0, 6, params, um)
        sol = solve(prob, DynamicTubeMPC(make_falcon()))
        @test sol.converged
        @test dtmpc_coupled_pairs(sol) == [(1, 2)]
        @test sol.s_trajectories[1] ≈ sol.s_trajectories[2]
    end

    # ── warm-start parity (deterministic: zero disturbance) ─────────────────
    @testset "Warm-start parity" begin
        prob, _, _ = make_dtmpc_problem(n_players = 2, N_win = 6, N_sim = 5)
        sol_ws  = solve(prob, DynamicTubeMPC(make_falcon(); warm_start = true))
        sol_nws = solve(prob, DynamicTubeMPC(make_falcon(); warm_start = false))
        @test sol_ws.converged && sol_nws.converged
        @test sol_ws.X ≈ sol_nws.X atol = 1e-3
    end

    # ── accessors ────────────────────────────────────────────────────────────
    @testset "Accessors" begin
        prob, _, _ = make_dtmpc_problem(n_players = 2, N_win = 6, N_sim = 4)
        sol = solve(prob, DynamicTubeMPC(make_falcon()))
        @test get_state_trajectory(sol) === sol.X
        @test get_control_trajectory(sol, 1) === sol.U[1]
        @test get_total_cost(sol, 2) == sol.costs[2]
        @test dtmpc_tube_trajectory(sol, 1) === sol.s_trajectories[1]
        @test dtmpc_simulation_steps(sol) == prob.n_sim_steps
        @test_throws ErrorException get_control_trajectory(sol, 3)
        @test_throws ErrorException get_total_cost(sol, 0)
    end

    # ── ancillary boundary-layer control (α now CARRIED from the closed loop) ──
    @testset "Boundary-layer ancillary control" begin
        T = Float64
        ū = T[0.5, 0.0]; x̄ = T[1.0, 0.5]; pos_dim = 1
        @test lsh_boundary_layer_control(ū, x̄, x̄, T(0.2), T(2.0), T(0.1), pos_dim) ≈ ū
        x = T[1.3, 0.5]
        u = lsh_boundary_layer_control(ū, x, x̄, T(0.2), T(2.0), T(0.1), pos_dim)
        @test u != ū
    end

    # ── free-initial mask (§4) ───────────────────────────────────────────────
    @testset "dtmpc_free_initial_mask" begin
        T = Float64
        p = DTMPCParameters(T, 2; rotational_tube = [true, false])
        mask = dtmpc_free_initial_mask(p, [4, 3])
        @test length(mask) == 2
        @test length(mask[1]) == 4 + 2*2   # physical(4) + 2 channels × (s,α)
        @test length(mask[2]) == 3 + 2*1   # physical(3) + 1 channel × (s,α)
        @test all(all, mask)
    end

    # ── verification helpers ────────────────────────────────────────────────
    @testset "Coupled-constraint helpers" begin
        T = Float64
        @test tighten_coupled_distance(T(1.0), T(0.1), T(0.2)) ≈ 1.3
        @test check_coupled_constraint_satisfaction(T[0.0], T[2.0], T(0.1), T(0.2), T(1.0))
        @test !check_coupled_constraint_satisfaction(T[0.0], T[1.2], T(0.1), T(0.2), T(1.0))
    end

    # ── shared-state game is rejected (requires PD-GNEP) ────────────────────
    @testset "Shared-state game rejected" begin
        T = Float64
        A = T[1.0 0.1; 0.0 1.0]; B1 = T[0.0; 0.1;;]; B2 = T[0.0; 0.1;;]
        Q = diagm(T[1.0, 0.1]); Qf = diagm(T[5.0, 1.0]); R = reshape(T[0.5], 1, 1)
        shared = LQGameProblem(A, [B1, B2], [Q, Q], [R, R], [Qf, Qf], zeros(T, 2), T(0.8); dt = 0.1)
        um = create_independent_uncertainty_models(x -> T(0.05), 2)
        params = make_params(2)
        @test_throws AssertionError DynamicTubeMPCProblem(shared, zeros(T, 2), 5, params, um)
    end

    # ── failure-path guard (session correctness fix) ─────────────────────────
    # A non-converged inner solve must never drive the ancillary boundary-layer
    # controllers or advance the measured state (previously confirmed to crash
    # downstream with a NaN/Inf matrix inside lsh_rotational_boundary_layer_control
    # after two consecutive non-converged windows on a live 100-window disturbed
    # run). Forces a STRUCTURALLY infeasible tightened actuator bound
    # (alpha_fixed*s0 >> actuator_max, infeasible even at zero control -- the
    # exact bug pattern this session's earlier s0-scale fix addressed at the
    # scenario level) so every window's inner solve fails deterministically,
    # both on the first attempt AND the fresh retry, and confirms the run
    # terminates cleanly (truncated, correctly-sized solution + diagnostics)
    # rather than crashing or silently propagating a bad trajectory.
    @testset "Failure-path guard: clean early termination, no crash" begin
        T = Float64
        n = 2
        params = make_params(n; s₀ = fill(T(10.0), n), α_fixed = fill(T(1.0), n),
                              tighten_actuator = true, actuator_max = fill(T(0.01), n))
        game = make_pdgnep_di(; n_players = n, N_win = 8)
        um = create_independent_uncertainty_models(x -> T(0.05), n)
        x0 = vcat([T[Float64(i), 0.0] for i in 1:n]...)
        N_sim = 5
        prob = DynamicTubeMPCProblem(game, x0, N_sim, params, um)

        sol = solve(prob, DynamicTubeMPC(make_falcon()); verbose = false)

        @test sol.solver_info[:terminated_early] == true
        @test sol.solver_info[:termination_reason] isa String
        n_completed = sol.solver_info[:n_steps_simulated]
        @test n_completed < N_sim   # genuinely truncated, not the full horizon
        @test size(sol.X, 2) == n_completed + 1
        @test all(size(sol.U[i], 2) == n_completed for i in 1:n)
        @test all(isfinite, sol.X)
        @test all(all(isfinite, sol.U[i]) for i in 1:n)
        @test sol.converged == false
    end

    @testset "recovery_callback: nothing behaves identically to no callback" begin
        T = Float64
        n = 2
        params = make_params(n; s₀ = fill(T(10.0), n), α_fixed = fill(T(1.0), n),
                              tighten_actuator = true, actuator_max = fill(T(0.01), n))
        game = make_pdgnep_di(; n_players = n, N_win = 8)
        um = create_independent_uncertainty_models(x -> T(0.05), n)
        x0 = vcat([T[Float64(i), 0.0] for i in 1:n]...)
        N_sim = 5
        prob = DynamicTubeMPCProblem(game, x0, N_sim, params, um)

        always_nothing(; kwargs...) = nothing
        sol = solve(prob, DynamicTubeMPC(make_falcon(); recovery_callback = always_nothing); verbose = false)

        @test sol.solver_info[:terminated_early] == true
        @test sol.converged == false
    end

    @testset "recovery_callback: a genuinely feasible recovery lets the window complete" begin
        # SAME structurally-infeasible setup as the failure-path guard test above (guaranteed
        # double-failure at the primary actuator_max=0.01), but the recovery_callback rebuilds
        # the SAME game against a relaxed (trivially feasible) DTMPCParameters and hands DTMPC
        # a genuinely converged solution -- exercising the full recovery path end-to-end.
        T = Float64
        n = 2
        hard_params = make_params(n; s₀ = fill(T(10.0), n), α_fixed = fill(T(1.0), n),
                                   tighten_actuator = true, actuator_max = fill(T(0.01), n))
        easy_params = make_params(n; tighten_actuator = false)
        game = make_pdgnep_di(; n_players = n, N_win = 8)
        um = create_independent_uncertainty_models(x -> T(0.05), n)
        x0 = vcat([T[Float64(i), 0.0] for i in 1:n]...)
        N_sim = 5
        prob = DynamicTubeMPCProblem(game, x0, N_sim, hard_params, um)

        function relax_and_retry(; step, current_x, s_cur, carried_widths, carried_alphas, warmstart, aug_game, game)
            aug_easy = build_augmented_game(game, easy_params, um, current_x, s_cur;
                                             discrete = true, carried_widths = carried_widths, carried_alphas = carried_alphas)
            return DifferentialGamesBase.solve(aug_easy, make_falcon(); warmstart = warmstart, verbose = false, check_compatibility = false)
        end
        sol = solve(prob, DynamicTubeMPC(make_falcon(); recovery_callback = relax_and_retry); verbose = false)

        @test sol.solver_info[:terminated_early] == false
        @test sol.solver_info[:n_steps_simulated] == N_sim
        @test sol.converged == true
        @test all(isfinite, sol.X)
    end

    # ── display ─────────────────────────────────────────────────────────────
    @testset "Display" begin
        prob, _, _ = make_dtmpc_problem(n_players = 2, N_win = 6, N_sim = 4)
        sol = solve(prob, DynamicTubeMPC(make_falcon()))
        buf = IOBuffer()
        show(buf, prob); show(buf, sol); show(buf, DynamicTubeMPC(make_falcon()))
        show(buf, MIME("text/plain"), prob); show(buf, MIME("text/plain"), sol)
        @test length(take!(buf)) > 0
    end

end # @testset "DTMPC"
