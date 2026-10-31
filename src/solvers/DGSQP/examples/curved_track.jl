# ============================================================================
# Curved-track racing study of Zhu & Borrelli (ICRA 2023), Sec. V-A
#
# Port of `scripts/DGSQP_ALGAMES_monte_carlo_curve.py` from
# https://github.com/zhu-edward/DGSQP (MIT licence, Copyright (c) 2022 Edward Zhu).
#
# Two cars, described by the kinematic bicycle model in a Frenet frame, race through a
# curved track segment with a θ° turn. Each car trades off input effort, input rate and
# progress along the track (cost (19) of the paper), subject to input bounds and rate
# bounds, the track boundary and a collision-avoidance constraint. The study draws random
# collision-free initial states, builds an initial guess by rolling out PID lane-following
# controllers, and records whether DGSQP converges to a point satisfying the KKT conditions
# (6) within 10⁻³.
#
# Everything is written against the open-loop interface `OpenLoopGame` so that the
# previous-input terms of the paper (input-rate cost and bounds) are expressed exactly as
# in the reference implementation, with Euler discretisation.
#
# Usage:
#     include("curved_track.jl")
#     study = run_curved_track_study(; thetas = [90], horizons = [10], trials = 20)
# ============================================================================

using LinearAlgebra
using Random
using Printf
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# ── Track ───────────────────────────────────────────────────────────────────
"""
    CurveTrack(; enter_straight = 1.0, curve_length = 8.0, swept_angle = π/2,
                 exit_straight = 5.0, width = 2.0)

A straight – circular arc (counter-clockwise) – straight track of centre-line segments
`[enter_straight, 0], [curve_length, curve_length/swept_angle], [exit_straight, 0]` starting at the
origin heading along +x, in Frenet coordinates `(s, e_y)` with `e_y ∈ [-width/2, width/2]`.
"""
Base.@kwdef struct CurveTrack
    enter_straight::Float64 = 1.0
    curve_length::Float64   = 8.0
    swept_angle::Float64    = π / 2
    exit_straight::Float64  = 5.0
    width::Float64          = 2.0
end

# Cumulative segment starts, curvatures and tangent angles at the segment starts.
function _segments(t::CurveTrack)
    lens = [t.enter_straight, t.curve_length, t.exit_straight]
    curv = [0.0, t.swept_angle / t.curve_length, 0.0]
    starts = [0.0; cumsum(lens)]
    ang = [0.0]
    for i in 1:3
        push!(ang, ang[end] + lens[i] * curv[i])
    end
    return lens, curv, starts, ang
end

"""
    curvature(track, s)

Piecewise-constant centre-line curvature at arc length `s` (as in the reference `get_curvature`).
"""
function curvature(t::CurveTrack, s)
    _, curv, starts, _ = _segments(t)
    i = clamp(searchsortedlast(starts[1:end-1], s), 1, 3)
    return curv[i]
end

"""
    tangent_angle(track, s)

Centre-line tangent angle at `s`, piecewise linear through the segment ends (reference `get_tangent_angle_casadi_fn`).
"""
function tangent_angle(t::CurveTrack, s)
    _, _, starts, ang = _segments(t)
    s <= starts[1] && return ang[1] + zero(s)
    s >= starts[end] && return ang[end] + zero(s)
    i = clamp(searchsortedlast(starts, s), 1, 3)
    return ang[i] + (ang[i+1] - ang[i]) * (s - starts[i]) / (starts[i+1] - starts[i])
end

"""
    local_to_global(track, s, e_y) -> (x, y)

Cartesian position of the Frenet point `(s, e_y)` (positive `e_y` to the left of the centre-line).
"""
function local_to_global(t::CurveTrack, s::Real, ey::Real)
    lens, curv, starts, ang = _segments(t)
    # origin and heading at each segment start
    pos = [(0.0, 0.0)]
    for i in 1:3
        x, y = pos[end]; ψ = ang[i]
        if curv[i] == 0
            push!(pos, (x + lens[i] * cos(ψ), y + lens[i] * sin(ψ)))
        else
            r = 1 / curv[i]; θ = lens[i] / r
            xc, yc = x - r * sin(ψ), y + r * cos(ψ)
            push!(pos, (xc + r * sin(ψ + θ), yc - r * cos(ψ + θ)))
        end
    end
    i = clamp(searchsortedlast(starts[1:end-1], s), 1, 3)
    d = s - starts[i]
    x0, y0 = pos[i]; ψ = ang[i]
    if curv[i] == 0
        return x0 + d * cos(ψ) - ey * sin(ψ), y0 + d * sin(ψ) + ey * cos(ψ)
    else
        r = 1 / curv[i]
        xc, yc = x0 - r * sin(ψ), y0 + r * cos(ψ)
        φ = ψ - π / 2 + d / r
        return xc + (r - ey) * cos(φ), yc + (r - ey) * sin(φ)
    end
end

# ── Vehicle model ───────────────────────────────────────────────────────────
"""
Parameters of the Frenet-frame kinematic bicycle of the reference implementation
(`CasadiKinematicBicycleCombined`): wheel base `L_f = L_r = 0.13`, mass 2.366,
quadratic drag 0.1 and slip drag 0.1.
"""
Base.@kwdef struct BicycleParams
    L_f::Float64 = 0.13
    L_r::Float64 = 0.13
    mass::Float64 = 2.366
    c_dr::Float64 = 0.1
    c_s::Float64 = 0.1
end

# State q = [x, y, v, e_ψ, s, e_y], input u = [a, δ]; one Euler step of length dt.
function bicycle_step(q, u, track::CurveTrack, p::BicycleParams, dt)
    x, y, v, eψ, s, ey = q
    a, δ = u
    L = p.L_f + p.L_r
    β = atan(tan(δ) * p.L_r, L)
    ψdot = v / p.L_r * sin(β)
    F_ext = -p.c_dr * v * abs(v) - p.c_s * ψdot^2
    c = curvature(track, s)
    ψt = tangent_angle(track, s)
    dx = v * cos(β + ψt + eψ)
    dy = v * sin(β + ψt + eψ)
    dv = a + F_ext / p.mass
    deψ = ψdot - c * v * cos(β + eψ) / (1 - ey * c)
    ds = v * cos(β + eψ) / (1 - ey * c)
    dey = v * sin(β + eψ)
    return q .+ dt .* [dx, dy, dv, deψ, ds, dey]
end

# ── The game ────────────────────────────────────────────────────────────────
"""
    CurvedTrackSetup(; track, N, dt = 0.1, ...)

Constants of the study: input bounds `|a| ≤ 2.1`, `|δ| ≤ 0.436`, rate bounds `|Δa| ≤ 10 dt`, `|Δδ| ≤ 4.5 dt`,
track half-width `1`, car radius `0.2`, cost weights `R = R_d = I`, `(w_p, w_c) = (10, 5)`.
"""
Base.@kwdef struct CurvedTrackSetup
    track::CurveTrack = CurveTrack()
    N::Int = 10
    dt::Float64 = 0.1
    bicycle::BicycleParams = BicycleParams()
    u_max::Vector{Float64} = [2.1, 0.436]
    du_max::Vector{Float64} = [10.0, 4.5]
    half_width::Float64 = 1.0
    car_radius::Float64 = 0.2
    input_weight::Vector{Float64} = [1.0, 1.0]
    rate_weight::Vector{Float64} = [1.0, 1.0]
    w_progress::Float64 = 10.0
    w_comp::Float64 = 5.0
end

# Joint Euler rollout of both cars from the initial states `q0 = [q¹; q²]`.
function _rollout_cars(S::CurvedTrackSetup, q0::AbstractVector, U1, U2)
    T = promote_type(eltype(q0), eltype(U1))
    q1 = Vector{T}(q0[1:6]); q2 = Vector{T}(q0[7:12])
    Q1 = [q1]; Q2 = [q2]
    for k in 1:S.N
        q1 = bicycle_step(q1, U1[:, k], S.track, S.bicycle, S.dt)
        q2 = bicycle_step(q2, U2[:, k], S.track, S.bicycle, S.dt)
        push!(Q1, q1); push!(Q2, q2)
    end
    return Q1, Q2
end

"""
    curved_track_game(setup, q0) -> OpenLoopGame

The two-car game of the paper from the joint initial state `q0` (12 numbers). Decision vector
`w = [vec(U¹); vec(U²)]`, `Uⁱ = [aᵢ,ₖ; δᵢ,ₖ] ∈ ℝ^{2×N}`; previous input `u₋₁ = 0`.
"""
function curved_track_game(S::CurvedTrackSetup, q0::AbstractVector)
    N = S.N
    split(w) = reshape(w[1:2N], 2, N), reshape(w[2N+1:4N], 2, N)
    dU(U) = U .- hcat(zeros(eltype(U), 2), U[:, 1:end-1])      # u_k − u_{k−1}, u_{−1} = 0

    function cost(i)
        return function (w)
            Us = split(w)
            Q1, Q2 = _rollout_cars(S, q0, Us[1], Us[2])
            Qi, Qj = i == 1 ? (Q1, Q2) : (Q2, Q1)
            U = Us[i]; ΔU = dU(U)
            J = 0.5 * sum(S.input_weight[r] * U[r, k]^2 + S.rate_weight[r] * ΔU[r, k]^2 for r in 1:2, k in 1:N)
            si, sj = Qi[end][5], Qj[end][5]
            return J - S.w_progress * si + S.w_comp * atan(sj - si)
        end
    end

    function constraints(w)
        Us = split(w)
        Q1, Q2 = _rollout_cars(S, q0, Us[1], Us[2])
        c = eltype(w)[]
        for i in 1:2
            U = Us[i]; ΔU = dU(U); Q = i == 1 ? Q1 : Q2
            for k in 1:N
                append!(c, [U[1, k] - S.u_max[1], -U[1, k] - S.u_max[1], U[2, k] - S.u_max[2], -U[2, k] - S.u_max[2]])
                append!(c, [ΔU[1, k] - S.dt * S.du_max[1], -ΔU[1, k] - S.dt * S.du_max[1],
                            ΔU[2, k] - S.dt * S.du_max[2], -ΔU[2, k] - S.dt * S.du_max[2]])
                ey = Q[k+1][6]
                append!(c, [ey - S.half_width, -ey - S.half_width])
            end
        end
        d_obs = 2 * S.car_radius
        for k in 1:N
            p1, p2 = Q1[k+1][1:2], Q2[k+1][1:2]
            push!(c, d_obs^2 - sum(abs2, p1 .- p2))
        end
        return c
    end

    return OpenLoopGame([2, 2], N, Function[cost(1), cost(2)], constraints)
end

# ── Warm start: PID lane followers ─────────────────────────────────────────
"""
    pid_warm_start(setup, q0_i, ref_lane, ref_speed) -> (U, Q)

Roll out a PID lane follower (steering PID on `5(e_y − ref) + e_ψ` with `Kp = 1`, `Ki = 0.005`; speed PID with
`Kp = 1`; rate- and magnitude-saturated as in the reference) for `N` steps. Returns the inputs `U` (`2×N`) and
states `Q` (`6×(N+1)`).
"""
function pid_warm_start(S::CurvedTrackSetup, q0::AbstractVector, ref_lane::Real, ref_speed::Real)
    N = S.N
    q = Vector{Float64}(q0)
    Q = zeros(6, N + 1); Q[:, 1] = q
    U = zeros(2, N)
    ei = 0.0; u_prev_s = 0.0; u_prev_a = 0.0
    for k in 1:N
        # speed PID: error v − v_ref, Kp = 1; rate then magnitude saturation (reference `PID.solve`)
        a = -(1.0 * (q[3] - ref_speed))
        a = clamp(a - u_prev_a, -S.du_max[1], S.du_max[1]) + u_prev_a
        a = clamp(a, -S.u_max[1], S.u_max[1])
        # steering PID on the combined lane/heading error, Kp = 1, Ki = 0.005
        e = 5.0 * (q[6] - ref_lane) + 1.0 * q[4]
        ei += e * S.dt
        δ = -(1.0 * e + 0.005 * ei)
        δ = clamp(δ - u_prev_s, -S.du_max[2], S.du_max[2]) + u_prev_s
        δ = clamp(δ, -S.u_max[2], S.u_max[2])
        u_prev_a, u_prev_s = a, δ
        U[:, k] = [a, δ]
        q = bicycle_step(q, [a, δ], S.track, S.bicycle, S.dt)
        Q[:, k+1] = q
    end
    return U, Q
end

# ── Monte Carlo study ───────────────────────────────────────────────────────
# Random collision-free initial states, as in the reference script.
function sample_initial_state(S::CurvedTrackSetup, rng::AbstractRNG)
    obs_d = 2 * S.car_radius
    first_seg = S.track.enter_straight
    while true
        s1 = max(0.1, rand(rng) * first_seg)
        e1 = rand(rng) * S.half_width * 2 - S.half_width
        v1 = rand(rng) + 2
        d = 2π * rand(rng)
        s2 = s1 + 1.2 * obs_d * cos(d)
        s2 < 0 && continue
        e2 = e1 + 1.2 * obs_d * sin(d)
        abs(e2) > S.half_width && continue
        v2 = rand(rng) + 2
        x1, y1 = local_to_global(S.track, s1, e1)
        x2, y2 = local_to_global(S.track, s2, e2)
        q1 = [x1, y1, v1, 0.0, s1, e1]
        q2 = [x2, y2, v2, 0.0, s2, e2]
        # PID initial guess must itself be collision free (as in the reference)
        U1, Q1 = pid_warm_start(S, q1, e1, v1)
        U2, Q2 = pid_warm_start(S, q2, e2, v2)
        any(k -> norm(Q1[1:2, k] - Q2[1:2, k]) < obs_d, 1:S.N+1) && continue
        return vcat(q1, q2), join_controls([U1, U2])
    end
end

"""
    run_curved_track_study(; thetas = [45, 75, 90], horizons = [10, 15, 20, 25], trials = 200, seed = 1,
                             solver = dgsqp_paper_options(), verbose = false) -> Vector{NamedTuple}

Monte Carlo study of Sec. V-A (Table I). For each turn angle `θ` (degrees) and horizon `N` it draws `trials`
random initial conditions and records, for DGSQP: `converged` (KKT within `p_tol`/`d_tol`), `fail` (diverged or
QP failure), `max_iters` (iteration limit reached), the mean number of QP solves and wall times.
The 200-trial, twelve-setting study of the paper takes a long time; use fewer `trials` for a smoke test.
"""
function run_curved_track_study(; thetas = [45, 75, 90], horizons = [10, 15, 20, 25], trials = 200, seed = 1,
                                  solver::DGSQP = dgsqp_paper_options(), verbose::Bool = false)
    rng = MersenneTwister(seed)
    rows = NamedTuple[]
    for θ in thetas, N in horizons
        S = CurvedTrackSetup(; track = CurveTrack(; swept_angle = θ * π / 180), N = N)
        conv = fail = maxit = 0
        qps = Int[]; times = Float64[]
        for t in 1:trials
            q0, w0 = sample_initial_state(S, rng)
            game = curved_track_game(S, q0)
            res = dgsqp_solve(game, w0; solver = solver)
            if res.converged
                conv += 1; push!(qps, res.qp_solves); push!(times, res.time)
            elseif res.status in (:diverged, :qp_fail)
                fail += 1
            else
                maxit += 1
            end
            verbose && @printf("θ=%d N=%d trial %d: %s (%d QP)\n", θ, N, t, res.status, res.qp_solves)
        end
        push!(rows, (theta = θ, N = N, trials = trials, converged = conv, fail = fail, max_iters = maxit,
                     mean_qp_solves = isempty(qps) ? NaN : sum(qps) / length(qps),
                     mean_time = isempty(times) ? NaN : sum(times) / length(times)))
    end
    return rows
end
