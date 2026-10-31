# ============================================================================
# Head-to-head racing through the first hairpin of the Austin track (Scenario 3 of
# Zhu & Borrelli, arXiv:2404.00186; the "very sharp corner" example of the FALCON suite)
#
# Two cars described by the 7-state dynamic bicycle model with Pacejka tyres (1/10 scale RC car
# parameters), discretised with RK4 and ten substeps per 0.1 s step as in the reference
# implementations, race through a 30 m segment of the Circuit of the Americas that contains the
# first hairpin. Each car minimises input effort and rate, a contouring/lag cost that ties the
# progress state `s` to its position, and a competitive terminal progress term;
# constraints are control bounds, track boundaries (hard, with the variable track width) and
# bounding-circle collision avoidance.
#
# This is the problem of `f1_austin_racing.jl` / `visualize_f1_austin.jl` in the FALCON repository
# (N = 70, Δt = 0.1, segment s ∈ [60, 90] m, q_l = 20), written as a `SeparableRolloutGame`
# and solved with the DGSQP solver in its racing-paper configuration.
#
# Requires DifferentialGamesBenchmarks (track data, dynamics) next to DifferentialGamesBaseSolvers.
#
#     include("austin_hairpin.jl")
#     game, w0, info = austin_hairpin_game(; N = 70)
#     res = dgsqp_solve(game, w0; solver = dgsqp_racing_options(), verbose = true)
# ============================================================================

using LinearAlgebra
using Printf
using DifferentialGamesBase
using DifferentialGamesBaseSolvers
using DifferentialGamesBenchmarks

const W_RATE_HAIRPIN = [1.0, 1.0, 1e-4]     # input-rate weights of the FALCON racing script

"""
    hairpin_track() -> RaceTrack

The 30 m segment `s ∈ [60, 90]` of the Austin tenth-scale track (re-based to 0) that contains the first hairpin.
"""
hairpin_track() = load_race_track(joinpath(pkgdir(DifferentialGamesBenchmarks), "data", "f1_austin_hairpin_60_90.csv"))

# One control step: RK4 with M substeps of the continuous bicycle model.
function bicycle_rk4_step(f, dt, M)
    h = dt / M
    return function (x, u)
        for _ in 1:M
            k1 = f(x, u, nothing, 0.0)
            k2 = f(x .+ (h / 2) .* k1, u, nothing, 0.0)
            k3 = f(x .+ (h / 2) .* k2, u, nothing, 0.0)
            k4 = f(x .+ h .* k3, u, nothing, 0.0)
            x = x .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        end
        return x
    end
end

# Arc length of the centre-line point closest to (x, y), searched within ±1.5 m of the guess `s_guess`.
function project_to_centerline(track::RaceTrack, x, y, s_guess)
    ss = range(max(0.0, s_guess - 1.5), min(track.length, s_guess + 1.5); length = 151)
    return ss[argmin([(x - track.x(s))^2 + (y - track.y(s))^2 for s in ss])]
end

# Closed-loop centre-line plan (curvature feed-forward, lateral/heading feedback, speed control) simulated with the
# same step map; returns the controls it applied (`3×N`).
function hairpin_warm_start(track::RaceTrack, P::RacingParameters, step, q0; v_ref = 2.5, ey_cap = 0.35, ey_fade = 10, brake_factor = 0.5)
    L = P.l_f + P.l_r
    κ(s) = begin
        h = 1e-3; s1, s2 = max(s - h, 0.0), min(s + h, track.length)
        t1, t2 = track_tangent(track, s1), track_tangent(track, s2)
        (t1[1] * t2[2] - t1[2] * t2[1]) / (s2 - s1)
    end
    _, _, ey0 = track_errors(track, q0[4], q0[5], q0[7])
    ey_hold = clamp(ey0, -ey_cap, ey_cap)       # keep at most `ey_cap` from the centre-line, as the track narrows in the hairpin
    q = copy(q0); U = zeros(3, P.N)
    for k in 1:P.N
        s = project_to_centerline(track, q[4], q[5], q[7])      # true projection, not the integrated progress state
        tx, ty = track_tangent(track, s)
        _, _, ey = track_errors(track, q[4], q[5], s)
        eψ = mod(q[6] - atan(ty, tx) + π, 2π) - π
        ey_ref = ey0 + (ey_hold - ey0) * min(1.0, (k - 1) / ey_fade)
        γ = clamp(atan(κ(s) * L) - 0.5 * (ey - ey_ref) - 1.0 * eψ, -P.u_max[2], P.u_max[2])
        κl = maximum(abs(κ(clamp(s + d, 0.0, track.length))) for d in 0.0:0.5:3.0)
        v_tgt = min(v_ref, max(0.8, brake_factor * sqrt(P.μ * 9.81 / (κl + 1e-3))))
        a = clamp(3.0 * (v_tgt - q[1]), P.u_min[1], P.u_max[1])
        u = [a, γ, 0.0]
        q = step(q, u)
        s_new = project_to_centerline(track, q[4], q[5], s)
        u[3] = clamp((s_new - s) / P.dt, P.u_min[3], P.u_max[3])    # progress rate that keeps s on the projection
        U[:, k] = u
        q[7] = s_new
    end
    return U
end



"""
    single_agent_plan(track, P, step, q0, U0; u_prev0) -> U

Each car solves its own problem alone, as the warm start of the FALCON script does: input effort and rate, the contouring/lag
cost, progress `−w_p s`, control bounds and the hard track boundaries, no opponent. Solved with DGSQP from the controller plan `U0`
(`3×N`); falls back to `U0` if the solve fails.
"""
function single_agent_plan(track::RaceTrack, P::RacingParameters, step, q0, U0; iters = 150)
    N = P.N
    tcost(q) = ((ec, el, _) = track_errors(track, q[4], q[5], q[7]); 0.5 * (P.q_c * ec^2 + P.q_l * el^2))
    cost = function (i, k, X, U, Up)
        k < N ? 0.5 * sum(P.w_input .* U[1] .^ 2) + 0.5 * sum(W_RATE_HAIRPIN .* (U[1] .- Up[1]) .^ 2) + tcost(X[1]) :
                -P.w_prog * X[1][7] + tcost(X[1])
    end
    cons = function (k, X, U, Up)
        c = eltype(X[1])[]
        k < N && (append!(c, U[1] .- P.u_max); append!(c, P.u_min .- U[1]))
        if k >= 1
            _, _, ey = track_errors(track, X[1][4], X[1][5], X[1][7])
            append!(c, [ey - track.left(X[1][7]), -ey - track.right(X[1][7])])
        end
        return c
    end
    g1 = SeparableRolloutGame([3], [7], N, [q0], [step], cost, cons; u_prev0 = [U0[:, 1]])
    r = dgsqp_solve(g1, vec(U0); solver = dgsqp_racing_options(; sqp_iters = iters, p_tol = 1e-4, d_tol = 2e-2))
    return r.status in (:diverged, :qp_fail) ? U0 : reshape(r.w, 3, N)
end

"""
    austin_hairpin_game(; N = 70, dt = 0.1, substeps = 10, q_l = 20.0, v_ref = 2.5) -> (game, w0, info)

The two-car hairpin racing game as a `SeparableRolloutGame`, an initial decision vector `w0` (PID-like centre-line
plans for both cars) and a `NamedTuple` with the track, parameters, step map and initial states.
"""
function austin_hairpin_game(; N::Int = 70, dt::Float64 = 0.1, substeps::Int = 10, q_l::Float64 = 20.0,
                               v_ref = 2.5, track::RaceTrack = hairpin_track(), q0 = nothing,
                               collision_aware::Bool = false, single_agent::Bool = false, warm_start = nothing)
    P = RacingParameters(; N = N, dt = dt, q_l = q_l)
    f = bicycle_dynamics(P)
    step = bicycle_rk4_step(f, dt, substeps)
    q0 = q0 === nothing ? collect(racing_starts(track, P)) : [Vector{Float64}(q0[1]), Vector{Float64}(q0[2])]
    d_min2 = racing_collision_distance(P)^2

    track_cost(q) = ((ec, el, _) = track_errors(track, q[4], q[5], q[7]); 0.5 * (P.q_c * ec^2 + P.q_l * el^2))

    cost = function (i, k, X, U, Up)
        j = 3 - i
        if k < N
            u, up = U[i], Up[i]
            return 0.5 * sum(P.w_input .* u .^ 2) + 0.5 * sum(W_RATE_HAIRPIN .* (u .- up) .^ 2) + track_cost(X[i])
        end
        return -P.w_prog * X[i][7] + P.w_comp * (X[j][7] - X[i][7]) + track_cost(X[i])
    end

    constraints = function (k, X, U, Up)
        T = eltype(X[1])
        c = T[]
        if k < N                                   # control bounds on u_k
            for i in 1:2
                append!(c, U[i] .- P.u_max)
                append!(c, P.u_min .- U[i])
            end
        end
        if k >= 1                                  # track boundary and collision on x_k
            for i in 1:2
                _, _, ey = track_errors(track, X[i][4], X[i][5], X[i][7])
                append!(c, [ey - track.left(X[i][7]), -ey - track.right(X[i][7])])
            end
            push!(c, d_min2 - ((X[1][4] - X[2][4])^2 + (X[1][5] - X[2][5])^2))
        end
        return c
    end

    if warm_start !== nothing         # controls supplied by the caller, e.g. the single-agent plans of the FALCON benchmark
        U1, U2 = Matrix{Float64}(warm_start[1]), Matrix{Float64}(warm_start[2])
        game = SeparableRolloutGame([3, 3], [7, 7], N, q0, [step, step], cost, constraints; u_prev0 = [U1[:, 1], U2[:, 1]])
        return game, join_controls([U1, U2]), (track = track, params = P, step = step, q0 = q0)
    end
    vr = v_ref isa Number ? [Float64(v_ref), Float64(v_ref)] : collect(Float64, v_ref)
    U1 = hairpin_warm_start(track, P, step, q0[1]; v_ref = vr[1])
    U2 = hairpin_warm_start(track, P, step, q0[2]; v_ref = vr[2])
    if collision_aware
        # Slow the car behind until the two plans do not collide (a crude stand-in for the reference's
        # collision-aware CA-LTV-MPC warm start).
        plan_separation(Ua, Ub) = begin
            xa, xb = q0[1], q0[2]; d = hypot(xa[4] - xb[4], xa[5] - xb[5])
            for k in 1:N
                xa = step(xa, Ua[:, k]); xb = step(xb, Ub[:, k])
                d = min(d, hypot(xa[4] - xb[4], xa[5] - xb[5]))
            end
            d
        end
        behind = q0[1][7] <= q0[2][7] ? 1 : 2
        for _ in 1:12
            plan_separation(U1, U2) >= sqrt(d_min2) && break
            vr[behind] = max(0.3, vr[behind] - 0.25)
            U1 = behind == 1 ? hairpin_warm_start(track, P, step, q0[1]; v_ref = vr[1]) : U1
            U2 = behind == 2 ? hairpin_warm_start(track, P, step, q0[2]; v_ref = vr[2]) : U2
        end
    end
    if single_agent
        U1 = single_agent_plan(track, P, step, q0[1], U1)
        U2 = single_agent_plan(track, P, step, q0[2], U2)
    end
    game = SeparableRolloutGame([3, 3], [7, 7], N, q0, [step, step], cost, constraints;
                                u_prev0 = [U1[:, 1], U2[:, 1]])      # the FALCON script's u₋₁ is the warm start's first control
    return game, join_controls([U1, U2]), (track = track, params = P, step = step, q0 = q0)
end

"""
    hairpin_metrics(game, w, info) -> NamedTuple

Quantities recomputed from the control sequences, independent of the solver's flags: final progress of each car,
smallest car separation, smallest distance of a car centre to a track edge (negative = off track), maximum
distance between a car and the centre-line point at its own progress.
"""
function hairpin_metrics(game::SeparableRolloutGame, w, info)
    N = game.N
    U = split_controls(game, w)
    Q = [Vector{Vector{Float64}}(undef, N + 1) for _ in 1:2]
    for i in 1:2
        Q[i][1] = game.x0[i]
        for k in 1:N
            Q[i][k+1] = info.step(Q[i][k], U[i][:, k])
        end
    end
    tr = info.track
    sep = minimum(hypot(Q[1][k][4] - Q[2][k][4], Q[1][k][5] - Q[2][k][5]) for k in 1:N+1)
    clear = minimum(begin
        _, _, ey = track_errors(tr, q[4], q[5], q[7]); min(tr.left(q[7]) - ey, tr.right(q[7]) + ey)
    end for i in 1:2 for q in Q[i])
    lag = maximum(hypot(q[4] - tr.x(q[7]), q[5] - tr.y(q[7])) for i in 1:2 for q in Q[i])
    return (progress = (Q[1][end][7], Q[2][end][7]), min_separation = sep, track_clearance = clear, max_lag = lag,
            states = Q)
end
