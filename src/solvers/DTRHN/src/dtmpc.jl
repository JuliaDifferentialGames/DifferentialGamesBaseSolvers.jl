# ============================================================================
# solvers/DTRHN/src/dtmpc.jl
#
# Per-Agent Dynamic Tube MPC (DTMPC) wrapper around a constrained PD-GNEP
# game solver (validated with ALGAMES; faithful path validated with FALCON).
#
# Follows the "Robust Dynamic-Tube FALCON" theory spine (v0.2) for the tube
# STRUCTURE, and the later unified manuscript (Outland & Arya 2026-07-23),
# Problem 2 / Eq. (9)-(10), for the tube BANDWIDTH — which the Spine and the
# manuscript disagree about, and where the manuscript is authoritative.
#
# Each agent `i` carries an ORDERED list of tube CHANNELS (translational,
# optionally rotational, optionally detection). Per channel `k`:
#
#   - width state  s_k  (scalar)
#   - bandwidth state α_k (scalar) with its own artificial control v_k:
#         α_k,{t+1} = α_k,t + dt · v_k,t          (manuscript Eq. 10, α̇=v)
#     bounded α_k ∈ [α_lo, α_max], |v_k| ≤ v_max. α is therefore CO-OPTIMIZED
#     with the trajectory, not a fixed design constant.
#       * α_lo > 0 is ESSENTIAL, not cosmetic: the tube ODE is unstable for
#         α ≤ 0 and the exact recursion below divides by α. A prior version of
#         this file that omitted the bound had the optimizer drive α negative
#         and s blow up.
#   - width recursion (equality; the tight/deterministic instance of the
#     theory's inequality (6), valid by Remark 1.1), now with the per-node α:
#         ρ_k,t     = exp(−α_k,t · dt)
#         s_k,{t+1} = ρ_k,t · s_k,t + (β_k/α_k,t)·(1 − ρ_k,t)
#     This is the exact integral of ṡ=−αs+β over one step. It is smooth but
#     genuinely NONLINEAR in the decision variable α (removable-singular as
#     α→0, where it tends to β·dt — `α_lo` is the numerical guard).
#   - NO objective curvature (Spine Thm 1): the augmented cost is the ORIGINAL
#     cost, unchanged. Manuscript Problem 2's objective is likewise J^i alone —
#     no w_Φ s² + w_α α² + w_v v² terms (unlike an older removed API). `s_k`,
#     `α_k` and `v_k` are regularized only by FALCON's own proximal term.
#
# Augmented layout (α appended AFTER all widths — this ordering is mandatory,
# not stylistic: it keeps the PHYSICAL block leading and every existing `s`
# index unchanged, which is what lets `build_augmented_cost`, the shared-
# constraint wrapper, the collision/containment tightenings and the games'
# own defensive `length(x) ≥ N ? x[N] : 0` tube reads survive untouched. An
# INTERLEAVED layout would silently make those read an α as if it were a width):
#     x̃ᵢ = [ xᵢ(nᵢ) ; s₁…s_nch ; α₁…α_nch ]
#     ũᵢ = [ uᵢ(mᵢ) ; v₁…v_nch ]
#
# The genuine "recursive-feasibility" freedom the theory needs (Lemma 3.1(b) /
# Prop. 3) lives ENTIRELY at each receding-horizon WINDOW's initial node: the
# initial pair (x̄_i,{t|t}, s_i,{t|t}) is a DECISION VARIABLE constrained only
# by initial containment ‖x_i,t − x̄_i,{t|t}‖ ≤ s_i,{t|t} (§4) — NOT a hard
# re-initialization. This requires an inner solver that supports FREEING the
# initial-state pin for those components (e.g. FALCON's `free_initial` field,
# via `dtmpc_free_initial_mask`); solvers without this notion (e.g. ALGAMES)
# keep the classical hard pin, for which the containment constraint is added
# but trivially satisfied (a valid, more conservative special case).
#
# Ancillary (execution-time) boundary-layer sliding controller, ONE per
# translational/rotational channel, gain K_k = α_k(t) · s_k(t) — now using the
# CARRIED per-step α from the closed loop, not a design constant:
#   translational: a = a* − K_t·sat(s̃/s_t),  s̃ = ṽ+λr̃
#   rotational:    ω = ω* − K_r·sat(p̃/s_r)
# Detection channel has NO ancillary (the accumulator isn't actuated).
#
# Tightening (§3): actuator/rate authority reserves K_k = α_k·s_k.
#   COST OF CO-OPTIMIZING α, stated plainly: with α a decision variable the
#   actuator reserve α_k·s_k is BILINEAR, so `‖a*‖ + α·s ≤ a_max` is no longer
#   convex and is marked `is_convex=false` — it moves off FALCON's hard,
#   always-linearized, no-AL-slack path onto the AL-penalized path. This is
#   precisely the trade `On-Orbit-Jamming-Game/docs/theory.md` §11 flagged as
#   the reason to keep α fixed; it was accepted deliberately when the
#   manuscript was made authoritative. If AL convergence degrades where
#   `tighten_actuator=true`, the intended fallback is a CONVEX conservative
#   surrogate α·s ≤ (α²/c + c·s²)/2 (exact at c = α/s) — NOT silently
#   reverting to a fixed α.
#   The rate tightening `‖ω*‖ + s_r ≤ rate_max` is α-FREE and stays convex, as
#   does detection `s_det+s_d≤limit`. Pointing-cone/link-floor/collision
#   tightenings are nonlinear geometric margins (Eq. 5's general Lipschitz
#   form) and were already AL-penalized (`is_convex=false`).
#
# References:
#   Lopez, Slotine, How (2019), "Dynamic Tube MPC for Nonlinear Systems" (ACC).
#   Mayne, Seron, Raković (2005) — free-initial-nominal-state RH-MPC convention.
#   Outland & Arya (2026), FALCON — PD-GNEP base algorithm (external);
#   Robust Dynamic-Tube FALCON theory spine v0.2 — extracted to
#     `On-Orbit-Jamming-Game/docs/theory_spine_v0.2.md` (it previously existed only in a
#     chat transcript, so this citation resolved to nothing).
#   Outland & Arya (2026-07-23), "Robust Game-theoretic Motion Planning over Extended Time
#     Horizons" — Problem 2 / Eq. (9)-(10). AUTHORITATIVE over the Spine on the tube
#     bandwidth: α is a CO-OPTIMIZED STATE (α̇=v, α∈[α_lo,α_max], |v|≤v_max), not the
#     Spine's fixed contraction ρ. See the α-refactor note in the layout section below.
#
# Solvers package imports already in scope:
#   GameSolver, solve, _solve, WarmstartData, solver_capabilities,
#   _shift_warmstart, _shift_feedback_strategy, _shift_open_loop_strategy
#   + all DifferentialGamesBase exports.
# ============================================================================

# ============================================================================
# Tube-channel specification
# ============================================================================

"""
    TubeChannel{T}

One scalar dynamic-tube channel for a single agent (Theory Spine v0.2, §1/§2).

# Fields
- `kind::Symbol`        : `:translational`, `:rotational`, or `:detection`
- `s₀::T`               : initial tube width
- `α₀::T`               : INITIAL bandwidth. α is a CO-OPTIMIZED STATE
  (manuscript Eq. 10, `α̇=v`), so this seeds `α_{t=0}` only; thereafter the
  solver chooses it subject to `α ∈ [α_lo, α_max]`, `|v| ≤ v_max`. The
  contraction is per-node, `ρ_t = exp(-α_t·dt)`, not a design constant.
- `D::T`, `η::T`        : constant disturbance bound / robustness margin
  (own-agent data only, per separability (S))
- `κ_speed::T`          : optional own-state-dependent slope: the majorant
  `δ̂(x̄,ū) = D + η + κ_speed·smooth_l1(v̄)` where `v̄` is the channel's own
  velocity-like sub-block (translational: `x̄[pos_dim+1:2pos_dim]`; zero for
  channels without a natural "own speed", e.g. detection/rotational — set
  `κ_speed=0` to disable). This realizes an own-state-dependent `δ_i`
  (§11 example) as an AFFINE-in-the-nonlinear-smooth-norm majorant; SCP
  re-linearizes it every iterate like any other nonlinear own-agent quantity
  (the standard SCvx exactness caveat, obligation V-3a).
- `state_offset::Int`, `state_dim::Int` : 0-based slice of the agent's
  PHYSICAL state this channel governs (translational: `[0,2pos_dim)`;
  rotational: `[rot_state_offset, rot_state_offset+rot_dim)`; detection:
  `[det_state_offset, det_state_offset+1)`).
- `ctrl_offset::Int`, `ctrl_dim::Int`   : 0-based slice of the agent's
  PHYSICAL control this channel's ancillary corrects (translational:
  `[0,pos_dim)`; rotational: `[rot_ctrl_offset,rot_ctrl_offset+rot_dim)`;
  detection: unused, `ctrl_dim=0`, no ancillary).
"""
struct TubeChannel{T}
    kind::Symbol
    s₀::T
    α₀::T
    D::T
    η::T
    κ_speed::T
    state_offset::Int
    state_dim::Int
    ctrl_offset::Int
    ctrl_dim::Int
end

"""
    ρ_of(α, dt) -> ρ

Per-node contraction `ρ_t = exp(-α_t·dt)`. `α` is now a DECISION VARIABLE read from the
augmented state, so this is evaluated per node rather than once per channel (it replaces the
old `ρ_fixed(c, dt)`, which closed over a design constant).
"""
@inline ρ_of(α, dt::Real) = exp(-α * dt)

# Smooth positive guard on α before dividing by it. α is hard-bounded below by `α_lo > 0` via
# StateBounds, so this should never bind on a feasible iterate; it exists because dividing by a
# near-zero α is catastrophic and an SCP trial point is not guaranteed to sit inside the box.
# √(α²+ε²) is C^∞ (unlike max(α,α_lo), whose kink would poison the linearization) and ≈α
# whenever α ≫ ε.
@inline _α_guard(α, α_lo::Real) = sqrt(α * α + (α_lo / 10)^2)

# Continuous-time RATE bound β(x̄,ū) = Δᵢ + D + η + κ_speed·smooth_l1(v̄), matching
# LSH's β = Δ + D + η (model-error bound Δ + disturbance bound D + robustness margin
# η — all continuous-time RATES, e.g. m/s²-equivalent, not yet integrated over a step).
# `Δ_i` is the SUPPLIED per-agent uncertainty model (own-state, separable — §S),
# evaluated by the caller and passed in as a scalar (zero if the channel has none, e.g.
# rotational/detection channels pass Δ_i≡0, their bound is pure D_r/D_d+η). `v̄` is the
# (own-state) sub-block `x̄[state_offset+1:state_offset+state_dim]` treated as a
# velocity-like quantity when κ_speed≠0 (translational channel).
function _tube_beta(c::TubeChannel{T}, x_i::AbstractVector, u_i::AbstractVector, Δ_i) where {T}
    base = Δ_i + c.D + c.η
    c.κ_speed == zero(T) && return base
    v̄ = @view x_i[c.state_offset+1 : c.state_offset+c.state_dim]
    ε = T(1e-6)
    l1 = sum(sqrt(v̄[j]^2 + ε^2) - ε for j in eachindex(v̄))
    return base + c.κ_speed * l1
end

"""
    tube_delta_hat(c, x_i, u_i, Δ_i, dt, α, α_lo) -> δ̂

The DISCRETE-step additive term matching the EXACT solution of the continuous
tube ODE `ṡ=−α·s+β` sampled over one step `dt` (LSH Eq. 13/27):
`s_{k+1} = ρ·s_k + (β/α)·(1−ρ)`, `ρ=e^{−α·dt}`. This is the "tight" (equality)
instance of the theory's inequality (6), valid by Remark 1.1 — NOT simply `β`
itself (which would under-integrate whenever `α·dt` is not small).

`α` is now passed IN per node (a decision variable read from the augmented state), rather than
taken from the channel as a design constant. `α_lo` supplies the smooth positive guard for the
division — see [`_α_guard`](@ref).
"""
function tube_delta_hat(c::TubeChannel{T}, x_i::AbstractVector, u_i::AbstractVector, Δ_i,
                        dt::Real, α, α_lo::Real) where {T}
    β  = _tube_beta(c, x_i, u_i, Δ_i)
    αg = _α_guard(α, α_lo)
    ρ  = ρ_of(αg, dt)
    return (β / αg) * (one(T) - ρ)
end

# ============================================================================
# Sampled-data tube recursion (session correction, TRANSLATIONAL channel only).
#
# The continuous-time recursion above (`ρ_of`/`tube_delta_hat`, exact for ṡ=-αs+β sampled
# on a ZOH interval) assumes the ancillary correction acts continuously between samples.
# The actual executed system is sampled-data: the ancillary boundary-layer correction
# (`lsh_boundary_layer_control`) is evaluated ONCE per physical step and held, the truth
# dynamics are rolled one discrete step, and the disturbance is injected as an ENDPOINT
# velocity kick (Δv=w·dt), not a continuously-acting perturbation. For a sampled linear
# system x_{k+1}=(1-α·h)x_k+d_k inside the boundary layer, the correct one-step contraction
# is |1-α·h|, not exp(-α·h) -- these differ substantially whenever α·h is not small
# (α=0.1, h=15s: exp(-1.5)=0.223 vs |1-1.5|=0.5), and DIVERGE (|1-α·h|>1) once α·h>2 -- the
# discrete-time instability boundary this sampled recursion imposes on α (see `α_max`'s own
# docstring update in straight_line_wolf_dtp).
"""
    sampled_tube_rho(alpha, dt) -> rho

Sampled-data (ZOH) one-step contraction magnitude `|1-alpha*dt|` for the translational
sliding-surface tube, replacing the continuous-time `exp(-alpha*dt)` for that channel only.
Stable (`rho<1`) iff `0 < alpha*dt < 2`.
"""
@inline sampled_tube_rho(alpha, dt::Real) = abs(one(alpha) - alpha * dt)

"""
    sampled_tube_increment(beta, dt) -> h*beta

Sampled-data one-step additive disturbance/margin term for the translational tube: the
rate bound `beta` (model error + disturbance bound + robustness margin, same quantity
`_tube_beta` already computes) integrated over one ZOH interval `dt` -- NOT divided by
alpha (that division was specific to the continuous exact-integral solution).
"""
@inline sampled_tube_increment(beta, dt::Real) = beta * dt

"""
    sampled_tube_next(phi, alpha, beta, dt) -> phi_{k+1}

One-step sampled-data translational tube-width recursion
`phi_{k+1} = |1-alpha*dt|*phi_k + dt*beta`, replacing the continuous
`rho_of(alpha,dt)*phi_k + tube_delta_hat(...)` for the translational channel only
(rotational/detection channels are UNCHANGED, still using the continuous recursion --
this correction is scoped exactly to the channel the componentwise sliding-surface
ancillary governs). No `_α_guard` division-by-alpha singularity here (unlike the
continuous recursion), so the raw (unguarded) `alpha` is used directly.
"""
@inline function sampled_tube_next(phi, alpha, beta, dt::Real)
    return sampled_tube_rho(alpha, dt) * phi + sampled_tube_increment(beta, dt)
end

# ============================================================================
# Parameters
# ============================================================================

"""
    DTMPCParameters{T}

Per-agent, per-channel Dynamic Tube MPC parameters (Theory Spine v0.2).

Every per-agent quantity is a length-`n_agents` vector; every per-agent field
suffixed `_r`/`_d` configures the OPTIONAL rotational/detection channel (only
meaningful where `rotational_tube[i]`/`detection_tube[i]` is true — per-agent,
so e.g. a jammer can carry 3 channels while a target carries 2).

# Core (translational) channel
- `s₀`, `α_fixed`, `D`, `η`, `κ_speed` : see [`TubeChannel`](@ref)
- `pos_dim` : leading state components treated as position (channel governs
  `[r;v]`, the first `2·pos_dim` physical state components)

# Optional rotational channel
- `rotational_tube::Vector{Bool}`, `s₀_r`, `α_fixed_r`, `D_r`, `η_r`
- `rot_dim`, `rot_state_offset`, `rot_ctrl_offset`

# Optional detection channel (no ancillary; pure tightening margin)
- `detection_tube::Vector{Bool}`, `s₀_d`, `α_fixed_d`, `D_d`, `η_d`
- `det_state_offset`, `det_limit`

# Tightening (Eq. 4/5)
- `tighten_actuator`, `actuator_max`       : `‖a*‖ + α·s ≤ a_max` (**nonconvex**:
  bilinear now that α is co-optimized — see the module header)
- `tighten_actuator_convex::Bool` : when `true`, replaces the bilinear `‖a*‖+scale·α·s≤a_max`
  tightening with the CONVEX conservative surrogate `‖a*‖+scale·(α²/c+c·s²)/2≤a_max`,
  `c=ᾱ/s̄` taken from that window's own carried α/s (exact at the reference point, see
  `build_tightened_actuator_constraint_convex`). Falls back to the bilinear form when no
  carried α/s is yet available (e.g. the very first window). Defaults to `false`
  (unchanged bilinear behavior for every existing caller) -- an opt-in numerics fix for
  when the bilinear form shows AL/conditioning distress despite comfortable physical
  headroom, not a physical relaxation.
- `trans_input_scale` : per-agent UNIT-CONVERSION factor for the translational channel.
  Both the tightening bound above and the runtime ancillary correction
  (`lsh_boundary_layer_control`) are written for a 2nd-order system with the control
  DIRECTLY equal to acceleration (`ẍ=u`, i.e. `B(x)=I`, mass=1). For a force-controlled
  agent (`ẍ=u/m`, i.e. `B(x)=I/m`), the physically correct forms are
  `‖F*‖ + m·α·s ≤ F_max` and `δF = -m·α·s·sat(·)` (the manuscript's `B(x)⁻¹` factor
  becomes `m·I` here). Defaults to `1.0` for every agent, exactly reproducing the
  original acceleration-control formulas byte-for-byte, so existing callers (e.g. an
  acceleration-controlled game) are completely unaffected. A force-controlled game
  passes `trans_input_scale=[m_1, m_2, ...]`.
- `tighten_rot_actuator`, `rot_actuator_max` : TORQUE tightening `‖τ*‖ + ‖J_rot‖₂·α_r·s_r
  ≤ rot_actuator_max` (**nonconvex**, ditto — session fix: previously missing the `J_rot`
  factor entirely, inconsistent with the runtime ancillary's own `J`-scaled correction, see
  `build_tightened_rot_actuator_constraint`'s docstring)
- `tighten_rot_rate`, `rot_rate_max` : RATE (state) tightening `‖ω*‖ + s_r ≤ rot_rate_max`
  (convex — see `build_tightened_rot_rate_constraint`). Was previously undocumented at this
  level; do not confuse with `tighten_rot_actuator` above (control vs. state, TORQUE vs.
  RATE — the stale prior wording here described `tighten_rot_actuator` using this rate
  formula, an internal documentation inconsistency, now fixed).
- `constraint_tightening`, `collision_min_distance`, `coupled_pairs` : collision
  margin `d_min − ‖r_i−r_j‖ + s_i + s_j ≤ 0` (nonconvex; Eq. 5 worked example)

# Ancillary (execution-time)
- `ancillary::Symbol` : `:lsh_boundary_layer` (translational + rotational
  boundary-layer sliding, LSH-faithful) or `:none`
- `λ` : sliding-surface slope (translational), per agent

# Co-optimized bandwidth (manuscript Eq. 10, `α̇=v`)
- `α_fixed`, `α_fixed_r`, `α_fixed_d` : despite the legacy names these are now the
  **INITIAL** bandwidths α₀ per channel, not design constants. The names are retained so
  existing `DTMPCParameters(...)` call sites keep working unchanged.
- `α_lo`, `α_max` : hard box on the α STATE, applied to every channel.
  **`α_lo > 0` is load-bearing**: the tube ODE is unstable for α ≤ 0 and the exact recursion
  divides by α. A prior implementation without this bound had the optimizer drive α negative
  and the width blow up.
- `v_max` : bound on the artificial control `|v| ≤ v_max` (ControlBounds).

# Numerical safety
- `s_max` : upper bound on any tube width (StateBounds; `s≥0` always enforced)
"""
struct DTMPCParameters{T}
    s₀::Vector{T}
    α_fixed::Vector{T}
    D::Vector{T}
    η::Vector{T}
    κ_speed::Vector{T}
    pos_dim::Int

    constraint_tightening::Bool
    collision_min_distance::T
    coupled_pairs::Vector{Tuple{Int, Int}}
    tighten_actuator::Bool
    actuator_max::Vector{T}
    trans_input_scale::Vector{T}
    tighten_actuator_convex::Bool

    ancillary::Symbol
    λ::Vector{T}
    s_max::T

    # Co-optimized bandwidth box (manuscript Problem 2). Scalars, applied to every
    # channel of every agent — matching the existing `s_max` precedent.
    α_lo::T
    α_max::T
    v_max::T

    rotational_tube::Vector{Bool}
    s₀_r::Vector{T}
    α_fixed_r::Vector{T}
    D_r::Vector{T}
    η_r::Vector{T}
    rot_dim::Int
    rot_state_offset::Int
    rot_ctrl_offset::Int
    tighten_rot_actuator::Bool
    rot_actuator_max::Vector{T}
    # SECOND-ORDER rotational channel (ω is a state, control is torque): λ_r is the
    # sliding-surface rate/attitude mixing coefficient (rotational analogue of the
    # translational `λ`); J_rot is each agent's diagonal principal inertia, needed by
    # the ancillary to convert a desired angular-acceleration correction into torque
    # (τ = J·ω̇_correction). tighten_rot_rate/rot_rate_max tighten the ω STATE (not a
    # control) against a physical max slew rate: ‖ω*‖+s_r ≤ rot_rate_max.
    λ_r::Vector{T}
    J_rot::Vector{Vector{T}}   # per-agent diagonal principal inertia [J1,J2,J3] (no
                               # StaticArrays dependency in this package — plain Vector)
    tighten_rot_rate::Bool
    rot_rate_max::Vector{T}

    detection_tube::Vector{Bool}
    s₀_d::Vector{T}
    α_fixed_d::Vector{T}
    D_d::Vector{T}
    η_d::Vector{T}
    det_state_offset::Int
    det_limit::T

    function DTMPCParameters{T}(
        n_agents::Int;
        s₀::Union{Nothing, Vector{T}} = nothing,
        α_fixed::Union{Nothing, Vector{T}} = nothing,
        D::Union{Nothing, Vector{T}} = nothing,
        η::Union{Nothing, Vector{T}} = nothing,
        κ_speed::Union{Nothing, Vector{T}} = nothing,
        pos_dim::Int = 2,
        constraint_tightening::Bool = true,
        collision_min_distance::T = T(1.0),
        coupled_pairs::Vector{Tuple{Int, Int}} = Tuple{Int, Int}[],
        tighten_actuator::Bool = false,
        actuator_max::Union{Nothing, Vector{T}} = nothing,
        trans_input_scale::Union{Nothing, Vector{T}} = nothing,
        tighten_actuator_convex::Bool = false,
        ancillary::Symbol = :none,
        λ::Union{Nothing, Vector{T}} = nothing,
        s_max::T = T(1e3),
        # Co-optimized bandwidth box (manuscript Problem 2). Defaults match the values the
        # previously-removed α-as-state implementation used (recorded in the project's
        # dtmpc_lsh_faithful notes), so behaviour is comparable to that known-working design.
        α_lo::T  = T(1e-2),
        α_max::T = T(5.0),
        v_max::T = T(1.0),
        rotational_tube::Union{Nothing, Vector{Bool}} = nothing,
        s₀_r::Union{Nothing, Vector{T}} = nothing,
        α_fixed_r::Union{Nothing, Vector{T}} = nothing,
        D_r::Union{Nothing, Vector{T}} = nothing,
        η_r::Union{Nothing, Vector{T}} = nothing,
        rot_dim::Int = 3,
        rot_state_offset::Int = 6,
        rot_ctrl_offset::Int = 3,
        tighten_rot_actuator::Bool = false,
        rot_actuator_max::Union{Nothing, Vector{T}} = nothing,
        λ_r::Union{Nothing, Vector{T}} = nothing,
        J_rot::Union{Nothing, Vector{Vector{T}}} = nothing,
        tighten_rot_rate::Bool = false,
        rot_rate_max::Union{Nothing, Vector{T}} = nothing,
        detection_tube::Union{Nothing, Vector{Bool}} = nothing,
        s₀_d::Union{Nothing, Vector{T}} = nothing,
        α_fixed_d::Union{Nothing, Vector{T}} = nothing,
        D_d::Union{Nothing, Vector{T}} = nothing,
        η_d::Union{Nothing, Vector{T}} = nothing,
        det_state_offset::Int = 9,
        det_limit::T = T(log(2)),
    ) where {T}
        @assert n_agents > 0 "n_agents must be positive"
        @assert ancillary in (:none, :lsh_boundary_layer) "unknown ancillary $ancillary"

        s₀      = s₀      !== nothing ? s₀      : fill(T(0.1),  n_agents)
        α_fixed = α_fixed !== nothing ? α_fixed : fill(T(1.0),  n_agents)
        D       = D       !== nothing ? D       : fill(T(0.1),  n_agents)
        η       = η       !== nothing ? η       : fill(T(0.01), n_agents)
        κ_speed = κ_speed !== nothing ? κ_speed : fill(T(0.0),  n_agents)
        λ       = λ       !== nothing ? λ       : fill(T(0.1),  n_agents)
        actuator_max = actuator_max !== nothing ? actuator_max : fill(T(Inf), n_agents)
        @assert !tighten_actuator || all(isfinite, actuator_max) "tighten_actuator requires finite actuator_max"
        trans_input_scale = trans_input_scale !== nothing ? trans_input_scale : fill(T(1.0), n_agents)
        @assert all(trans_input_scale .> 0) "trans_input_scale must be > 0"
        @assert all(α_fixed .> 0) "α_fixed must be > 0 (tube ODE stable iff α>0)"

        rotational_tube = rotational_tube !== nothing ? rotational_tube : fill(false, n_agents)
        s₀_r      = s₀_r      !== nothing ? s₀_r      : fill(T(0.05), n_agents)
        α_fixed_r = α_fixed_r !== nothing ? α_fixed_r : fill(T(0.5),  n_agents)
        D_r       = D_r       !== nothing ? D_r       : fill(T(1e-3),n_agents)
        η_r       = η_r       !== nothing ? η_r       : fill(T(1e-4),n_agents)
        rot_actuator_max = rot_actuator_max !== nothing ? rot_actuator_max : fill(T(Inf), n_agents)
        @assert !tighten_rot_actuator || all(isfinite, rot_actuator_max) "tighten_rot_actuator requires finite rot_actuator_max"
        λ_r   = λ_r   !== nothing ? λ_r   : fill(T(0.1), n_agents)
        J_rot = J_rot !== nothing ? J_rot : [fill(T(1.0), 3) for _ in 1:n_agents]
        rot_rate_max = rot_rate_max !== nothing ? rot_rate_max : fill(T(Inf), n_agents)
        @assert !tighten_rot_rate || all(isfinite, rot_rate_max) "tighten_rot_rate requires finite rot_rate_max"
        @assert length(J_rot) == n_agents && all(length(j) == 3 for j in J_rot) "J_rot must be n_agents vectors of length 3"

        detection_tube = detection_tube !== nothing ? detection_tube : fill(false, n_agents)
        s₀_d      = s₀_d      !== nothing ? s₀_d      : fill(T(0.02), n_agents)
        α_fixed_d = α_fixed_d !== nothing ? α_fixed_d : fill(T(0.5),  n_agents)
        D_d       = D_d       !== nothing ? D_d       : fill(T(1e-3),n_agents)
        η_d       = η_d       !== nothing ? η_d       : fill(T(1e-4),n_agents)

        for v in (s₀, α_fixed, D, η, κ_speed, λ, actuator_max, trans_input_scale,
                  s₀_r, α_fixed_r, D_r, η_r, rot_actuator_max, λ_r, rot_rate_max,
                  s₀_d, α_fixed_d, D_d, η_d)
            @assert length(v) == n_agents "all per-agent parameter vectors must have length n_agents"
        end
        @assert length(rotational_tube) == n_agents && length(detection_tube) == n_agents
        @assert all(s₀ .≥ 0) "initial tube widths must be non-negative"
        @assert pos_dim ≥ 1 "pos_dim must be ≥ 1"
        any(rotational_tube) && @assert(all(α_fixed_r[rotational_tube] .> 0), "α_fixed_r must be > 0 where rotational_tube")
        any(detection_tube)  && @assert(all(α_fixed_d[detection_tube]  .> 0), "α_fixed_d must be > 0 where detection_tube")
        # Co-optimized-bandwidth box. α_lo>0 is not a nicety: α≤0 makes the tube ODE unstable
        # and the exact width recursion divides by α.
        @assert α_lo > 0 "α_lo must be > 0 (tube ODE unstable and δ̂ divides by α for α≤0)"
        @assert α_max > α_lo "α_max must exceed α_lo"
        @assert v_max > 0 "v_max must be > 0"
        # The INITIAL α of every enabled channel must lie inside the box, else the window's
        # very first state is infeasible by construction.
        @assert all(α_lo .≤ α_fixed .≤ α_max) "α_fixed (initial α) must lie in [α_lo, α_max]"
        any(rotational_tube) && @assert(all(α_lo .≤ α_fixed_r[rotational_tube] .≤ α_max),
                                        "α_fixed_r (initial α) must lie in [α_lo, α_max]")
        any(detection_tube)  && @assert(all(α_lo .≤ α_fixed_d[detection_tube] .≤ α_max),
                                        "α_fixed_d (initial α) must lie in [α_lo, α_max]")

        for (i, j) in coupled_pairs
            @assert 1 ≤ i ≤ n_agents && 1 ≤ j ≤ n_agents && i ≠ j "invalid coupled pair"
        end

        new{T}(s₀, α_fixed, D, η, κ_speed, pos_dim,
               constraint_tightening, collision_min_distance, coupled_pairs,
               tighten_actuator, actuator_max, trans_input_scale, tighten_actuator_convex,
               ancillary, λ, s_max,
               α_lo, α_max, v_max,
               rotational_tube, s₀_r, α_fixed_r, D_r, η_r,
               rot_dim, rot_state_offset, rot_ctrl_offset,
               tighten_rot_actuator, rot_actuator_max,
               λ_r, J_rot, tighten_rot_rate, rot_rate_max,
               detection_tube, s₀_d, α_fixed_d, D_d, η_d,
               det_state_offset, det_limit)
    end
end

DTMPCParameters(::Type{T}, n_agents::Int; kwargs...) where {T} =
    DTMPCParameters{T}(n_agents; kwargs...)

n_agents(p::DTMPCParameters) = length(p.s₀)

"""
    tube_channels(params, i) -> Vector{TubeChannel}

Ordered tube channels for agent `i`: always `[translational, (rotational),
(detection)]`, matching augmented-state layout order.
"""
function tube_channels(p::DTMPCParameters{T}, i::Int) where {T}
    chans = TubeChannel{T}[TubeChannel{T}(:translational, p.s₀[i], p.α_fixed[i], p.D[i], p.η[i],
                                          p.κ_speed[i], 0, 2p.pos_dim, 0, p.pos_dim)]
    if p.rotational_tube[i]
        # SECOND-ORDER rotational channel: state_dim=2·rot_dim covers the joint [p;ω]
        # block (attitude + body rate, ω now a STATE not a control); ctrl_dim=rot_dim
        # is unchanged (torque, same dimension as the old ω control).
        push!(chans, TubeChannel{T}(:rotational, p.s₀_r[i], p.α_fixed_r[i], p.D_r[i], p.η_r[i],
                                    zero(T), p.rot_state_offset, 2*p.rot_dim, p.rot_ctrl_offset, p.rot_dim))
    end
    if p.detection_tube[i]
        push!(chans, TubeChannel{T}(:detection, p.s₀_d[i], p.α_fixed_d[i], p.D_d[i], p.η_d[i],
                                    zero(T), p.det_state_offset, 1, 0, 0))
    end
    return chans
end

n_tube_channels_agent(p::DTMPCParameters, i::Int) = 1 + Int(p.rotational_tube[i]) + Int(p.detection_tube[i])
uses_multi_tube(p::DTMPCParameters) = any(p.rotational_tube) || any(p.detection_tube)

# ============================================================================
# Augmented layout:
#     x̃ᵢ = [xᵢ(nᵢ); s₁…s_nch; α₁…α_nch]      ũᵢ = [uᵢ(mᵢ); v₁…v_nch]
# α is appended AFTER all widths so the physical block stays leading and every
# `s` index is unchanged from the pre-α layout (see the module header).
# ============================================================================

"""
    augmented_state_dims(state_dims, nt_vec) -> Vector{Int}

Per-agent augmented state dims `nᵢ + 2·nch_i`: physical state, one scalar WIDTH `s` per tube
channel, then one scalar BANDWIDTH `α` per channel (co-optimized, manuscript Eq. 10).
"""
augmented_state_dims(state_dims::Vector{Int}, nt_vec::Vector{Int}) =
    [state_dims[i] + 2 * nt_vec[i] for i in eachindex(state_dims)]
augmented_state_dims(state_dims::Vector{Int}, n_tubes::Int) =
    [n + 2 * n_tubes for n in state_dims]

"""
    augmented_control_dims(control_dims, nt_vec) -> Vector{Int}

Per-agent augmented control dims `mᵢ + nch_i`: physical control plus one artificial
bandwidth-rate `v` per tube channel (`α̇=v`, manuscript Eq. 10).
"""
augmented_control_dims(control_dims::Vector{Int}, nt_vec::Vector{Int}) =
    [control_dims[i] + nt_vec[i] for i in eachindex(control_dims)]

block_offsets(dims::Vector{Int}) = [0; cumsum(dims)[1:end-1]]

"""
    agent_physical_block(x_joint, state_dims, i) -> xᵢ

Slice agent `i`'s physical state block from the joint PHYSICAL state vector.
"""
function agent_physical_block(x_joint::AbstractVector, state_dims::Vector{Int}, i::Int)
    off = block_offsets(state_dims)
    return x_joint[off[i] + 1 : off[i] + state_dims[i]]
end

_as_scalar(x::Number) = x
_as_scalar(x::AbstractVector) = x[1]

"""
    extract_agent_augmented(joint_aug, aug_offsets, state_dims, nt_vec, i) -> (xᵢ, s, α)

Extract agent `i`'s physical state, per-channel widths `s`, and per-channel bandwidths `α`
from a JOINT augmented state vector. (Returns a 3-tuple: `α` was added when the bandwidth
became a co-optimized state.)
"""
function extract_agent_augmented(joint_aug::AbstractVector, aug_offsets::Vector{Int},
                                 state_dims::Vector{Int}, nt_vec::Vector{Int}, i::Int)
    o = aug_offsets[i]; n = state_dims[i]; nch = nt_vec[i]
    x_i = joint_aug[o + 1 : o + n]
    s = [joint_aug[o + n + k] for k in 1:nch]
    α = [joint_aug[o + n + nch + k] for k in 1:nch]
    return (x_i, s, α)
end

# z_full layout translation (augmented → physical). BOTH strides now differ from the physical
# ones: the state carries 2·nch extra scalars (s block then α block) and the control carries
# nch extra (the bandwidth rates v).
#
# This function feeds the LIVE CROSS-PLAYER COUPLING (`z_full`) that both games read opponent
# positions/attitudes out of. An arithmetic slip here does not error — it silently hands the
# stage costs garbage, which then shows up only as a mysteriously wrong equilibrium. It is
# covered directly by a round-trip test rather than left to downstream symptoms.
function _physical_zfull(aug_zfull, sdims::Vector{Int}, cdims::Vector{Int}, N::Int,
                         nt_vec::Vector{Int})
    np   = length(sdims)
    Tz   = eltype(aug_zfull[1])
    phys = Vector{Vector{Tz}}(undef, np)
    for j in 1:np
        n = sdims[j]; m = cdims[j]; nch = nt_vec[j]
        n_aug = n + 2 * nch          # physical state + s block + α block
        m_aug = m + nch              # physical control + v block
        astr = n_aug + m_aug
        pstr = n + m
        zj   = aug_zfull[j]
        out  = Vector{Tz}(undef, pstr * N + n)
        @inbounds for k in 0:(N - 1)
            ab = k * astr; pb = k * pstr
            out[pb+1   : pb+n]   = @view zj[ab+1       : ab+n]
            out[pb+n+1 : pb+n+m] = @view zj[ab+n_aug+1 : ab+n_aug+m]
        end
        ab = N * astr; pb = N * pstr
        @inbounds out[pb+1 : pb+n] = @view zj[ab+1 : ab+n]
        phys[j] = out
    end
    return phys
end

@inline function _with_physical_zfull(p, sdims::Vector{Int}, cdims::Vector{Int}, N::Int,
                                      nt_vec::Vector{Int})
    (p isa NamedTuple && hasproperty(p, :z_full) && p.z_full !== nothing) || return p
    return merge(p, (; z_full = _physical_zfull(p.z_full, sdims, cdims, N, nt_vec)))
end

# ============================================================================
# Augmented dynamics (discrete): physical map unchanged; each channel's width
# advances by the EXACT tight recursion s_{k+1}=ρ_k·s_k+δ̂_k(x̄,ū) (Remark 1.1:
# the equality is the tight/dominating instance of the theory's inequality
# (6) — NO SCP linearization error attaches to a constant/affine δ̂; a
# κ_speed≠0 majorant is smooth-nonlinear and is re-linearized every SCP
# iterate like any other own-agent nonlinear quantity, obligation V-3a).
# ============================================================================

"""
    build_augmented_player_dynamics_discrete(orig_f, n_i, m_i, chans, x_other_frozen,
        dt, sdims, cdims, N, nt_vec; Δ_funcs) -> Function

Discrete next-state augmented field: `x_next = orig_f(x_i,u_i,p,t)` (physical,
z_full-translated for cross-player coupling), plus one width recursion per
channel in `chans` (ordered `[translational, rotational?, detection?]`).
`Δ_funcs[k](x_i,x_other) -> scalar` supplies each channel's own-state model-error
bound (LSH's Δ; §S separable); defaults to `nothing` per channel ⇒ Δ≡0 (pure
`D+η` bound). The translational channel's Δ is normally the SAME per-agent
uncertainty model passed to `DynamicTubeMPCProblem` (see `build_augmented_game`).
"""
function build_augmented_player_dynamics_discrete(
    orig_f::Function, n_i::Int, m_i::Int, chans::Vector{TubeChannel{T}},
    x_other_frozen::AbstractVector,
    dt::Real, sdims::Vector{Int}, cdims::Vector{Int}, N::Int, nt_vec::Vector{Int};
    Δ_funcs::Vector{<:Union{Nothing,Function}} = fill(nothing, length(chans)),
    α_lo::Real = 1e-2,
) where {T}
    nch = length(chans)
    return function f_aug(x_aug, u_aug, p, t)
        x_i = x_aug[1:n_i]
        u_i = u_aug[1:m_i]
        p_phys = _with_physical_zfull(p, sdims, cdims, N, nt_vec)
        x_next = orig_f(x_i, u_i, p_phys, t)
        Tz = promote_type(eltype(x_aug), eltype(u_aug))
        s_next = Vector{Tz}(undef, nch)
        α_next = Vector{Tz}(undef, nch)
        @inbounds for k in 1:nch
            c   = chans[k]
            s_k = x_aug[n_i + k]              # width block
            α_k = x_aug[n_i + nch + k]        # bandwidth block (appended after all widths)
            v_k = u_aug[m_i + k]              # artificial bandwidth-rate control
            Δ_i = Δ_funcs[k] === nothing ? zero(T) : _as_scalar(Δ_funcs[k](x_i, x_other_frozen))
            αg  = _α_guard(α_k, α_lo)
            if c.kind === :translational && false
                # Sampled-data recursion (session correction): the componentwise sliding-
                # surface ancillary this channel governs is a ZOH sampled controller, not a
                # continuously-acting one -- see `sampled_tube_next`'s docstring.
                #
                # TEMPORARILY DISABLED (`&& false`, session diagnostic -- not a revert of the
                # underlying work): this in-progress, still-uncommitted change to the
                # translational tube recursion is what makes the certified pitch15_wolf3pct
                # Stage-5 solar scenario (tuned against the OLD continuous-time recursion)
                # go infeasible at window 1. Bypassed here only long enough to reproduce that
                # certified result and pull real Monte-Carlo/FALCON-timing numbers from it;
                # the sampled-data recursion itself is untouched below and should be
                # re-enabled (this `&& false` removed) once pitch15_wolf3pct -- or whichever
                # scenario needs it -- is re-tuned/re-certified against it.
                β_i = _tube_beta(c, x_i, u_i, Δ_i)
                s_next[k] = sampled_tube_next(s_k, α_k, β_i, dt)
            else
                # Rotational/detection channels UNCHANGED: exact one-step integral of
                # ṡ=−αs+β, now with the per-node co-optimized α.
                s_next[k] = ρ_of(αg, dt) * s_k + tube_delta_hat(c, x_i, u_i, Δ_i, dt, α_k, α_lo)
            end
            # α̇ = v (manuscript Eq. 10), exact under ZOH v.
            α_next[k] = α_k + dt * v_k
        end
        return vcat(x_next, s_next, α_next)
    end
end

# ============================================================================
# Augmented cost: THE ORIGINAL cost, unchanged (Thm 1 — tube widths carry NO
# objective curvature; they are regularized only by the inner solver's own
# proximal term). Only the physical (x,u) slice + z_full translation.
# ============================================================================

"""
    build_augmented_cost(orig_obj, n_i, m_i; sdims, cdims, N, nt_vec) -> PlayerObjective

Wrap a player's NonlinearStageCost/TerminalCost to slice the physical `(x,u)`
out of the augmented `[x; s₁;…;s_nch]` state and translate `p.z_full` to
physical layout. Adds NOTHING to the cost (Theorem 1).
"""
function build_augmented_cost(
    orig_obj::PlayerObjective, n_i::Int, m_i::Int;
    sdims::Vector{Int}, cdims::Vector{Int}, N::Int, nt_vec::Vector{Int},
)
    sc = orig_obj.stage_cost
    @assert sc isa NonlinearStageCost "build_augmented_cost requires a NonlinearStageCost (FALCON path)"
    orig_sc = orig_obj.stage_cost.func
    orig_tc = orig_obj.terminal_cost.func

    aug_stage = function (x_aug, u_aug, p, t)
        x = @view x_aug[1:n_i]; u = @view u_aug[1:m_i]
        p2 = _with_physical_zfull(p, sdims, cdims, N, nt_vec)
        return orig_sc(x, u, p2, t)
    end
    aug_term = function (x_aug, p)
        x = @view x_aug[1:n_i]
        p2 = _with_physical_zfull(p, sdims, cdims, N, nt_vec)
        return orig_tc(x, p2)
    end
    return PlayerObjective(orig_obj.player_id,
                           NonlinearStageCost(aug_stage),
                           NonlinearTerminalCost(aug_term),
                           orig_obj.scaling)
end

# ============================================================================
# Carry the game's own shared constraints through augmentation. Without this,
# `build_augmented_game` silently DROPS every user-supplied shared constraint —
# `shared = Any[]` below used to be repopulated ONLY with the tightened collision
# constraint, so any other coupling constraint on the original game (e.g. a coupled
# detection/visibility cap) would exist in an open-loop FALCON solve but vanish the
# moment the same game was wrapped in DTMPC. Private constraints and objectives are
# unaffected — they are already carried through via `priv_i`/`build_augmented_cost`.
# ============================================================================

"""
    build_augmented_shared_constraint(sc, aug_offsets, sdims, cdims, nt_vec) -> shared constraint

Wrap a user-supplied shared constraint (any subtype of `AbstractSharedInequality` /
`AbstractSharedEquality` — `ProximityConstraint`, `SharedNonlinearInequality`, …) so it
evaluates against the PHYSICAL joint state and control sliced out of the augmented ones.
Each agent's augmented state block is `[physical(sdims[i]); s(nch); α(nch)]` and its
augmented control block is `[physical(cdims[i]); v(nch)]`, so BOTH arguments need
re-slicing (the control did not, before the bandwidth became co-optimized — forgetting it
would silently feed a user constraint the `v` values as if they were physical controls).

The per-agent tube widths and bandwidths are additionally passed through `p` as
`tube_widths` / `tube_alphas` (vectors of per-channel vectors). A constraint that wants to
tighten itself against its own tube can read them defensively with
`hasproperty(p, :tube_widths)`; constraints that ignore them are unaffected. Passing them
via `p` rather than exposing them as decision variables is exact here, not a shortcut:
`s` and `α` evolve autonomously within a window (their recursion reads only its own
previous value, the nominal, and `v`), so from the constraint's viewpoint they are
per-node constants.

Dispatches through `evaluate_constraint` rather than a `.func` field so it works for
constraint types that compute their value programmatically (`ProximityConstraint`) as
well as ones that store a `func` closure. The re-slice `x_aug ↦ x_phys` is a LINEAR
projection, so wrapping with `jacobian=nothing` (ForwardDiff over the wrapper) is exact,
not an approximation — it just forgoes any analytic jacobian the original constraint
type had (e.g. `ProximityConstraint`'s), which only costs performance, not correctness.
"""
function build_augmented_shared_constraint(sc, aug_offsets::Vector{Int}, sdims::Vector{Int},
                                           cdims::Vector{Int}, nt_vec::Vector{Int})
    np = length(sdims)
    cu_offsets = block_offsets(augmented_control_dims(cdims, nt_vec))
    aug_func = function (x_aug, u_aug, p, t)
        x_phys = vcat([x_aug[aug_offsets[j]+1 : aug_offsets[j]+sdims[j]] for j in 1:np]...)
        u_phys = vcat([u_aug[cu_offsets[j]+1 : cu_offsets[j]+cdims[j]] for j in 1:np]...)
        widths = [[x_aug[aug_offsets[j]+sdims[j]+k]            for k in 1:nt_vec[j]] for j in 1:np]
        alphas = [[x_aug[aug_offsets[j]+sdims[j]+nt_vec[j]+k]  for k in 1:nt_vec[j]] for j in 1:np]
        p2 = p isa NamedTuple ? merge(p, (; tube_widths = widths, tube_alphas = alphas)) : p
        return evaluate_constraint(sc, x_phys, u_phys, p2, t)
    end
    players = get_players(sc)
    dim = constraint_output_dim(sc)
    return sc isa AbstractSharedEquality ? SharedEquality(players; func=aug_func, dim=dim) :
                                            SharedInequality(players; func=aug_func, dim=dim)
end

# ============================================================================
# Decision-dependent tightened constraints (§3, Eq. 4/5)
# ============================================================================

"""
    build_tightened_collision_constraint(i, j, aug_offsets, n_i, n_j, pos_dim, d_min; ε)
        -> SharedNonlinearInequality

`c_ij(x̃) = d_min − ‖rᵢ−rⱼ‖ + sᵢ + sⱼ ≤ 0` (Eq. 5 worked example). Reads each
agent's TRANSLATIONAL channel width (augmented index `nᵢ+1`). Nonconvex
(`−‖·‖` is concave) — AL-penalized, as any collision-avoidance constraint.
"""
function build_tightened_collision_constraint(
    i::Int, j::Int, aug_offsets::Vector{Int}, n_i::Int, n_j::Int, pos_dim::Int, d_min;
    ε = 1e-6,
)
    oi = aug_offsets[i]; oj = aug_offsets[j]
    si_idx = oi + n_i + 1; sj_idx = oj + n_j + 1
    func = function (x, u, p, t)
        r_i = x[oi + 1 : oi + pos_dim]; r_j = x[oj + 1 : oj + pos_dim]
        dist = sqrt(sum(abs2, r_i .- r_j) + ε^2)
        return [d_min - dist + x[si_idx] + x[sj_idx]]
    end
    return SharedNonlinearInequality([i, j], func, nothing, 1)
end

"""
    build_tightened_actuator_constraint(i, n_i, nch, chan_idx, pos_dim, a_max; ε)
        -> PrivateInequality

Corollary-2 actuator tightening `‖a*‖ + α·s − a_max ≤ 0`, reading BOTH the width `s` and
the bandwidth `α` from the augmented state.

**NONCONVEX**, and deliberately so: with α co-optimized (manuscript Eq. 10) the reserve
`α·s` is BILINEAR in two decision variables, so this can no longer be marked
`is_convex=true` and moves onto FALCON's AL-penalized path instead of its hard,
always-linearized, no-slack path. `docs/theory.md` §11 identified exactly this as the cost
of co-optimizing α; it is accepted, not overlooked. If AL convergence degrades where
`tighten_actuator=true`, the intended fallback is the convex conservative surrogate
`α·s ≤ (α²/c + c·s²)/2` (exact at `c = α/s`), NOT reverting to a fixed α.
"""
function build_tightened_actuator_constraint(i::Int, n_i::Int, nch::Int, chan_idx::Int,
                                             pos_dim::Int, a_max, scale = 1; ε = 1e-9)
    s_idx = n_i + chan_idx
    α_idx = n_i + nch + chan_idx
    func = function (x, u, p, t)
        T = promote_type(eltype(x), eltype(u))
        a2 = zero(T)
        @inbounds for k in 1:pos_dim
            a2 += u[k]^2
        end
        anorm = sqrt(a2 + T(ε)^2)
        # `scale` converts the co-optimized bandwidth/width product (an ACCELERATION-
        # authority reserve, `α·s`) into the control's own units -- 1.0 for an
        # acceleration-controlled agent (B(x)=I), or the agent's mass for a force-
        # controlled agent (B(x)=I/m, so B(x)⁻¹=m·I) -- see `trans_input_scale` docstring.
        return [anorm + T(scale) * x[α_idx] * x[s_idx] - T(a_max)]
    end
    return PrivateInequality(i; func = func, dim = 1, is_convex = false)
end

"""
    build_tightened_actuator_constraint_convex(i, n_i, nch, chan_idx, pos_dim, a_max, scale, c;
        ε) -> PrivateInequality

CONVEX conservative surrogate for the actuator tightening above, per that function's own
docstring: `α·s ≤ (α²/c + c·s²)/2` for any fixed `c>0`, by AM-GM (equality iff `α/s=c`). The
caller supplies `c` as the CURRENT/reference ratio `ᾱ/s̄` (e.g. this window's carried α/s),
making the surrogate exact at that point while remaining a valid (conservative) bound
elsewhere and, crucially, CONVEX in `(α,s)` for fixed `c` — unlike the bilinear `α·s` form,
this can be marked `is_convex=true` and takes FALCON's hard, always-linearized, no-AL-
penalty path. Intended for use when `tighten_actuator=true`'s bilinear form shows AL/
conditioning distress (large `ΔJ_ζ`, degraded convergence) despite comfortable physical
headroom (`r_total=(‖F*‖+scale·α·s)/a_max` well below 1) -- a numerics fix, not a physical
relaxation: at `c=ᾱ/s̄` the two forms agree exactly, so this does not change what the
constraint MEANS at the reference point, only how tightly FALCON can enforce it away from
it.
"""
function build_tightened_actuator_constraint_convex(i::Int, n_i::Int, nch::Int, chan_idx::Int,
                                                     pos_dim::Int, a_max, scale, c; ε = 1e-9)
    s_idx = n_i + chan_idx
    α_idx = n_i + nch + chan_idx
    func = function (x, u, p, t)
        T = promote_type(eltype(x), eltype(u))
        a2 = zero(T)
        @inbounds for k in 1:pos_dim
            a2 += u[k]^2
        end
        anorm = sqrt(a2 + T(ε)^2)
        surrogate = (x[α_idx]^2 / T(c) + T(c) * x[s_idx]^2) / 2
        return [anorm + T(scale) * surrogate - T(a_max)]
    end
    return PrivateInequality(i; func = func, dim = 1, is_convex = true)
end

"""
    build_tightened_rot_actuator_constraint(i, n_i, nch, chan_idx, rot_ctrl_offset, rot_dim,
        τ_max, Jdiag; ε) -> PrivateInequality

Rotational actuator (TORQUE) counterpart: `‖τ*‖ + ‖J‖₂·α_r·s_r − τ_max ≤ 0`, reading both
the rotational width and bandwidth from the augmented state. **NONCONVEX** for the same
bilinearity reason as [`build_tightened_actuator_constraint`](@ref) — see its docstring.

`Jdiag` is REQUIRED (session correctness fix — previously omitted entirely, an internal
inconsistency with the runtime ancillary `lsh_boundary_layer_control`'s rotational
correction `δτ = -J·α_r·Φ_r·sat(s̃/Φ_r)`, which explicitly converts the co-optimized
angular-ACCELERATION reserve `α_r·Φ_r` into torque via `J`; this tightening bound must use
the SAME conversion to stay consistent with what the ancillary can actually apply at
runtime). The physical actuator constraint is Euclidean, `‖τ*+δτ‖₂≤τ_max`, and each
per-axis saturation component has magnitude ≤1, so `‖δτ‖₂ ≤ ‖J‖₂·α_r·s_r` where
`‖J‖₂=sqrt(J₁²+J₂²+J₃²)` for the diagonal inertia vector — a clean, sufficient (if slightly
conservative for non-uniform `J`) tightening. `Jgain=1.0` recovers the old (dimensionally
wrong) behavior exactly, so this is opt-in only insofar as a caller passes `Jdiag=[1,1,1]`;
every real caller should pass its actual diagonal inertia.
"""
function build_tightened_rot_actuator_constraint(
    i::Int, n_i::Int, nch::Int, chan_idx::Int, rot_ctrl_offset::Int, rot_dim::Int,
    τ_max, Jdiag; ε = 1e-9,
)
    s_idx = n_i + chan_idx
    α_idx = n_i + nch + chan_idx
    Jgain = norm(Jdiag)
    func = function (x, u, p, t)
        T = promote_type(eltype(x), eltype(u))
        w2 = zero(T)
        @inbounds for k in 1:rot_dim
            w2 += u[rot_ctrl_offset + k]^2
        end
        wnorm = sqrt(w2 + T(ε)^2)
        return [wnorm + T(Jgain) * x[α_idx] * x[s_idx] - T(τ_max)]
    end
    return PrivateInequality(i; func = func, dim = 1, is_convex = false)
end

"""
    build_tightened_rot_rate_constraint(i, n_i, chan_idx, rot_state_offset, rot_dim,
        rate_max; ε) -> PrivateInequality

STATE tightening for the body rate `ω` (now a state, not a control, under the
second-order rotational channel): `‖ω*‖ + s_r − rate_max ≤ 0` (Eq. 4 form, `P=I`
over the channel's full tube width — the worst case where the entire tube
budget sits in the ω sub-block). Convex, `is_convex=true`. A physical max-slew-
rate bound (e.g. a wheel/actuator design limit), independent of the torque
tightening above.
"""
function build_tightened_rot_rate_constraint(
    i::Int, n_i::Int, chan_idx::Int, rot_state_offset::Int, rot_dim::Int,
    rate_max; ε = 1e-9,
)
    s_idx = n_i + chan_idx
    func = function (x, u, p, t)
        T = eltype(x)
        w2 = zero(T)
        @inbounds for k in 1:rot_dim
            w2 += x[rot_state_offset + rot_dim + k]^2
        end
        wnorm = sqrt(w2 + T(ε)^2)
        return [wnorm + x[s_idx] - T(rate_max)]
    end
    return PrivateInequality(i; func = func, dim = 1, is_convex = true)
end

"""
    build_tightened_detection_constraint(i, det_state_idx, s_d_idx, limit) -> PrivateInequality

`s_det + s_d − limit ≤ 0` (Eq. 4, `P=1`). Linear ⇒ `is_convex=true`.
"""
function build_tightened_detection_constraint(i::Int, det_state_idx::Int, s_d_idx::Int, limit)
    func = function (x, u, p, t)
        T = eltype(x)
        return [x[det_state_idx] + x[s_d_idx] - T(limit)]
    end
    return PrivateInequality(i; func = func, dim = 1, is_convex = true)
end

"""
    build_initial_containment_constraint(i, n_i, chan_idx, state_offset, state_dim; ε)
        -> PrivateInequality

Initial-containment constraint (§4): `‖x_measured[offset:offset+dim] −
x̄_i[offset:offset+dim]‖ ≤ s_i,k`, active ONLY at the window's first node
(`p.timestep==0`; trivially satisfied, `≤ -1`, elsewhere). Genuinely convex
(norm ball) — `is_convex=true`. `x_measured` is the NUMERIC measured value at
this channel's state slice, baked into the closure.

With a solver that frees the initial-state pin for the corresponding
components (e.g. FALCON's `free_initial`, via [`dtmpc_free_initial_mask`](@ref)),
this is the ONLY constraint tying the window's initial nominal/tube-width
pair to reality — realizing the theory's optimize-over-initial-nominal-state
convention. With a solver that keeps the classical hard pin (e.g. ALGAMES),
`x̄_i,{t|t}` IS the measurement, so this constraint reduces to `0 ≤ s_i,{t|t}`
and is trivially satisfied (a valid, more conservative special case).
"""
function build_initial_containment_constraint(
    i::Int, n_i::Int, chan_idx::Int, state_offset::Int, state_dim::Int,
    x_measured_slice::AbstractVector; ε = 1e-9,
)
    s_idx = n_i + chan_idx
    x_meas = collect(x_measured_slice)
    func = function (x, u, p, t)
        T = eltype(x)
        k = hasproperty(p, :timestep) ? p.timestep : -1
        k != 0 && return T[-one(T)]
        e = @view(x[state_offset+1 : state_offset+state_dim]) .- T.(x_meas)
        enorm = sqrt(sum(abs2, e) + T(ε)^2)
        return [enorm - x[s_idx]]
    end
    return PrivateInequality(i; func = func, dim = 1, is_convex = true)
end

"""
    build_initial_containment_constraint_sliding(i, n_i, chan_idx, pos_dim, lambda,
        x_measured_slice; ε) -> PrivateInequality

Session correction, TRANSLATIONAL channel only: the componentwise boundary-layer ancillary
(`lsh_boundary_layer_control`) regulates the sliding surface `s_j = v_j + lambda*r_j`
AXIS-BY-AXIS against one scalar `Phi` (`|s_j| <= Phi` per axis, independently for
j=1:pos_dim) -- an l_infty-style box on the sliding vector, NOT a Euclidean ball on the
concatenated `[r;v]` state slice (`build_initial_containment_constraint` above, still used
unchanged for rotational/detection channels). Enforcing
`||[r;v]-[r*;v*]||_2 <= Phi` is a DIFFERENT, dimensionally mixed (meters and m/s summed
inside one norm) certificate that neither matches what the ancillary defends nor what the
tube ODE's `beta` is derived against. This constraint instead ties the window's initial
nominal/tube-width pair to reality using the SAME sliding-surface quantity the ancillary
and the sampled-data recursion (`sampled_tube_next`) both use:
`max_j |v_j - v*_j + lambda*(r_j - r*_j)| <= Phi`.

`pos_dim` position/velocity axes are assumed (channel state = `[r(pos_dim); v(pos_dim)]`,
`state_offset`/`state_dim=2*pos_dim` as usual). Uses a smooth max (log-sum-exp is avoided in
favor of a plain smoothed max via `sqrt(x^2+eps^2)`-based softabs then a smooth max
surrogate) is NOT used here -- `maximum(abs, .)` is left non-smooth-relaxed deliberately
matched to a `PrivateNonlinearInequality`-compatible max-of-abs via `dim=pos_dim` (one row
per axis) rather than a single smoothed scalar, avoiding a smoothing bias in the constraint
that the ancillary itself does not have (its own per-axis saturation is likewise exact,
not smoothed, at the boundary `|s_j|=Phi`).

Session correction: `||v_tilde+lambda*r_tilde|| <= Phi` ALONE does not bound position
independently -- `v_tilde=-lambda*r_tilde` gives `s=0` for ARBITRARILY LARGE `r_tilde`, so
the free initial nominal position can cancel against velocity on this manifold while the
sliding constraint is perfectly satisfied (confirmed empirically: lowering
`lambda_trans` 1.0->0.02 collapsed Stage-4 to near-total passivity, NOT because a smaller
lambda is wrong for the ancillary -- the exact discrete counterfactual confirmed
`lambda_trans=0.02` DOES make the ancillary reduce sliding error -- but because it exposed
that this constraint alone never bounded position). `pos_pin_tol` (default `Inf`, meaning
"add no position constraint", the prior behavior): when finite, additionally pins the
window's free initial nominal POSITION to the measurement within `pos_pin_tol` (a
temporary diagnostic device -- pass `0.0`-ish for an near-exact pin, e.g. `1e-6`; the
permanent fix is a carried position-radius `R_r` via `R_{k+1}=exp(-lambda*dt)*R_k +
(1-exp(-lambda*dt))/lambda*Phi_k`, not implemented here yet).
"""
function build_initial_containment_constraint_sliding(
    i::Int, n_i::Int, chan_idx::Int, pos_dim::Int, lambda::Real,
    x_measured_slice::AbstractVector; pos_pin_tol::Real = Inf,
)
    s_idx = n_i + chan_idx
    x_meas = collect(x_measured_slice)
    has_pin = isfinite(pos_pin_tol)
    ndim = has_pin ? 4pos_dim : 2pos_dim
    func = function (x, u, p, t)
        T = eltype(x)
        k = hasproperty(p, :timestep) ? p.timestep : -1
        k != 0 && return fill(-one(T), ndim)
        r = @view x[1:pos_dim]
        v = @view x[pos_dim+1:2pos_dim]
        r_meas = @view x_meas[1:pos_dim]
        v_meas = @view x_meas[pos_dim+1:2pos_dim]
        Φ = x[s_idx]
        pos_rows = [T(lambda) * (r[j] - T(r_meas[j])) + (v[j] - T(v_meas[j])) - Φ for j in 1:pos_dim]
        neg_rows = [-(T(lambda) * (r[j] - T(r_meas[j])) + (v[j] - T(v_meas[j]))) - Φ for j in 1:pos_dim]
        if !has_pin
            return vcat(pos_rows, neg_rows)
        end
        tol = T(pos_pin_tol)
        pin_pos_rows = [(r[j] - T(r_meas[j])) - tol for j in 1:pos_dim]
        pin_neg_rows = [-(r[j] - T(r_meas[j])) - tol for j in 1:pos_dim]
        return vcat(pos_rows, neg_rows, pin_pos_rows, pin_neg_rows)
    end
    return PrivateInequality(i; func = func, dim = ndim, is_convex = true)
end

# ============================================================================
# Free-initial mask (§4) — solver-agnostic helper.
# ============================================================================

"""
    dtmpc_free_initial_mask(params, sdims) -> Vector{Vector{Bool}}

Per-agent boolean mask over the AUGMENTED initial-state components (physical
`+` tube widths `+` tube bandwidths) marking every component `true` (free — not
hard-pinned). Per the theory's §4 convention, the initial nominal physical state and the
initial tube widths are decision variables (constrained only by initial containment,
[`build_initial_containment_constraint`](@ref)); the initial bandwidth is likewise free,
constrained only by its `[α_lo, α_max]` box.

Pass the result as `FALCONSolver(...; free_initial=dtmpc_free_initial_mask(params,sdims))`
when constructing the inner FALCON solver (this package has no FALCON
dependency, so the mask is computed here but applied by the caller). Inner
solvers without this notion (e.g. ALGAMES) simply ignore it — the classical
hard pin remains, a valid conservative special case.
"""
function dtmpc_free_initial_mask(params::DTMPCParameters{T}, sdims::Vector{Int}) where {T}
    np = length(sdims)
    return [fill(true, sdims[i] + 2 * n_tube_channels_agent(params, i)) for i in 1:np]
end

# ============================================================================
# Augmented game assembly
# ============================================================================

"""
    build_augmented_game(game, params, uncertainty_models, x_phys_joint, s_cur;
                         discrete=true) -> GameProblem (PD-GNEP)

Assemble the augmented `PDGNEProblem` solved at one receding-horizon window.
`x_phys_joint` is the MEASURED joint physical state (used both to seed the
per-agent frozen-neighbor references for coupled uncertainty AND as the
numeric target of the initial-containment constraints). `s_cur[i]` is agent
`i`'s CARRIED translational-channel width (warm-start seed only — NOT a hard
pin; see [`dtmpc_free_initial_mask`](@ref)). `uncertainty_models[i]` supplies
Δᵢ(xᵢ,x_other) for the TRANSLATIONAL channel only (own-state κ_speed handles
the rest, per channel, via [`TubeChannel`](@ref)).
"""
function build_augmented_game(
    game::GameProblem{T},
    params::DTMPCParameters{T},
    uncertainty_models::Vector{<:Function},
    x_phys_joint::Vector{T},
    s_cur::Vector{T};
    discrete::Bool = true,
    carried_widths::Union{Nothing, Vector{Vector{T}}} = nothing,
    carried_alphas::Union{Nothing, Vector{Vector{T}}} = nothing,
) where {T}
    @assert game.dynamics isa SeparableDynamics "DynamicTubeMPC requires a PD-GNEP (SeparableDynamics) horizon game"
    @assert discrete "DTMPC (theory-faithful path) requires the discrete inner solver (FALCON)"
    np = n_players(game)
    sdims = [state_dim(game, i) for i in 1:np]
    cdims = [control_dim(game, i) for i in 1:np]
    pdyn  = game.dynamics.player_dynamics
    N     = n_steps(game)
    dt    = game.time_horizon.dt

    nt_vec = [n_tube_channels_agent(params, i) for i in 1:np]
    aug_offsets = block_offsets(augmented_state_dims(sdims, nt_vec))

    players = PlayerSpec{T}[]
    for i in 1:np
        x_i = agent_physical_block(x_phys_joint, sdims, i)
        x_other = isempty([j for j in 1:np if j != i]) ? T[] :
                  vcat([agent_physical_block(x_phys_joint, sdims, j) for j in 1:np if j != i]...)
        chans = tube_channels(params, i)
        nch = length(chans)

        # Channel-1 (translational) uses the supplied uncertainty model as its
        # Δ term (LSH β=Δ+D+η); other channels have Δ≡0 (their bound is pure D_r/D_d+η).
        Δ_funcs = Union{Nothing,Function}[k == 1 ? uncertainty_models[i] : nothing for k in 1:nch]
        f_aug = build_augmented_player_dynamics_discrete(
            pdyn[i], sdims[i], cdims[i], chans, x_other, dt, sdims, cdims, N, nt_vec;
            Δ_funcs = Δ_funcs, α_lo = params.α_lo)
        obj_aug = build_augmented_cost(get_objective(game, i), sdims[i], cdims[i];
                                       sdims = sdims, cdims = cdims, N = N, nt_vec = nt_vec)

        priv_i = Any[c for c in game.private_constraints if c.player == i]

        # Per-channel state bounds: 0 ≤ s_k ≤ s_max (width block).
        for k in 1:nch
            push!(priv_i, StateBounds(i; state_offset = sdims[i] + k - 1, state_dim = 1,
                                      lower = [zero(T)], upper = [params.s_max]))
        end
        # Per-channel state bounds: α_lo ≤ α_k ≤ α_max (bandwidth block, appended after the
        # widths). α_lo > 0 keeps the tube ODE stable and the δ̂ division well-posed.
        for k in 1:nch
            push!(priv_i, StateBounds(i; state_offset = sdims[i] + nch + k - 1, state_dim = 1,
                                      lower = [params.α_lo], upper = [params.α_max]))
        end
        # Artificial bandwidth-rate control bound |v_k| ≤ v_max (manuscript Eq. 10). This is
        # the first ControlBounds this assembler emits — before α was co-optimized there was
        # no augmented control at all.
        push!(priv_i, ControlBounds(i; control_offset = cdims[i], control_dim = nch,
                                    lower = fill(-params.v_max, nch),
                                    upper = fill( params.v_max, nch)))
        # Per-channel initial-containment constraint (§4), against the MEASURED slice.
        # Session correction: the TRANSLATIONAL channel uses the sliding-surface l_infty
        # form (matching the componentwise ancillary and the sampled-data tube recursion
        # this channel now uses -- see `build_initial_containment_constraint_sliding`'s
        # docstring); rotational/detection channels are UNCHANGED (Euclidean norm ball,
        # still the correct semantics there -- no componentwise sliding-surface ancillary
        # governs those channels).
        for (k, c) in enumerate(chans)
            if c.kind === :translational
                x_meas_slice = x_i[c.state_offset+1 : c.state_offset+c.state_dim]
                # Session finding: a hard position pin here (pos_pin_tol finite) was tested and
                # confirmed to suppress FALCON's own NOMINAL solution to near-passivity (not an
                # ancillary-execution effect -- the nominal u* itself collapsed), independent of
                # lambda_slide. Reverted to pos_pin_tol=Inf (no position constraint beyond the
                # sliding-surface one) -- see build_initial_containment_constraint_sliding's
                # docstring for the still-open theoretical question (the sliding constraint alone
                # does not independently bound position) and the pending proper fix (a carried
                # position-radius constraint), NOT this hard pin.
                push!(priv_i, build_initial_containment_constraint_sliding(
                    i, sdims[i], k, params.pos_dim, params.λ[i], x_meas_slice))
            else
                x_meas_slice = x_i[c.state_offset+1 : c.state_offset+c.state_dim]
                push!(priv_i, build_initial_containment_constraint(
                    i, sdims[i], k, c.state_offset, c.state_dim, x_meas_slice))
            end
        end
        # Per-channel tightening (Corollary 2 actuator/rate authority, detection floor).
        for (k, c) in enumerate(chans)
            if c.kind === :translational && params.tighten_actuator
                # Convex surrogate (opt-in, default off): linearize the bilinear alpha*s
                # reserve about this window's own carried alpha/s (c=carried_alphas/
                # carried_widths, channel 1 = translational) -- exact at that point, convex
                # everywhere. Falls back to the exact bilinear form when no carried α/s is
                # yet available (e.g. the very first window).
                use_convex = params.tighten_actuator_convex &&
                             carried_widths !== nothing && carried_alphas !== nothing
                if use_convex
                    c_i = max(carried_alphas[i][1], T(1e-8)) / max(carried_widths[i][1], T(1e-8))
                    push!(priv_i, build_tightened_actuator_constraint_convex(
                        i, sdims[i], nch, k, params.pos_dim, params.actuator_max[i],
                        params.trans_input_scale[i], c_i))
                else
                    push!(priv_i, build_tightened_actuator_constraint(
                        i, sdims[i], nch, k, params.pos_dim, params.actuator_max[i],
                        params.trans_input_scale[i]))
                end
            elseif c.kind === :rotational
                # Second-order channel: torque (control) tightening and rate (state)
                # tightening are independent constraints — both may apply.
                if params.tighten_rot_actuator
                    push!(priv_i, build_tightened_rot_actuator_constraint(
                        i, sdims[i], nch, k, params.rot_ctrl_offset, params.rot_dim,
                        params.rot_actuator_max[i], params.J_rot[i]))
                end
                if params.tighten_rot_rate
                    push!(priv_i, build_tightened_rot_rate_constraint(
                        i, sdims[i], k, params.rot_state_offset, params.rot_dim,
                        params.rot_rate_max[i]))
                end
            elseif c.kind === :detection
                push!(priv_i, build_tightened_detection_constraint(
                    i, params.det_state_offset + 1, sdims[i] + k, params.det_limit))
            end
        end

        s_seed = carried_widths !== nothing ? carried_widths[i] : T[c.s₀ for c in chans]
        length(s_seed) == nch || (s_seed = T[c.s₀ for c in chans])
        s_seed = copy(s_seed); s_seed[1] = s_cur[i]
        # Bandwidth seed: the carried closed-loop α if available, else the channel's α₀.
        # Clamped into the box so the window's initial state is never infeasible by
        # construction (the seed comes from a previous window's solution, which respects the
        # box, but a caller-supplied carry could not be assumed to).
        α_seed = carried_alphas !== nothing && length(carried_alphas[i]) == nch ?
                 T[clamp(a, params.α_lo, params.α_max) for a in carried_alphas[i]] :
                 T[clamp(c.α₀, params.α_lo, params.α_max) for c in chans]
        x0_i = vcat(x_i, s_seed, α_seed)
        push!(players, PlayerSpec(i, sdims[i] + 2 * nch, cdims[i] + nch,
                                  x0_i, f_aug, obj_aug, priv_i))
    end

    shared = Any[]
    # User shared constraints, carried through augmentation (see
    # `build_augmented_shared_constraint` above). Skip a user-supplied `ProximityConstraint`
    # on a pair already covered by `constraint_tightening` below, so the collision floor
    # isn't enforced twice (once untightened, once tightened).
    tightened_pairs = params.constraint_tightening ? Set(Set.(params.coupled_pairs)) : Set{Set{Int}}()
    for sc in game.shared_constraints
        sc isa ProximityConstraint && Set(get_players(sc)) in tightened_pairs && continue
        push!(shared, build_augmented_shared_constraint(sc, aug_offsets, sdims, cdims, nt_vec))
    end
    if params.constraint_tightening
        for (i, j) in params.coupled_pairs
            push!(shared, build_tightened_collision_constraint(
                i, j, aug_offsets, sdims[i], sdims[j], params.pos_dim, params.collision_min_distance))
        end
    end

    th = game.time_horizon
    @assert th isa DiscreteTime "DynamicTubeMPC requires a DiscreteTime horizon"
    return PDGNEProblem(players, shared, T(th.tf), T(th.dt))
end

# ============================================================================
# Ancillary boundary-layer sliding control (execution-time)
# ============================================================================

"""
    lsh_boundary_layer_control(ū_i, x_i, x̄_i, s_i, α_i, λ_i, pos_dim; scale) -> u_total_i

LSH boundary-layer sliding ancillary for a 2nd-order translational subsystem.
`α_i` is the channel's CO-OPTIMIZED bandwidth at this step (manuscript Eq. 10) — the caller
passes the carried closed-loop value (`ch_α[i][1][t]`), not a design constant:
`K = scale·α_i·s_i`. Sliding variable `s̃=ṽ+λr̃`; `sat(·)` clamped to `[−1,1]`.
`scale` is `trans_input_scale[i]` (default 1, i.e. acceleration control / B(x)=I) --
the mass for a force-controlled agent (B(x)=I/m), converting the desired acceleration
correction into a force correction (`δF=m·δa`), mirroring how the rotational ancillary
scales by `J_r` to convert an angular-acceleration correction into torque.
"""
function lsh_boundary_layer_control(
    ū_i::AbstractVector{T}, x_i::AbstractVector{T}, x̄_i::AbstractVector{T},
    s_i::T, α_fixed_i::T, λ_i::T, pos_dim::Int; ε::T = T(1e-6), scale::T = one(T),
) where {T}
    u = collect(ū_i)
    Φ = max(s_i, ε)
    K = scale * α_fixed_i * s_i
    @inbounds for k in 1:pos_dim
        r̃ = x_i[k] - x̄_i[k]
        ṽ = x_i[pos_dim+k] - x̄_i[pos_dim+k]
        s̃ = ṽ + λ_i * r̃
        u[k] -= K * clamp(s̃ / Φ, -one(T), one(T))
    end
    return u
end

"""
    lsh_rotational_boundary_layer_control(u_i, x_i, x̄_i, s_r, α_fixed_r, λ_r, J_r,
        rot_state_offset, rot_dim, rot_ctrl_offset) -> u_total_i

SECOND-ORDER boundary-layer ancillary for the rotational channel: `ω` is now a
STATE (not a control) and the control `u` is body TORQUE. The channel's error
state is `[p̃;ω̃]` (attitude MRP error, rate error), covering
`x_i[rot_state_offset+1 : rot_state_offset+2·rot_dim]`.

Mirrors the translational second-order ancillary's sliding surface `s̃=ṽ+λr̃`,
but `p̃` and `ω̃` live in different spaces related by the (attitude-dependent)
MRP kinematics `ṗ=¼B(p)ω` — so `p̃` is first converted to an "equivalent rate"
`ω_equiv = 4·B(p)⁻¹·p̃` (the exact — not small-angle — kinematics; valid at any
attitude, reusing the inversion from the first-order ancillary this replaces)
before combining: `s̃ = ω̃ + λ_r·ω_equiv`. The correction is a torque
`δτ[k] = -J_r[k]·K·sat(s̃[k]/Φ)`, `K=α_r·s_r` (`α_r` the channel's CO-OPTIMIZED bandwidth
at this step, carried in from the closed loop, not a design constant) — scaling by the
(diagonal) inertia converts the desired angular-ACCELERATION correction into torque
(`τ=J·ω̇`), matching how the translational ancillary's `u[k] -= K·sat(s̃/Φ)`
IS the acceleration control directly (mass=1 there; here it must be scaled).

Only `rot_dim==3` (the physical MRP dimension) uses the exact `B(p)` inversion;
other `rot_dim` values fall back to treating `p̃` directly as the rate-error
partner (`s̃=ω̃+λ_r·p̃`), since the `B(p)` formula is MRP-specific.
"""
function lsh_rotational_boundary_layer_control(
    u_i::AbstractVector{T}, x_i::AbstractVector{T}, x̄_i::AbstractVector{T},
    s_r::T, α_fixed_r::T, λ_r::T, J_r::AbstractVector{T},
    rot_state_offset::Int, rot_dim::Int, rot_ctrl_offset::Int;
    ε::T = T(1e-6),
) where {T}
    u = collect(u_i)
    Φ = max(s_r, ε)
    K = α_fixed_r * s_r
    ω̃ = ntuple(k -> x_i[rot_state_offset+rot_dim+k] - x̄_i[rot_state_offset+rot_dim+k], rot_dim)
    if rot_dim == 3
        p1 = x_i[rot_state_offset+1]; p2 = x_i[rot_state_offset+2]; p3 = x_i[rot_state_offset+3]
        p2n = p1^2 + p2^2 + p3^2
        # MRP SHADOW-SET CONSISTENCY: x̄_i's p may sit on the opposite MRP branch from x_i's p
        # (the shadow map p→−p/‖p‖² represents the IDENTICAL physical attitude whenever
        # ‖p‖ nears/crosses 1) even though the two trajectories are physically close — naive
        # subtraction then gives a huge SPURIOUS p̃ (observed: ‖p̃‖ up to 1.42 against a tube of
        # 0.028, tracking a periodic near-180° "flip" pattern, not a real divergence). Pick
        # whichever of x̄'s p or its shadow is CLOSER to x_i's p before differencing.
        p̄1 = x̄_i[rot_state_offset+1]; p̄2 = x̄_i[rot_state_offset+2]; p̄3 = x̄_i[rot_state_offset+3]
        p̄n2 = p̄1^2 + p̄2^2 + p̄3^2
        d_raw = (p1-p̄1)^2 + (p2-p̄2)^2 + (p3-p̄3)^2
        if p̄n2 > ε
            p̄s1 = -p̄1/p̄n2; p̄s2 = -p̄2/p̄n2; p̄s3 = -p̄3/p̄n2
            d_shadow = (p1-p̄s1)^2 + (p2-p̄s2)^2 + (p3-p̄s3)^2
            if d_shadow < d_raw
                p̄1, p̄2, p̄3 = p̄s1, p̄s2, p̄s3
            end
        end
        p̃vec = T[p1-p̄1, p2-p̄2, p3-p̄3]
        # B(p) = (1−p'p)·I + 2·[p]_× + 2·p·p'  (so that ṗ = ¼B(p)ω, matching mrp_pdot)
        B = Matrix{T}(undef, 3, 3)
        B[1,1] = (1-p2n)+2p1*p1; B[1,2] = -2p3+2p1*p2;   B[1,3] =  2p2+2p1*p3
        B[2,1] =  2p3+2p1*p2;    B[2,2] = (1-p2n)+2p2*p2; B[2,3] = -2p1+2p2*p3
        B[3,1] = -2p2+2p1*p3;    B[3,2] =  2p1+2p2*p3;    B[3,3] = (1-p2n)+2p3*p3
        ω_equiv = 4 .* (B \ p̃vec)
        @inbounds for k in 1:3
            s̃k = ω̃[k] + λ_r * ω_equiv[k]
            u[rot_ctrl_offset + k] -= J_r[k] * K * clamp(s̃k / Φ, -one(T), one(T))
        end
    else
        @inbounds for k in 1:rot_dim
            p̃k = x_i[rot_state_offset+k] - x̄_i[rot_state_offset+k]
            s̃k = ω̃[k] + λ_r * p̃k
            u[rot_ctrl_offset + k] -= J_r[k] * K * clamp(s̃k / Φ, -one(T), one(T))
        end
    end
    return u
end

# ============================================================================
# Uncertainty model helpers
# ============================================================================

"""
    create_independent_uncertainty_models(Δ_func, n_agents) -> Vector{Function}

Wrap a single own-state uncertainty bound `Δ(xᵢ) -> scalar` into per-agent
models `Δᵢ(xᵢ,x_other)` that ignore the neighbor argument (separable).
"""
function create_independent_uncertainty_models(Δ_func::Function, n_agents::Int)
    return Function[(x_i, x_other) -> Δ_func(x_i) for _ in 1:n_agents]
end
function create_independent_uncertainty_models(Δ_funcs::Vector{<:Function})
    return Function[(x_i, x_other) -> Δ_funcs[i](x_i) for i in 1:length(Δ_funcs)]
end

"""
    create_coupled_uncertainty_models(Δ_coupled, n_agents) -> Vector{Function}

Per-agent bounds `Δᵢ(xᵢ,x_other)` depending on the (frozen) neighbor
reference (obligation T-1: separability via frozen-neighbor approximation).
"""
function create_coupled_uncertainty_models(Δ_coupled::Function, n_agents::Int)
    return Function[(x_i, x_other) -> Δ_coupled(x_i, x_other) for _ in 1:n_agents]
end

"""
    default_drag_uncertainty(C_d) -> Function

Own-state drag bound `Δ(x) = C_d·‖v‖²`, last two state components = velocity.
NOTE: since this is a NONLINEAR (quadratic) bound, it is NOT affine and hence
incurs the standard SCvx exactness caveat when used as the translational
channel's ADDITIONAL `uncertainty_models` term (re-linearized every SCP
iterate, exact at the converged limit — obligation V-3a). Prefer `κ_speed`
(§11 example own-state δ, smooth-L1 majorant) when a convex-affine-in-the-
smoothed-norm bound suffices.
"""
function default_drag_uncertainty(C_d::T) where {T}
    return function (x)
        n = length(x)
        v = n ≥ 2 ? norm(x[n-1:n]) : norm(x)
        return C_d * v^2
    end
end

function default_coupled_drag_uncertainty(C_d::T) where {T}
    return function (x_i, x_other)
        ni = length(x_i)
        vi = ni ≥ 2 ? norm(x_i[ni-1:ni]) : norm(x_i)
        vo = length(x_other) ≥ 2 ? norm(x_other[end-1:end]) :
             (length(x_other) > 0 ? norm(x_other) : zero(T))
        return C_d * (vi^2 + T(0.5) * vo^2)
    end
end

# ============================================================================
# Coupled-constraint verification helpers (post-hoc checks, NOT the mechanism)
# ============================================================================

tighten_coupled_distance(d_min::T, s_i::T, s_j::T) where {T} = d_min + s_i + s_j

function check_coupled_constraint_satisfaction(
    r_i::AbstractVector{T}, r_j::AbstractVector{T}, s_i::T, s_j::T, d_min::T,
) where {T}
    return norm(r_i .- r_j) ≥ d_min + s_i + s_j
end

# ============================================================================
# Problem / Solution / Solver types
# ============================================================================

"""
    DynamicTubeMPCProblem{T}

Receding-horizon DTMPC problem.

# Fields
- `horizon_game`       : PD-GNEP (`SeparableDynamics`) template for the prediction window
- `x0`                 : initial joint PHYSICAL state
- `n_sim_steps`        : number of closed-loop steps to simulate
- `params`             : `DTMPCParameters`
- `uncertainty_models` : per-agent `Δᵢ(xᵢ,x_other) -> scalar` (translational channel)
- `disturbance`        : optional closed-loop disturbance `d(t)->Vector` added to the
                         physical state each step (default `nothing` ⇒ deterministic)
"""
struct DynamicTubeMPCProblem{T}
    horizon_game::GameProblem{T}
    x0::Vector{T}
    n_sim_steps::Int
    params::DTMPCParameters{T}
    uncertainty_models::Vector{Function}
    disturbance::Union{Nothing, Function}
    truth_dynamics::Union{Nothing, SeparableDynamics}

    function DynamicTubeMPCProblem(
        horizon_game::GameProblem{T}, x0::Vector{T}, n_sim_steps::Int,
        params::DTMPCParameters{T}, uncertainty_models::Vector{<:Function};
        disturbance::Union{Nothing, Function} = nothing,
        truth_dynamics::Union{Nothing, SeparableDynamics} = nothing,
    ) where {T}
        np = n_players(horizon_game)
        @assert horizon_game.dynamics isa SeparableDynamics "DynamicTubeMPC requires a PD-GNEP (SeparableDynamics) horizon game"
        @assert length(x0) == state_dim(horizon_game) "x0 must match total physical state dim"
        @assert n_sim_steps > 0 "n_sim_steps must be positive"
        @assert length(uncertainty_models) == np "need one uncertainty model per agent"
        @assert n_agents(params) == np "params must match number of players"
        new{T}(horizon_game, x0, n_sim_steps, params,
               Vector{Function}(uncertainty_models), disturbance, truth_dynamics)
    end
end

"""
    DynamicTubeMPCSolution{T}

Closed-loop DTMPC solution.

# Fields
- `problem`        : the solved problem
- `X`              : joint physical state trajectory `(n × N_sim+1)`
- `U`              : per-player nominal physical controls `[ (mᵢ × N_sim) ]`
- `s_trajectories` : per-agent TRANSLATIONAL tube-width trajectories `[ (N_sim+1,) ]`
- `costs`          : per-player accumulated physical costs
- `converged`      : all inner solves converged
- `solve_time`     : total inner solve time
- `solver_info`    : diagnostics (per-channel trajectories, rotational/detection widths, etc.)
"""
struct DynamicTubeMPCSolution{T}
    problem::DynamicTubeMPCProblem{T}
    X::Matrix{T}
    U::Vector{Matrix{T}}
    s_trajectories::Vector{Vector{T}}
    costs::Vector{T}
    converged::Bool
    solve_time::Float64
    solver_info::Dict{Symbol, Any}
end

"""
    DynamicTubeMPC{S<:GameSolver} <: GameSolver

Dynamic Tube MPC wrapper. Generic over any **constrained** inner `GameSolver`.
The theory-faithful free-initial convention (§4) requires an inner solver that
supports freeing the initial-state pin (FALCON's `free_initial`, configured by
the caller via [`dtmpc_free_initial_mask`](@ref)); other solvers (e.g. ALGAMES)
keep the classical hard pin (a valid, conservative special case).

`step_callback` (default `nothing`, in which case behaviour is unchanged) is invoked once
per outer receding-horizon step as

    step_callback(; step, t_abs, x, u, s, α, converged)

letting a caller carry its own state between mission steps — a consumed budget, a mission
clock, a peer position for a fixed-point constraint — WITHOUT decomposing the run into
repeated single-step `solve()` calls, which would lose warm-start propagation and force the
tube widths to be hand-reseeded. See the call site for the full argument contract; it is
observational only.

`recovery_callback` (default `nothing`) fires ONLY when a window's primary solve AND the
built-in fresh-warmstart retry have both failed, as

    recovery_callback(; step, current_x, s_cur, carried_widths, carried_alphas, warmstart,
                        aug_game, game) -> Union{Nothing, <:Any}

handed the SAME `warmstart` this window actually entered with (never a failed attempt's own
output) plus enough to rebuild an augmented game via [`build_augmented_game`](@ref) against a
DIFFERENT `game`/`DTMPCParameters` if desired. Return `nothing` to defer to the existing
clean-termination guard, or a solution object satisfying the same duck-typed interface as an
inner solve (`.converged`, [`get_trajectory`](@ref)) — re-validated here with the same
independent finiteness check applied to any primary solve, since a recovery callback is
exactly as untrusted as the inner solver until proven otherwise. This package has NO opinion
on what a recovery attempt should look like (a parameter-homotopy schedule, a relaxed
formulation, ...) — that logic belongs entirely to the caller.

`warmstart_override_callback` (default `nothing`) fires EVERY window, immediately BEFORE the
inner solve (and before `pre_solve_callback`, so that hook's existing "the TRUE warm start
being used this window" contract still holds), as

    warmstart_override_callback(; step, current_x, s_cur, carried_widths, carried_alphas,
                                   default_warmstart) -> Union{Nothing, WarmstartData}

`default_warmstart` is whatever this window would otherwise use (the caller-supplied `WS0`
at step 1, or `warm_start=true`'s own shifted-previous-solution seed thereafter). Returning
`nothing` keeps that default unchanged; returning a `WarmstartData` REPLACES it for this
window's solve only (the package's own `warm_start=true` shift-forward logic still runs
afterward, based on whichever `WarmstartData` was actually used to solve). This exists so a
caller can seed each window from an EXTERNAL reference trajectory (e.g. re-anchoring at the
measured state and forward-rolling a frozen open-loop trajectory's own controls through the
current window dynamics) instead of the package's default "shift the previous window's own
solution" warm start — a genuinely different question from `recovery_callback` (which only
fires on failure) and orthogonal to `warm_start` (which this callback does not disable; it
overrides what feeds the SAME machinery). This package has no opinion on how the replacement
warm start is constructed — that logic belongs entirely to the caller.

`ancillary_callback` (default `nothing`, session addition): fires once PER AGENT at EVERY
physical step (replan or not), immediately after that agent's ancillary correction is
computed, as

    ancillary_callback(; step, agent, is_replan, node, x, xbar, ustar, u_total, s, alpha)

`step` is the 1-based physical step, `agent` the player index, `is_replan` whether THIS step
triggered a fresh solve, `node` the local window node just used (`j` -- 1 on every replan
step), `x`/`xbar` the measured/nominal physical state this step's correction was computed
from, `ustar` the window's own planned control at that node, `u_total` the control ACTUALLY
applied (ustar + ancillary correction, if any), and `s`/`alpha` the carried tube
width/bandwidth used for the correction. Exists to audit, per physical step, whether the
ancillary boundary-layer correction is genuinely recomputed from the CURRENT measured state
on every step (including non-replan steps under `replan_stride>1`) or is stale/zero --
purely observational, mutating the arguments has no effect on the loop.

`replan_stride` (session addition, default `1` — every existing call site is unchanged and
behaves byte-identically): re-solve the window only every `replan_stride` physical steps
instead of every step. On a non-replan step, the LAST solved window's own prediction is
reused at the appropriate local node (no re-optimization), while the ancillary
`:lsh_boundary_layer` feedback correction and the tube width/bandwidth propagation
(`ch_s`/`ch_α`) still run EVERY physical step, exactly as they do on a replan step — nothing
about the executed closed-loop trajectory is open-loop or frozen between replans, only the
underlying nominal plan is held fixed for `replan_stride` steps at a time. When a replan
does occur, the next warm start is built via the SAME `_shift_warmstart`/
`_fix_dtmpc_warmstart_tail` machinery used at `replan_stride=1`, generalized to shift by
`replan_stride` nodes (a genuine multi-step dynamics rollout of the held-last-control tail,
not a naive repeat — see `_shift_trajectories`'s `shift` kwarg). `replan_stride` must not
exceed the window length (`n_steps(horizon_game)`); this is checked at solve time, once the
game is available.
"""
struct DynamicTubeMPC{S <: GameSolver} <: GameSolver
    inner_solver::S
    warm_start::Bool
    verbose_inner::Bool
    params::Union{Nothing, DTMPCParameters}
    discrete_dynamics::Bool
    step_callback::Union{Nothing, Function}
    pre_solve_callback::Union{Nothing, Function}
    recovery_callback::Union{Nothing, Function}
    window_callback::Union{Nothing, Function}
    warmstart_override_callback::Union{Nothing, Function}
    replan_stride::Int
    ancillary_callback::Union{Nothing, Function}
    post_solve_callback::Union{Nothing, Function}

    function DynamicTubeMPC(
        inner_solver::S;
        warm_start::Bool = true,
        verbose_inner::Bool = false,
        params::Union{Nothing, DTMPCParameters} = nothing,
        discrete_dynamics::Bool = true,
        step_callback::Union{Nothing, Function} = nothing,
        pre_solve_callback::Union{Nothing, Function} = nothing,
        recovery_callback::Union{Nothing, Function} = nothing,
        window_callback::Union{Nothing, Function} = nothing,
        warmstart_override_callback::Union{Nothing, Function} = nothing,
        replan_stride::Int = 1,
        ancillary_callback::Union{Nothing, Function} = nothing,
        post_solve_callback::Union{Nothing, Function} = nothing,
    ) where {S <: GameSolver}
        @assert replan_stride >= 1 "replan_stride must be >= 1, got $replan_stride"
        new{S}(inner_solver, warm_start, verbose_inner, params, discrete_dynamics, step_callback, pre_solve_callback, recovery_callback, window_callback, warmstart_override_callback, replan_stride, ancillary_callback, post_solve_callback)
    end
end

function solver_capabilities(::Type{<:DynamicTubeMPC})
    return [:RecedingHorizon, :RobustMPC, :TubeBased, :MultiAgent, :GNEP, :SharedConstraints]
end

# ============================================================================
# Inner-solution joint-state adapter
# ============================================================================

function _joint_augmented_state(inner_sol, np::Int)
    inner_sol.state_trajectory !== nothing && return inner_sol.state_trajectory
    return reduce(vcat, (get_trajectory(inner_sol, i).states for i in 1:np))
end

# Discrete one-step rollout: player dynamics are discrete next-state maps
# `x_{k+1}=f(x_k,u_k,p,k)` (FALCON convention). `k` is the 0-based absolute
# timestep so time-dependent terms align. Passes `x_full` (current joint
# physical state) so a player's TRUTH dynamics can read the ACTUAL opponent
# state (e.g. detection integrating against the real maneuvering target,
# not a stale reference) when it supports doing so.
function _rollout_discrete(dyn::SeparableDynamics, x::AbstractVector,
                           u::AbstractVector, k::Int)
    sdims = dyn.state_dims; cdims = dyn.control_dims
    so = [0; cumsum(sdims)[1:end-1]]
    co = [0; cumsum(cdims)[1:end-1]]
    x_next = similar(x)
    for i in eachindex(sdims)
        xi = @view x[so[i]+1 : so[i]+sdims[i]]
        ui = @view u[co[i]+1 : co[i]+cdims[i]]
        x_next[so[i]+1 : so[i]+sdims[i]] =
            dyn.player_dynamics[i](xi, ui, (player=i, timestep=k, x_full=x), k)
    end
    return x_next
end

# ============================================================================
# Warm-start tail correction (session correctness fix)
# ============================================================================

"""
    _fix_dtmpc_warmstart_tail(ws, inner_sol, cdims, nt_vec, dt, α_lo, α_max) -> WarmstartData

`_shift_warmstart` (rhn.jl `_shift_trajectories`) extends the warm-start tail by one node
using the classical MPC hold-last-control convention: repeat the outgoing window's
terminal control verbatim, then roll the dynamics forward one step. That convention is
correct for PHYSICAL controls, but WRONG for the artificial per-channel bandwidth-rate
control `v` that drives the co-optimized tube bandwidth (manuscript Eq. 10,
`alpha_{k+1}=alpha_k+dt*v_k`): holding a nonzero `v` and re-extending the tail by one more
node on EVERY receding-horizon step compounds `alpha_tail += dt*v_held` across many
shifts, eventually driving the warm-start's OWN late nodes outside the hard box
`alpha_lo<=alpha<=alpha_max` — directly confirmed via a per-node postmortem (window 43/44
of a real solar-game run: warm-start alpha up to 6.8 against alpha_max=5.0) — even though
the CONVERGED solution the warm start derives from already satisfies the box everywhere,
by construction.

Fix: for each agent's translational/rotational tube channels, zero the tail's artificial
rate-control rows (`v_tail=0`, so `alpha_tail=alpha_last`, automatically feasible since the
preceding converged `alpha_last` already is) BEFORE re-rolling the tail state through the
SAME per-player augmented dynamics `_shift_trajectories` already uses — keeping the tube
width `s` dynamically consistent with the corrected `alpha`/physical rows (`s` is NOT
independently clipped or reset; it is whatever the augmented dynamics produce from the
corrected control). Physical control rows are left untouched (still held, per the
classical convention — only the LAST `nt_vec[i]` rows of each agent's control vector are
the artificial rate controls).

Also serves as a standing correctness check: asserts every returned trajectory's
bandwidth-state rows lie in `[α_lo-ε, α_max+ε]` at every node — a DTMPC-generated warm
start should never hand FALCON one of its own co-optimized states outside its own hard
box. This is a permanent invariant, not a debug-only toggle: cheap (O(N) per window) and
exactly the class of defect that motivated this fix.

`shift` (session addition, default `1` — every existing call site is unchanged): the number
of NEWLY-APPENDED tail nodes to fix (i.e. however many nodes `_shift_warmstart` was itself
called with `shift=` for). For `shift>1` (DTMPC's `replan_stride` option), the appended tail
holds the SAME artificial rate-control `v` for all `shift` appended nodes before this fix
runs — zeroing only the last one would leave `alpha` compounding uncorrected on the other
`shift-1` appended nodes. Fixes all `shift` appended nodes by zeroing `v` on each and
re-rolling the tail states sequentially from the last NON-appended node, so `alpha_tail`
stays pinned at `alpha_last` (not drifting) across the whole coarser-cadence gap.
"""
function _fix_dtmpc_warmstart_tail(ws::WarmstartData{T}, inner_sol, cdims::Vector{Int},
                                    nt_vec::Vector{Int}, dt::T, α_lo::T, α_max::T;
                                    ε::T = T(1e-6), shift::Int = 1) where {T}
    ws.trajectories === nothing && return ws
    dyn = inner_sol.game.dynamics
    fixed_trajs = map(ws.trajectories) do tr
        N = size(tr.controls, 2)
        i = findfirst(obj -> obj.player_id == tr.player_id, inner_sol.game.objectives)
        nch = nt_vec[i]; cd = cdims[i]
        m = clamp(shift, 0, N)   # number of appended tail nodes actually present to fix
        if N > 1 && nch > 0 && m > 0
            new_controls = copy(tr.controls)
            for k in (N-m+1):N
                new_controls[cd+1:cd+nch, k] .= zero(T)   # v_tail = 0 -> alpha_tail = alpha_last
            end
            if dyn isa SeparableDynamics
                fᵢ = dyn.player_dynamics[i]
                new_states = copy(tr.states)
                for k in (N-m+1):N
                    new_states[:, k+1] = fᵢ(new_states[:, k], new_controls[:, k], nothing, (k - 1) * dt)
                end
            else
                new_states = tr.states
            end
            tr = Trajectory{T}(tr.player_id, new_states, new_controls, copy(tr.times), tr.cost)
        end
        n_i = size(tr.states, 1) - 2 * nch   # physical state dim for this agent
        for k in 1:nch
            α_row = n_i + nch + k
            α_vals = @view tr.states[α_row, :]
            @assert all(v -> α_lo - ε <= v <= α_max + ε, α_vals) "DTMPC warm-start tail fix: " *
                "player $(tr.player_id) channel $k alpha out of box [$(α_lo),$(α_max)]: " *
                "$(extrema(α_vals))"
        end
        tr
    end
    return WarmstartData{T}(fixed_trajs, ws.strategy, ws.dual_variables, ws.solver_specific)
end

# ============================================================================
# Main receding-horizon loop
# ============================================================================

function _solve_dtmpc(
    prob::DynamicTubeMPCProblem{T},
    solver::DynamicTubeMPC,
    warmstart::Union{Nothing, WarmstartData},
    verbose::Bool,
) where {T}
    game = prob.horizon_game
    @assert game.dynamics isa SeparableDynamics "DynamicTubeMPC requires a PD-GNEP (SeparableDynamics) horizon game"
    @assert solver.replan_stride <= n_steps(game) "replan_stride ($(solver.replan_stride)) must not exceed the window length ($(n_steps(game)))"

    np    = n_players(game)
    n     = state_dim(game)
    N_sim = prob.n_sim_steps
    params = solver.params !== nothing ? solver.params : prob.params

    sdims = [state_dim(game, i) for i in 1:np]
    cdims = [control_dim(game, i) for i in 1:np]
    nt_vec = [n_tube_channels_agent(params, i) for i in 1:np]
    chans_all = [tube_channels(params, i) for i in 1:np]
    aug_offsets = block_offsets(augmented_state_dims(sdims, nt_vec))
    phys_offsets = block_offsets(sdims)

    th = game.time_horizon
    times_w = collect(range(zero(T), T(th.tf), length = n_steps(game) + 1))

    X = Matrix{T}(undef, n, N_sim + 1)
    X[:, 1] = prob.x0
    U_players = [Matrix{T}(undef, cdims[i], N_sim) for i in 1:np]

    # Separately-tracked TRULY-APPLIED control (nominal ū_i plus the ancillary
    # boundary-layer correction κ_i, when `ancillary=:lsh_boundary_layer`; identical
    # to `U_players` otherwise). `U_players` intentionally stores the NOMINAL
    # (pre-ancillary) control ū_i — matching the cost accounting below, which is
    # evaluated on the planned trajectory, not the realized one — so this is a
    # separate array rather than a change to `U_players`'s existing meaning.
    U_applied = [Matrix{T}(undef, cdims[i], N_sim) for i in 1:np]

    # Separately-tracked NOMINAL joint state, used as the WARM-START SEED for
    # each window's (free) initial nominal — NOT a hard pin (§4). With
    # nominal-tracking ancillary (`ancillary=:lsh_boundary_layer`) the plan is
    # seeded from this carried nominal so the ancillary correction is genuinely
    # nonzero (x−x̄≠0); replanning from the raw measurement would make it a
    # no-op. The CONTAINMENT constraint (added in `build_augmented_game`)
    # anchors x̄_{t|t} to the true measurement regardless of this seed.
    X̄ = Matrix{T}(undef, n, N_sim + 1)
    X̄[:, 1] = prob.x0
    nominal_mode = params.ancillary === :lsh_boundary_layer

    # Per-agent, per-channel tube-width trajectories `ch_s[i][k]`.
    ch_s = [[Vector{T}(undef, N_sim + 1) for _ in 1:nt_vec[i]] for i in 1:np]
    for i in 1:np, k in 1:nt_vec[i]
        ch_s[i][k][1] = chans_all[i][k].s₀
    end
    s_traj = [ch_s[i][1] for i in 1:np]   # translational-channel alias

    # Per-agent, per-channel tube-BANDWIDTH trajectories `ch_α[i][k]` (co-optimized,
    # manuscript Eq. 10 — mirrors `ch_s` exactly, seeded from each channel's α₀).
    ch_α = [[Vector{T}(undef, N_sim + 1) for _ in 1:nt_vec[i]] for i in 1:np]
    for i in 1:np, k in 1:nt_vec[i]
        ch_α[i][k][1] = chans_all[i][k].α₀
    end
    α_traj = [ch_α[i][1] for i in 1:np]   # translational-channel alias

    total_costs = zeros(T, np)
    coupled_violations = Tuple{Int, Int, Int}[]

    total_time = 0.0
    n_converged = 0
    all_converged = true
    terminated_early = false
    termination_reason = nothing
    N_completed = N_sim

    verbose && @info "DTMPC: $(N_sim)-step closed loop, $(np) agents, augmented PD-GNEP" inner = typeof(solver.inner_solver)

    # `replan_stride` bookkeeping (session addition, no-op path when replan_stride==1): a
    # step is a REPLAN step when `(t - last_replan_t) % replan_stride == 0` (true for t=1,
    # then every `replan_stride`-th step thereafter). On a non-replan step, `cached_inner_sol`
    # / `cached_joint_states` (the last REPLAN's own converged solution/prediction) are reused
    # at local window node `j = t - last_replan_t + 1` instead of the hardcoded node 1 a
    # replan step always reads — the ancillary correction and tube propagation below are
    # otherwise byte-identical code, just parameterized by `j` (which equals 1 on every
    # replan step, so `replan_stride=1` reproduces the exact original per-step behavior).
    cached_inner_sol = nothing
    cached_joint_states = nothing
    last_replan_t = 0

    for t in 1:N_sim
        current_x = X[:, t]                              # measured state
        is_replan_step = last_replan_t == 0 || (t - last_replan_t) % solver.replan_stride == 0
        plan_x    = nominal_mode ? X̄[:, t] : current_x    # warm-start seed for the window
        s_cur = [ch_s[i][1][t] for i in 1:np]
        carried = [T[ch_s[i][k][t] for k in 1:nt_vec[i]] for i in 1:np]
        carried_α = [T[ch_α[i][k][t] for k in 1:nt_vec[i]] for i in 1:np]

        local inner_sol, joint_states
        if is_replan_step
            aug_game = build_augmented_game(game, params, prob.uncertainty_models,
                                            current_x, s_cur;    # containment anchors to the MEASUREMENT
                                            discrete = solver.discrete_dynamics,
                                            carried_widths = carried,
                                            carried_alphas = carried_α)

            # ── Warm-start override hook ─────────────────────────────────────────────────────
            # Fires BEFORE pre_solve_callback (so that hook still observes the TRUE warm start
            # about to be used) and before the inner solve. Lets a caller replace this window's
            # warm start entirely (e.g. seed from an external reference trajectory) instead of
            # the package's own default (caller-supplied WS0 at step 1, shifted-previous-
            # solution thereafter when warm_start=true). See DynamicTubeMPC's docstring.
            if solver.warmstart_override_callback !== nothing
                override = solver.warmstart_override_callback(;
                    step = t, current_x = current_x, s_cur = s_cur,
                    carried_widths = carried, carried_alphas = carried_α,
                    default_warmstart = warmstart,
                )
                override !== nothing && (warmstart = override)
            end

            # ── Pre-solve diagnostic hook ────────────────────────────────────────────────────
            # Fires immediately BEFORE the inner solve, with the exact object about to be
            # handed to it: the augmented game, the TRUE warm start being used this window
            # (already shifted from the previous window's own converged solution when
            # warm_start=true -- NOT reconstructable from any saved post-hoc trajectory), the
            # measured state, and the carried tube widths/bandwidths. Exists so a caller can
            # inspect a specific window's actual input (e.g. compute chi_seed via
            # _compute_infeasibility on the real warm start, not an approximation) without
            # needing to re-derive it. Purely observational, mirrors step_callback's contract.
            if solver.pre_solve_callback !== nothing
                solver.pre_solve_callback(;
                    step = t, current_x = current_x, s_cur = s_cur,
                    carried_widths = carried, carried_alphas = carried_α,
                    warmstart = warmstart, aug_game = aug_game,
                )
            end

            t0 = time()
            inner_sol = solve(aug_game, solver.inner_solver;
                              warmstart = warmstart, verbose = solver.verbose_inner,
                              check_compatibility = false)
            total_time += time() - t0

            # ── MANDATORY failure-path guard (session correctness fix) ─────────────────────────
            # A non-converged inner_sol must NEVER drive the ancillary boundary-layer controllers
            # or advance the measured state: doing so unconditionally previously crashed downstream
            # (lsh_rotational_boundary_layer_control hit a matrix containing NaN/Inf) after a window
            # failed to converge -- a non-converged solve's returned trajectory carries NO
            # correctness guarantee (it may be the warm start echoed back unchanged, a partially
            # -stepped trust-region iterate, or worse). `inner_sol.converged` is checked together
            # with an INDEPENDENT `allfinite` check on the actual states/controls about to be used,
            # since no inner solver's internal bookkeeping is assumed trustworthy on its own.
            _finite_sol(sol) = all(i -> all(isfinite, get_trajectory(sol, i).controls) &&
                                         all(isfinite, get_trajectory(sol, i).states), 1:np)
            ok = inner_sol.converged && _finite_sol(inner_sol)
            if !ok
                verbose && @warn "DTMPC: inner solve did not converge/was non-finite at step $t — retrying once from a fresh default trial (no warm start, reset trust region)"
                inner_sol = solve(aug_game, solver.inner_solver;
                                  warmstart = nothing, verbose = solver.verbose_inner,
                                  check_compatibility = false)
                ok = inner_sol.converged && _finite_sol(inner_sol)
            end
            # ── Recovery hook (session extension, generic — DTMPC stays ignorant of WHY a caller
            # can recover a failed window): fires ONLY after both the primary solve and the fresh-
            # warmstart retry above have failed. Handed the SAME `warmstart` this window actually
            # entered with (the last-known-good shifted warm start — never the failed attempt's own
            # output) plus the augmented game and measured/tube state, so a caller can retry with
            # e.g. a parameter-homotopy wrapper (a terminal-cost continuation schedule, a relaxed
            # formulation, etc.) built entirely OUTSIDE this package. Must return either `nothing`
            # (defer to the existing clean-termination guard below) or an object satisfying the
            # SAME duck-typed interface as `inner_sol` (`.converged`, `get_trajectory(sol, i)`) —
            # re-checked here with the SAME independent finiteness test, since a recovery callback
            # is exactly as untrusted as the primary solver until proven otherwise. ────────────────
            if !ok && solver.recovery_callback !== nothing
                recovered = solver.recovery_callback(;
                    step = t, current_x = current_x, s_cur = s_cur,
                    carried_widths = carried, carried_alphas = carried_α,
                    warmstart = warmstart, aug_game = aug_game, game = game,
                )
                if recovered !== nothing && recovered.converged && _finite_sol(recovered)
                    inner_sol = recovered
                    ok = true
                end
            end
            if !ok
                all_converged = false
                terminated_early = true
                N_completed = t - 1
                termination_reason = "inner solve failed (non-converged or non-finite) at step $t, even after a fresh retry"
                @warn "DTMPC: $termination_reason — terminating the simulation cleanly at step $N_completed rather than propagating an untrustworthy trajectory into the ancillary controller or the next window"
                break
            end

            n_converged += 1
            last_replan_t = t

            # ── Post-solve diagnostic hook (session addition, purely observational) ─────────────
            # Fires once per REPLAN step (never on a cached/reused non-replan step), immediately
            # after this window's inner_sol is confirmed converged/finite. Exists so a caller can
            # record per-step solver performance (SCP iteration count, wall-clock solve time --
            # both already computed by the inner solver and carried on `inner_sol`, e.g.
            # `.iterations`/`.solve_time` on a FALCON GNEPSolution, but otherwise discarded by this
            # loop) without needing to fork or re-implement the receding-horizon loop itself.
            if solver.post_solve_callback !== nothing
                solver.post_solve_callback(; step = t, inner_sol = inner_sol)
            end

            joint_states = _joint_augmented_state(inner_sol, np)
            cached_inner_sol = inner_sol
            cached_joint_states = joint_states
        else
            # Non-replan step (`replan_stride>1` only): reuse the last REPLAN's own converged
            # prediction, no re-solve. The ancillary correction/tube propagation below reads
            # local window node `j` (computed next) instead of node 1 -- everything else in
            # this branch is identical to the replan branch's post-solve bookkeeping.
            inner_sol = cached_inner_sol
            joint_states = cached_joint_states
        end
        j = t - last_replan_t + 1   # local window node for THIS physical step (1 on every replan step)

        u_total_agents = Vector{Vector{T}}(undef, np)
        window_controls_local = Vector{Matrix{T}}(undef, np)   # session addition, see window_callback below
        for i in 1:np
            traj = get_trajectory(inner_sol, i)
            window_controls_local[i] = traj.controls[1:cdims[i], :]
            ū_i = traj.controls[1:cdims[i], j]
            x_i = agent_physical_block(current_x, sdims, i)
            nch = nt_vec[i]; chans = chans_all[i]

            x̄_i, _, _         = extract_agent_augmented(joint_states[:, j], aug_offsets, sdims, nt_vec, i)
            x̄_next_i, s2v, α2v = extract_agent_augmented(joint_states[:, j+1], aug_offsets, sdims, nt_vec, i)
            for k in 1:nch
                ch_s[i][k][t + 1] = max(s2v[k], zero(T))
                ch_α[i][k][t + 1] = clamp(α2v[k], params.α_lo, params.α_max)
            end

            if params.ancillary === :lsh_boundary_layer
                # Ancillary gain K=α·s now uses the CARRIED closed-loop α (co-optimized),
                # not the channel's fixed design constant.
                u_tmp = lsh_boundary_layer_control(ū_i, x_i, x̄_i, ch_s[i][1][t],
                                                   ch_α[i][1][t], params.λ[i], params.pos_dim;
                                                   scale = params.trans_input_scale[i])
                if nch ≥ 2 && chans[2].kind === :rotational
                    u_tmp = lsh_rotational_boundary_layer_control(
                        u_tmp, x_i, x̄_i, ch_s[i][2][t], ch_α[i][2][t], params.λ_r[i], params.J_rot[i],
                        params.rot_state_offset, params.rot_dim, params.rot_ctrl_offset)
                end
                u_total_agents[i] = u_tmp
            else
                u_total_agents[i] = ū_i
            end

            if solver.ancillary_callback !== nothing
                solver.ancillary_callback(; step = t, agent = i, is_replan = is_replan_step,
                                           node = j, x = x_i, xbar = x̄_i, xbar_next = x̄_next_i, ustar = ū_i,
                                           u_total = u_total_agents[i], s = ch_s[i][1][t], s_next = max(s2v[1], zero(T)),
                                           alpha = ch_α[i][1][t])
            end

            U_players[i][:, t] = ū_i
            U_applied[i][:, t] = u_total_agents[i]
            obj = get_objective(game, i)
            total_costs[i] += obj.scaling * evaluate_stage_cost(obj.stage_cost, x_i, ū_i, nothing, t)
        end

        u_total_joint = vcat(u_total_agents...)
        truth_dyn = prob.truth_dynamics === nothing ? game.dynamics : prob.truth_dynamics
        x_next = _rollout_discrete(truth_dyn, current_x, u_total_joint, t - 1)
        if prob.disturbance !== nothing
            # Backward-compatible: a `d(t)` disturbance is a pure function of time;
            # a `d(x, t)` disturbance additionally sees the current joint physical
            # state, so it can inject STATE-DEPENDENT perturbations (e.g. an
            # attitude disturbance shaped by the full MRP kinematics δp=¼B(p)δω·dt).
            dstep = applicable(prob.disturbance, current_x, t) ?
                        prob.disturbance(current_x, t) : prob.disturbance(t)
            x_next = x_next .+ dstep
        end
        X[:, t + 1] = x_next

        # Advance the carried NOMINAL SEED: the window's own PREDICTED next-step physical
        # state at local node `j+1` (disturbance-free), used only to warm-start the NEXT
        # window (never a hard pin — see the loop-top comment). Node `j+1` (not hardcoded 2)
        # so this stays correct on a non-replan step, which reads further into the SAME
        # cached window's own prediction than node 2.
        if nominal_mode
            for i in 1:np
                x̄n_i = joint_states[aug_offsets[i]+1 : aug_offsets[i]+sdims[i], j+1]
                X̄[phys_offsets[i]+1 : phys_offsets[i]+sdims[i], t + 1] = x̄n_i
            end
        else
            X̄[:, t + 1] = x_next
        end

        for (i2, j2) in params.coupled_pairs
            r_i = agent_physical_block(X[:, t + 1], sdims, i2)[1:params.pos_dim]
            r_j = agent_physical_block(X[:, t + 1], sdims, j2)[1:params.pos_dim]
            if !check_coupled_constraint_satisfaction(r_i, r_j, ch_s[i2][1][t+1], ch_s[j2][1][t+1],
                                                      params.collision_min_distance)
                push!(coupled_violations, (t, i2, j2))
            end
        end

        # ── Per-step callback hook ────────────────────────────────────────────────────────
        # Fires ONCE per outer receding-horizon step, after the measured state `X[:,t+1]` and
        # the applied controls are final, and before the warm start is shifted forward.
        #
        # It exists because this loop is otherwise a closed box: a caller that needs to carry
        # its own bookkeeping between mission steps — a consumed detection budget, a mission
        # clock, a peer position for a fixed-point constraint — previously had no hook and had
        # to decompose the run into repeated `n_sim_steps=1` `solve()` calls instead. That
        # decomposition is expensive and lossy: `inner_sol` is not exposed on the public
        # solution struct, so warm-start propagation (below) cannot be reconstructed across
        # separate calls, and the evolved tube widths in `ch_s`/`ch_α` have to be hand-reseeded
        # from `DTMPCParameters` each time. Both of those live correctly inside this loop.
        #
        # Signature — keyword arguments, so it can be extended without breaking callers:
        #   step_callback(; step, t_abs, x, u, s, α, converged)
        # `step` is the 1-based outer step just completed, `t_abs = (step-1)·dt` the absolute
        # mission time at the START of that step (matching the convention `_rollout_discrete`
        # is handed), `x` the measured joint physical state AFTER the step, `u` the per-player
        # applied physical controls, `s`/`α` the per-agent per-channel tube widths/bandwidths
        # carried into the next step, and `converged` this window's inner-solve status.
        #
        # Purely observational by contract: mutating `x`, `u`, `s`, or `α` here does NOT feed
        # back into the loop, which reads its own arrays. Side effects on the CALLER's state
        # (Refs the game's constraints close over, logs) are the intended use.
        if solver.step_callback !== nothing
            solver.step_callback(;
                step      = t,
                t_abs     = (t - 1) * game.time_horizon.dt,
                x         = X[:, t + 1],
                u         = [U_players[i][:, t] for i in 1:np],
                s         = [T[ch_s[i][k][t + 1] for k in 1:nt_vec[i]] for i in 1:np],
                α         = [T[ch_α[i][k][t + 1] for k in 1:nt_vec[i]] for i in 1:np],
                converged = inner_sol.converged,
            )
        end

        # ── Per-window FULL-prediction callback hook (session addition, SEPARATE from
        # step_callback so no existing caller's fixed-keyword-list callback function can
        # break — several jamming-game scripts define `step_callback`/`on_step` WITHOUT a
        # catch-all `kwargs...`, so adding a new kwarg to THAT invocation would have been a
        # real regression). Fires with this window's OWN just-solved FULL trajectory, per
        # agent, physical rows only (`n_x_i × (N_win+1)`, local node 1 = the measurement this
        # window solved from, local node j = x*_{t+j-1|t}) — sliced directly from
        # `joint_states` (already computed above for the node-1/node-2 bookkeeping; no
        # duplicate solve). Exists because callers previously had no way to inspect a
        # window's own PREDICTED trajectory beyond one step ahead (`:nominal_X` only ever
        # stitches node 2 of each window) — needed to answer "does THIS window's own
        # prediction still contain a future event the stitched one-step-ahead view can't
        # show" (e.g. a delayed reacquisition maneuver) without any additional solve.
        # `window_controls` (session addition): this window's OWN just-solved FULL control
        # sequence, per agent, physical rows only (`n_u_i × N_win`) -- captured from
        # `get_trajectory(inner_sol,i).controls` during the per-agent loop above (no
        # duplicate solve). Needed to compare what an EARLIER window planned to apply at a
        # LATER local node against what a subsequent window actually re-optimized and
        # applied at that same absolute step (Δu_k = u*_{k|k} - u*_{k|k'} for k'<k).
        #
        # `window_aug_node1` (session addition): this window's OWN converged AUGMENTED
        # (tube-width `s` + bandwidth `alpha`) node-1 values, per agent, in that order
        # (`s_chan1,...,s_chanNch,alpha_chan1,...,alpha_chanNch`) — i.e. Phi_{i,t|t} and
        # alpha_{i,t|t}, the window's own FREE initial-containment decision (per
        # `dtmpc_free_initial_mask`, EVERY augmented initial-state component is free, not
        # hard-pinned; only `build_initial_containment_constraint` ties it to the
        # measurement). Neither `window_states` (physical rows only) nor `step_callback`'s
        # `s`/`alpha` (both anchored to node 2, i.e. gamma_{t+1|t}/alpha_{t+1|t} — see that
        # callback's own call site) expose this node-1 quantity; needed to test the
        # INITIAL-containment constraint's own geometry directly, separately from the
        # already-exposed one-step-ahead prediction.
        #
        # Purely observational by contract, same as `step_callback`.
        if solver.window_callback !== nothing
            window_states = [joint_states[aug_offsets[i]+1 : aug_offsets[i]+sdims[i], :] for i in 1:np]
            window_aug_node1 = [joint_states[aug_offsets[i]+sdims[i]+1 : aug_offsets[i]+sdims[i]+2*nt_vec[i], 1] for i in 1:np]
            solver.window_callback(; step = t, window_states = window_states,
                                     window_controls = window_controls_local,
                                     window_aug_node1 = window_aug_node1,
                                     converged = inner_sol.converged)
        end

        # Only propagate the shifted warm start when the inner solve actually CONVERGED.
        # A non-converged `inner_sol` is, by definition, not a trusted trajectory — feeding
        # it forward as the next window's starting point can compound rather than recover
        # (direct diagnosis on a live game: once one window failed to converge, EVERY
        # subsequent window failed too, for the rest of the run, whereas an isolated re-solve
        # of a similar state from FALCON's own fresh default trial often still succeeded).
        # Falling back to `nothing` lets FALCON rebuild its own default trial (rollout from
        # the newly MEASURED state) each time a window fails, giving the next window a
        # clean, well-defined starting point instead of inheriting whatever the failed SCP
        # run left behind.
        # Only rebuild the warm start on an ACTUAL replan step (session addition for
        # `replan_stride`): on a non-replan step `inner_sol` is the CACHED window, unchanged
        # since the last replan, so recomputing a shift-by-1 from it here every skip-step
        # would be wasted work AND wrong (it would overwrite the correctly shift-by-
        # `replan_stride` warm start built at the last replan with a stale shift-by-1 of the
        # SAME cached solution). `replan_stride=1` makes `is_replan_step` true every step,
        # so this reproduces the original per-step behavior exactly.
        if is_replan_step
            if solver.warm_start && inner_sol.converged
                warmstart = _fix_dtmpc_warmstart_tail(_shift_warmstart(inner_sol; shift=solver.replan_stride), inner_sol,
                                                       cdims, nt_vec, T(game.time_horizon.dt),
                                                       params.α_lo, params.α_max; shift=solver.replan_stride)
            else
                warmstart = nothing
            end
        end
    end

    # ── Clean truncation on early termination (failure-path guard above) ───────────────────
    # All the per-step arrays were pre-allocated for the FULL N_sim; on an early `break` the
    # columns/entries beyond N_completed+1 are uninitialized. Slice everything down to what
    # was actually simulated so the returned solution never exposes garbage, and rebind N_sim
    # so the unchanged tail code below (which reads N_sim/`X[:, end]`/etc.) is correct as-is.
    if terminated_early
        N_sim = N_completed
        X  = X[:, 1:N_sim+1]
        X̄  = X̄[:, 1:N_sim+1]
        U_players = [U_players[i][:, 1:N_sim] for i in 1:np]
        U_applied = [U_applied[i][:, 1:N_sim] for i in 1:np]
        ch_s = [[ch_s[i][k][1:N_sim+1] for k in 1:nt_vec[i]] for i in 1:np]
        ch_α = [[ch_α[i][k][1:N_sim+1] for k in 1:nt_vec[i]] for i in 1:np]
        s_traj = [ch_s[i][1] for i in 1:np]
        α_traj = [ch_α[i][1] for i in 1:np]
        coupled_violations = filter(v -> v[1] <= N_sim, coupled_violations)
    end

    for i in 1:np
        obj = get_objective(game, i)
        x_i_final = agent_physical_block(X[:, end], sdims, i)
        total_costs[i] += obj.scaling * evaluate_terminal_cost(obj.terminal_cost, x_i_final, nothing)
    end

    verbose && @info "DTMPC complete" total_time = round(total_time, digits = 3) converged = all_converged coupled_violations = length(coupled_violations)

    rot_idx(i) = (nt_vec[i] ≥ 2 && chans_all[i][2].kind === :rotational) ? 2 : 0
    det_idx(i) = params.detection_tube[i] ? nt_vec[i] : 0
    s_r_traj = [rot_idx(i) > 0 ? ch_s[i][rot_idx(i)] : fill(zero(T), N_sim + 1) for i in 1:np]
    s_d_traj = [det_idx(i) > 0 ? ch_s[i][det_idx(i)] : fill(zero(T), N_sim + 1) for i in 1:np]
    α_r_traj = [rot_idx(i) > 0 ? ch_α[i][rot_idx(i)] : fill(zero(T), N_sim + 1) for i in 1:np]

    return DynamicTubeMPCSolution{T}(
        prob, X, U_players, s_traj, total_costs,
        all_converged, total_time,
        Dict{Symbol, Any}(
            :n_steps_simulated => N_sim,
            :n_inner_converged => n_converged,
            :terminated_early => terminated_early,
            :termination_reason => termination_reason,
            :inner_solver_type => typeof(solver.inner_solver),
            :tube_parameters => params,
            :coupled_pairs => params.coupled_pairs,
            :coupled_constraint_violations => coupled_violations,
            :all_coupled_constraints_satisfied => isempty(coupled_violations),
            :nominal_X => X̄,
            :nominal_mode => nominal_mode,
            :applied_U => U_applied,
            :rotational_tube => copy(params.rotational_tube),
            :detection_tube => copy(params.detection_tube),
            :s_r_trajectories => s_r_traj,
            :s_d_trajectories => s_d_traj,
            :channel_s_trajectories => ch_s,
            :n_tube_channels => nt_vec,
            # Co-optimized bandwidth (manuscript Eq. 10, α̇=v) — additive, mirrors ch_s/s_r_traj.
            :α_trajectories => α_traj,
            :α_r_trajectories => α_r_traj,
            :channel_α_trajectories => ch_α,
        ),
    )
end

"""
    dtmpc_bandwidth_trajectory(sol, i; channel=1) -> Vector

Convenience accessor for agent `i`'s co-optimized bandwidth trajectory `α(t)` on the given
channel (`1`=translational, `2`=rotational if enabled, …), i.e.
`sol.solver_info[:channel_α_trajectories][i][channel]`.
"""
dtmpc_bandwidth_trajectory(sol, i::Int; channel::Int = 1) =
    sol.solver_info[:channel_α_trajectories][i][channel]

# ============================================================================
# solve / _solve interface
# ============================================================================

function DifferentialGamesBase.solve(
    prob::DynamicTubeMPCProblem{T}, solver::DynamicTubeMPC;
    warmstart::Union{Nothing, WarmstartData} = nothing, verbose::Bool = false, kwargs...,
) where {T}
    return _solve_dtmpc(prob, solver, warmstart, verbose)
end

function _solve(
    prob::DynamicTubeMPCProblem{T}, solver::DynamicTubeMPC,
    warmstart::Union{Nothing, WarmstartData}, verbose::Bool,
) where {T}
    return _solve_dtmpc(prob, solver, warmstart, verbose)
end

# ============================================================================
# Accessors
# ============================================================================

dtmpc_state_trajectory(sol::DynamicTubeMPCSolution) = sol.X
dtmpc_control_trajectory(sol::DynamicTubeMPCSolution, i::Int) = sol.U[i]
dtmpc_tube_trajectory(sol::DynamicTubeMPCSolution, i::Int = 1) = sol.s_trajectories[i]
dtmpc_total_cost(sol::DynamicTubeMPCSolution, i::Int) = sol.costs[i]
dtmpc_n_players(sol::DynamicTubeMPCSolution) = n_players(sol.problem.horizon_game)
dtmpc_state_dim(sol::DynamicTubeMPCSolution) = state_dim(sol.problem.horizon_game)
dtmpc_agent_state_dim(sol::DynamicTubeMPCSolution, i::Int) = state_dim(sol.problem.horizon_game, i)
dtmpc_simulation_steps(sol::DynamicTubeMPCSolution) = size(sol.X, 2) - 1
dtmpc_coupled_pairs(sol::DynamicTubeMPCSolution) = sol.problem.params.coupled_pairs

DifferentialGamesBase.get_state_trajectory(sol::DynamicTubeMPCSolution) = sol.X
function DifferentialGamesBase.get_control_trajectory(sol::DynamicTubeMPCSolution, i::Int)
    np = n_players(sol.problem.horizon_game)
    (1 ≤ i ≤ np) || error("player index $i out of range 1:$np")
    return sol.U[i]
end
function DifferentialGamesBase.get_total_cost(sol::DynamicTubeMPCSolution, i::Int)
    np = n_players(sol.problem.horizon_game)
    (1 ≤ i ≤ np) || error("player index $i out of range 1:$np")
    return sol.costs[i]
end

function validate_coupled_constraints(sol::DynamicTubeMPCSolution{T}) where {T}
    return get(sol.solver_info, :all_coupled_constraints_satisfied, true)
end
all_coupled_constraints_satisfied(sol::DynamicTubeMPCSolution) =
    get(sol.solver_info, :all_coupled_constraints_satisfied, true)
get_coupled_constraint_violations(sol::DynamicTubeMPCSolution) =
    get(sol.solver_info, :coupled_constraint_violations, Tuple{Int,Int,Int}[])

# ============================================================================
# Display
# ============================================================================

function Base.show(io::IO, prob::DynamicTubeMPCProblem{T}) where {T}
    np = n_players(prob.horizon_game)
    print(io, "DynamicTubeMPCProblem{$T}:")
    print(io, "\n  n_players: ", np)
    print(io, "\n  total_state_dim: ", state_dim(prob.horizon_game))
    print(io, "\n  n_sim_steps: ", prob.n_sim_steps)
    print(io, "\n  constraint_tightening: ", prob.params.constraint_tightening)
    print(io, "\n  collision_min_distance: ", prob.params.collision_min_distance)
    print(io, "\n  n_coupled_pairs: ", length(prob.params.coupled_pairs))
    for i in 1:min(np, 3)
        print(io, "\n  Agent $i: s₀=", prob.params.s₀[i], " α_fixed=", prob.params.α_fixed[i],
              " D=", prob.params.D[i], " η=", prob.params.η[i])
    end
    np > 3 && print(io, "\n  ... ($(np) agents)")
end

function Base.show(io::IO, sol::DynamicTubeMPCSolution{T}) where {T}
    np = n_players(sol.problem.horizon_game)
    print(io, "DynamicTubeMPCSolution{$T}:")
    print(io, "\n  converged: ", sol.converged)
    print(io, "\n  solve_time: ", round(sol.solve_time, digits = 4), " s")
    print(io, "\n  n_sim_steps: ", size(sol.X, 2) - 1)
    print(io, "\n  n_players: ", np)
    print(io, "\n  costs: ", sol.costs')
    print(io, "\n  coupled_constraints_satisfied: ",
          get(sol.solver_info, :all_coupled_constraints_satisfied, true))
    for i in 1:min(np, 3)
        print(io, "\n  Agent $i: final s=", round(sol.s_trajectories[i][end], digits = 4),
              " cost=", round(sol.costs[i], digits = 4))
    end
    np > 3 && print(io, "\n  ... ($(np) agents)")
end

Base.show(io::IO, ::MIME"text/plain", x::DynamicTubeMPCProblem) = show(io, x)
Base.show(io::IO, ::MIME"text/plain", x::DynamicTubeMPCSolution) = show(io, x)

function Base.show(io::IO, solver::DynamicTubeMPC{S}) where {S}
    print(io, "DynamicTubeMPC{$S}:")
    print(io, "\n  inner_solver: ", typeof(solver.inner_solver))
    print(io, "\n  warm_start: ", solver.warm_start)
    print(io, "\n  params: ", solver.params === nothing ? "from problem" : "override")
end
Base.show(io::IO, ::MIME"text/plain", x::DynamicTubeMPC) = show(io, x)
