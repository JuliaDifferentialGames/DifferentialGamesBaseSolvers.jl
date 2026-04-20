# ============================================================================
# utils.jl  —  ALGAMES internal utilities
#
# ── Residual / y-vector layout ───────────────────────────────────────────────
#
#   Derived directly from the reference implementation (ALGAMES.jl,
#   src/core/newton_core.jl, vertical_indices / horizontal_indices).
#
#   G (residual) and y (primal-dual vector) share the same row layout:
#
#     For each player i ∈ 1:np:
#       For each step k ∈ 1:N:
#         x-block(i,k+1)   n  entries  ← stationarity ∂L^i/∂x(k+1)
#         u-block(i,k)     mᵢ entries  ← stationarity ∂L^i/∂uᵢ(k)
#     For each step k ∈ 1:N:
#       dyn-block(k)       n  entries  ← dynamics residual D_k
#
#   Note: each player has its OWN x-rows.  Player 1's x(k+1) row and
#   player 2's x(k+1) row are different entries of G.  They both contain
#   the same dynamics Jacobian costate term but different cost gradients.
#
#   Total length of y:
#     np × N × (n + mᵢ) + N × n
#   where mᵢ varies per player.
#
# ── Stationarity per (i, k) ──────────────────────────────────────────────────
#
#   Reading from global_quantities.jl (residual! function):
#
#   x(k+1) row for player i:
#     G^i_{x(k+1)} = ∂J^i/∂x(k+1)                 [cost gradient at x(k+1)]
#                  + (∂f_k/∂x_k)ᵀ λᵢₖ              [Aₖᵀ λᵢₖ  — costate backprop]
#                  - λᵢₖ                             [−I term from ∂D_k/∂x(k+1)]
#                  + constraint gradient term
#
#   Wait — let me be precise by reading the stamps:
#
#     add2sub(res_sub[(:opt,i,:x,1,k)],   ∇dyn[:,∂x]' * λᵢₖ)   ← Aₖᵀ λᵢₖ at x(k)
#     add2sub(res_sub[(:opt,i,:u,i,k)],   ∇dyn[:,∂uᵢ]' * λᵢₖ)  ← Bᵢₖᵀ λᵢₖ at uᵢ(k)
#     add2sub(res_sub[(:opt,i,:x,1,k+1)], -λᵢₖ)                 ← −λᵢₖ at x(k+1)
#
#   So the structure at step k is:
#
#     G^i_{x(k)}   += ∂J^i/∂x(k)  +  Aₖᵀ λᵢₖ        (cost + costate forward)
#     G^i_{u^i(k)} += ∂J^i/∂u^i(k) + Bᵢₖᵀ λᵢₖ       (cost + costate)
#     G^i_{x(k+1)} += −λᵢₖ                            (costate backward link)
#
#   Plus terminal cost: ∂J^i/∂x(N+1) is added to G^i_{x(N+1)} (the k=N
#   x-block). In the reference, the cost gradient loop runs for k=1:N covering
#   all state knot points including the terminal.
#
#   Dynamics residual block:
#     G_{dyn,k} = f(x_k, u_k) − x(k+1)    for k = 1:N
#
# ── Constraint term ──────────────────────────────────────────────────────────
#
#   Added by constraint_residual! after the per-player loop. For each player i:
#     G^i_{x(k)}   += ∂C/∂x(k)ᵀ · pen_k
#     G^i_{u^i(k)} += ∂C/∂u^i(k)ᵀ · pen_k
#   where pen_k = Iᵨ·(λ + ρ·C(x_k, u_k)).
#
# ── y-vector ─────────────────────────────────────────────────────────────────
#
#   y mirrors G:  y contains the decision variables that correspond to each G row.
#
#     For player i, step k:
#       x-block(i,k+1): x(k+1)      ← shared state, but replicated per player
#       u-block(i,k):   u^i(k)      ← player i's controls
#     For step k:
#       dyn-block(k):   λᵢₖ for i=1 (dynamics multiplier — one set per player,
#                       but in the reference they are per-player in the column index)
#
#   Actually: the reference horizontal_indices assigns x(k+1) as a SHARED
#   column (not per-player), and λᵢₖ as a per-player column. This means H
#   is rectangular when viewed as (G rows) × (y columns).  The row (G) layout
#   is per-player for the stationarity conditions; the column (y) layout is
#   shared for x and per-player for u and λ.
#
#   For our dense implementation we use the simplest possible layout that
#   is correct and easy to verify:
#
#   y = [x(2), …, x(N+1),              ← N × n  (shared state, not per-player)
#        U[1](:,1), …, U[1](:,N),      ← m₁·N
#        …,
#        U[np](:,1), …, U[np](:,N),    ← mₙₚ·N
#        Λ_dyn[1](:,1), …, Λ_dyn[1](:,N),  ← n·N per player
#        …,
#        Λ_dyn[np](:,1), …, Λ_dyn[np](:,N)]
#
#   G = [G^1_{x(2)}, …, G^1_{x(N+1)},    ← player 1's x-rows
#        G^1_{u^1(1)}, …, G^1_{u^1(N)},   ← player 1's u-rows
#        G^2_{x(2)}, …, G^2_{x(N+1)},    ← player 2's x-rows
#        G^2_{u^2(1)}, …, G^2_{u^2(N)},   ← player 2's u-rows
#        …,
#        G_{dyn,1}, …, G_{dyn,N}]         ← dynamics rows (shared)
#
#   This layout is NOT the same as y.  G and y have the same length but
#   different ordering.  H = ∂G/∂y is a square matrix in this layout.
#   Since we compute H via ForwardDiff.jacobian(y -> G(y)), the layout
#   is automatically consistent.
#
# ============================================================================

# ── Namespace note ────────────────────────────────────────────────────────────
# Include order in DifferentialGamesBaseSolvers.jl:
#   include("solvers/ALGAMES/src/algames.jl")   ← defines ALGAMESWorkspace, ALGAMES
#   include("solvers/ALGAMES/src/utils.jl")     ← this file; uses those types
#
# algames.jl MUST be included first. Method signatures here reference
# ALGAMESWorkspace and ALGAMES; Julia evaluates signatures at definition time,
# so both types must exist when this file is parsed.
#
# All other names used here (LinearAlgebra, ForwardDiff, DifferentialGamesBase
# types) are provided by the enclosing module.
# Do not add `using` or `import` statements to this file.

# ============================================================================
# ALGAMESBuffers
# ============================================================================

"""
    ALGAMESBuffers{T}

Preallocated buffers for one Newton iteration.
`G` (residual) and `H` (Jacobian) are overwritten each call.
"""
mutable struct ALGAMESBuffers{T}
    G::Vector{T}
    H::Matrix{T}
end

ALGAMESBuffers(::Type{T}, n_y::Int) where {T} =
    ALGAMESBuffers{T}(zeros(T, n_y), zeros(T, n_y, n_y))

# ============================================================================
# y-vector dimensions and index helpers
# ============================================================================

"""
    _y_dim(ws) -> Int

Total length of the primal-dual vector y (= length of G).

Layout:
  n·N           (states x(2)…x(N+1), shared — not replicated per player)
  Σᵢ mᵢ · N    (per-player controls)
  np · n · N    (per-player dynamics multipliers Λ_dyn)
"""
function _y_dim(ws::ALGAMESWorkspace)
    return ws.n * ws.N +
           sum(ws.control_dims) * ws.N +
           ws.np * ws.n * ws.N
end

# ── y-vector index helpers ────────────────────────────────────────────────────
# These map (variable, player, step) to a range in y.
# G uses a different layout (see _G_dim and G-index helpers below).

"""Index of x(k+1) in y.  k ∈ 1:N → x(2)…x(N+1)."""
@inline function _x_idx(ws::ALGAMESWorkspace, k::Int)
    return (k-1)*ws.n + 1 : k*ws.n
end

"""Index of U[i][:,k] in y."""
@inline function _u_idx(ws::ALGAMESWorkspace, i::Int, k::Int)
    base = ws.n * ws.N + sum(ws.control_dims[1:i-1]) * ws.N + (k-1)*ws.control_dims[i]
    return base + 1 : base + ws.control_dims[i]
end

"""Index of Λ_dyn[i][:,k] in y."""
@inline function _λdyn_idx(ws::ALGAMESWorkspace, i::Int, k::Int)
    base = ws.n * ws.N +
           sum(ws.control_dims) * ws.N +
           (i-1) * ws.n * ws.N +
           (k-1) * ws.n
    return base + 1 : base + ws.n
end

# ── G-vector index helpers ────────────────────────────────────────────────────
# G layout (matches the description in the file header):
#   For i ∈ 1:np, k ∈ 1:N:
#     G^i_x(k+1)   n entries
#     G^i_u^i(k)   mᵢ entries
#   For k ∈ 1:N:
#     G_dyn(k)     n entries

"""
    _G_dim(ws) -> Int

Total length of the residual vector G.
For unconstrained problems = _y_dim(ws).
"""
function _G_dim(ws::ALGAMESWorkspace)
    stat = sum(ws.np * (ws.n + ws.control_dims[i]) for i in 1:ws.np) * ws.N
    # Simplified: same as _y_dim since G and y have the same total size.
    # stat = np * n * N (x-rows) + Σmᵢ * N (u-rows) + n * N (dyn-rows)
    return ws.np * ws.n * ws.N +
           sum(ws.control_dims) * ws.N +
           ws.n * ws.N
end

"""
    _Gx_idx(ws, i, k) -> UnitRange

Row range in G for player i's x(k+1) stationarity block.
k ∈ 1:N  →  covers x(2) through x(N+1).

G layout: for player i, each step k has (n + mᵢ) entries: x-block first, then u-block.
Player i's block starts after all blocks of players 1…i-1.
"""
@inline function _Gx_idx(ws::ALGAMESWorkspace, i::Int, k::Int)
    # Bytes before player i's block
    base_i = sum((ws.n + ws.control_dims[ii]) for ii in 1:i-1; init=0) * ws.N
    # Within player i: step k's x-block
    stride = ws.n + ws.control_dims[i]
    base   = base_i + (k-1) * stride
    return base + 1 : base + ws.n
end

"""
    _Gu_idx(ws, i, k) -> UnitRange

Row range in G for player i's u^i(k) stationarity block.
"""
@inline function _Gu_idx(ws::ALGAMESWorkspace, i::Int, k::Int)
    base_i = sum((ws.n + ws.control_dims[ii]) for ii in 1:i-1; init=0) * ws.N
    stride = ws.n + ws.control_dims[i]
    base   = base_i + (k-1) * stride + ws.n   # skip x-block
    return base + 1 : base + ws.control_dims[i]
end

"""
    _Gdyn_idx(ws, k) -> UnitRange

Row range in G for the dynamics residual at step k.
"""
@inline function _Gdyn_idx(ws::ALGAMESWorkspace, k::Int)
    stat_total = sum(ws.n + ws.control_dims[i] for i in 1:ws.np; init=0) * ws.N
    base       = stat_total + (k-1) * ws.n
    return base + 1 : base + ws.n
end

# ============================================================================
# _pack_y / _unpack_y
# ============================================================================

"""
    _pack_y(ws) -> Vector{T}

Pack the current iterate into y (primal-dual decision variables).
x₀ = x(1) is NOT in y (it is the fixed initial condition).
"""
function _pack_y(ws::ALGAMESWorkspace{T}) where {T}
    n_y = _y_dim(ws)
    y   = Vector{T}(undef, n_y)
    for k in 1:ws.N;          y[_x_idx(ws, k)]        .= ws.X[:, k+1];      end
    for i in 1:ws.np, k in 1:ws.N; y[_u_idx(ws, i, k)]  .= ws.U[i][:, k];   end
    for i in 1:ws.np, k in 1:ws.N; y[_λdyn_idx(ws,i,k)] .= ws.Λ_dyn[i][:,k]; end
    return y
end

"""
    _unpack_y(ws_ref, y) -> ALGAMESWorkspace{S}

Unpack y into a workspace.  Fixed quantities (x₀, λ, ρ, dims) come from ws_ref.
S may differ from T (e.g. ForwardDiff dual numbers during AD).
"""
function _unpack_y(ws_ref::ALGAMESWorkspace{T}, y::AbstractVector{S}) where {T, S}
    n, np, N = ws_ref.n, ws_ref.np, ws_ref.N
    cdim     = ws_ref.control_dims
    coff     = ws_ref.control_offsets

    X_new       = Matrix{S}(undef, n, N+1)
    X_new[:, 1] .= ws_ref.X[:, 1]    # x₀ fixed
    for k in 1:N
        X_new[:, k+1] .= y[_x_idx(ws_ref, k)]
    end

    U_new = [Matrix{S}(undef, cdim[i], N) for i in 1:np]
    for i in 1:np, k in 1:N
        U_new[i][:, k] .= y[_u_idx(ws_ref, i, k)]
    end

    Λ_new = [Matrix{S}(undef, n, N) for _ in 1:np]
    for i in 1:np, k in 1:N
        Λ_new[i][:, k] .= y[_λdyn_idx(ws_ref, i, k)]
    end

    # λ and ρ are outer-loop variables — not differentiated w.r.t.
    λ_new = convert(Vector{S}, ws_ref.λ)
    ρ_new = convert(Vector{S}, ws_ref.ρ)

    return ALGAMESWorkspace{S}(
        X_new, U_new, Λ_new, λ_new, ρ_new,
        n, np, N, ws_ref.nc_step, cdim, coff
    )
end

# ============================================================================
# _build_residual!
# ============================================================================

"""
    _build_residual!(buf, game, ws)

Assemble G into buf.G by calling _G_into!.
"""
function _build_residual!(buf::ALGAMESBuffers{T},
                           game::GameProblem{T},
                           ws  ::ALGAMESWorkspace{T}) where {T}
    fill!(buf.G, zero(T))
    _G_into!(buf.G, game, ws)
end

"""
    _G_into!(G, game, ws)

Core residual assembly.  G and ws may use a different element type S
(e.g. ForwardDiff duals) from the problem type T.

Structure (see file header for full derivation):

  For each step k ∈ 1:N and player i ∈ 1:np:
    G^i_{x(k+1)} = ∂J^i/∂x(k+1)                       [cost gradient at x(k+1)]
                 + Aₖᵀ λᵢₖ                              [costate: ∂D_k/∂x(k) * λᵢₖ]
                 − λᵢₖ                                  [∂D_{k}/∂x(k+1) = −I, × λᵢₖ]
                 + Jcxᵀ pen_k                            [constraint AL gradient]

    G^i_{u^i(k)} = ∂J^i/∂u^i(k)                        [cost gradient]
                 + Bᵢₖᵀ λᵢₖ                              [costate: ∂D_k/∂u^i(k) * λᵢₖ]
                 + Jcuᵢᵀ pen_k                           [constraint AL gradient]

  For each step k ∈ 1:N:
    G_{dyn,k} = f(x_k, u_k) − x(k+1)                   [dynamics feasibility]

TRANSPOSE AUDIT (most common silent bug):
  Every occurrence of Aₖ, Bᵢₖ in the residual uses the TRANSPOSE:  Aₖᵀ, Bᵢₖᵀ.
  This is the PMP costate update ∂D/∂x · λ = Aₖᵀ λ (NOT Aₖ λ).
  The Jacobian's off-diagonal blocks are also transposed accordingly.

Terminal cost:
  At k = N the term ∂J^i/∂x(k+1) = ∂J^i/∂x(N+1) is the TERMINAL cost gradient.
  The stage cost contributes to x(k) = x(N), not x(N+1).
  Implementation: for k = 1…N, cost_gradient at x(k) → G^i_{x(k)};
                  terminal_cost_gradient at x(N+1) → G^i_{x(N+1)}.

Boundary condition at k=1:
  The term ∂J^i/∂x(k) at k=1 would write to G^i_{x(1)} = G^i_{x₀},
  but x₀ is fixed, so this row does NOT exist in G.  We skip k=1 for
  stage-cost-at-xₖ contributions; the x(k+1) = x(2) row gets the
  Aₖᵀ λᵢ₁ and −λᵢ₁ terms from step k=1.
"""
function _G_into!(G ::AbstractVector{S},
                  game::GameProblem{T},
                  ws  ::ALGAMESWorkspace{S}) where {T, S}
    n, np, N = ws.n, ws.np, ws.N
    cdim     = ws.control_dims
    coff     = ws.control_offsets
    U_flat   = _concat_controls(ws)    # (m_total × N)

    for k in 1:N
        x_k  = ws.X[:, k]
        x_k1 = ws.X[:, k+1]
        u_k  = U_flat[:, k]

        # ── Discrete dynamics step and residual ──────────────────────────
        # D_k = x_next(x_k, u_k) − x(k+1)
        # For LinearDynamics: x_next = A x + B u  (exact discrete map).
        # For continuous dynamics: x_next = RK4(f, x_k, u_k, dt).
        # Using the same _discrete_step in the Jacobian (via _dyn_jac_discrete)
        # ensures the residual and its derivative are consistent.
        dt   = S(game.time_horizon.dt)
        f_k  = _discrete_step(game.dynamics, x_k, u_k, nothing, k, dt)
        D_k  = f_k .- x_k1

        # ── Discrete dynamics Jacobian ∂x_next/∂x and ∂x_next/∂u ─────────
        # TRANSPOSE AUDIT: Ak = ∂x_next/∂x, Bk = ∂x_next/∂u.
        # Only Aₖᵀ and Bᵢₖᵀ appear in the residual — never plain Aₖ.
        Ak, Bk = _dyn_jac_discrete(game.dynamics, x_k, u_k, n, k, dt)

        # ── Constraint evaluation and AL penalty ──────────────────────────
        C_k, Jcx_k, Jcu_k = _con_eval(game, ws, x_k, u_k, k)
        pen_k = _al_pen(C_k, ws, k, game)

        # ── Per-player stationarity ───────────────────────────────────────
        for i in 1:np
            λᵢₖ     = ws.Λ_dyn[i][:, k]      # dynamics multiplier μᵛ_k : (n,)
            mi       = cdim[i]
            u_rng_i  = coff[i]+1 : coff[i]+mi
            obj_i    = get_objective(game, i)

            # Cost gradients w.r.t. full joint (x_k, u_k).
            #
            # Two cost conventions coexist in DGB:
            #
            # 1. NonlinearStageCost (from minimize(CompositeCostTerm...)):
            #    gradient = (x,u,p,t) -> cost_term_gradient(stage, x, u, p, t)
            #    cost_term_gradient concatenates z=[x;u] and differentiates through
            #    evaluate_cost_term, which calls player_slice(x, offset, dim) and
            #    player_slice(u, offset, dim) internally.
            #    → Must receive full joint (x_k, u_k); returns (∇x of length n, ∇u of length m_total).
            #
            # 2. LQStageCost (from LQGameProblem / direct construction):
            #    stage_cost_gradient does Q*x + M*u + q with no slicing.
            #    Q is sized (n×n) for shared-state games, (nᵢ×nᵢ) for separable games.
            #    → For shared-state games: pass full (x_k, u_k[u_rng_i]).
            #    → For separable games: pass private slice (x_k[s_rng_i], u_k[u_rng_i]).
            #
            # We dispatch on stage cost type to handle both cases correctly.
            gx_ik, gu_ik_full = if obj_i.stage_cost isa LQStageCost
                # LQStageCost: no offset fields; must receive correctly-sized inputs.
                if length(game.metadata.state_dims) == np
                    # Separable: pass private state slice and private control
                    sdim_i = game.metadata.state_dims[i]
                    soff_i = game.metadata.state_offsets[i]
                    s_rng_i = soff_i+1 : soff_i+sdim_i
                    gx_local, gu_local = stage_cost_gradient(
                        obj_i.stage_cost, x_k[s_rng_i], u_k[u_rng_i], nothing, k
                    )
                    # Embed private gradient back into full-n vector
                    gx_full = zeros(S, n)
                    gx_full[s_rng_i] .= gx_local
                    gu_full = zeros(S, sum(cdim))
                    gu_full[u_rng_i] .= gu_local
                    gx_full, gu_full
                else
                    # Shared state: Q is (n×n), pass full x and private u
                    gx_full, gu_local = stage_cost_gradient(
                        obj_i.stage_cost, x_k, u_k[u_rng_i], nothing, k
                    )
                    gu_full = zeros(S, sum(cdim))
                    gu_full[u_rng_i] .= gu_local
                    gx_full, gu_full
                end
            else
                # NonlinearStageCost / CompositeCostTerm: pass full joint (x, u).
                # Returns (∇x of length n, ∇u of length m_total).
                stage_cost_gradient(obj_i.stage_cost, x_k, u_k, nothing, k)
            end
            # Player i's control gradient is the u_rng_i slice of gu_ik_full
            gu_ik = gu_ik_full[u_rng_i]

            # ── G^i_{u^i(k)} row ────────────────────────────────────────────
            # = ∂J^i/∂u^i(k) + Bᵢₖᵀ λᵢₖ + Jcuᵢᵀ pen_k
            # TRANSPOSE: Bᵢₖ is (n×mᵢ), so Bᵢₖᵀ is (mᵢ×n), times λᵢₖ (n,) → (mᵢ,) ✓
            Bᵢₖ  = Bk[:, u_rng_i]       # (n × mᵢ)
            Jcuᵢ = Jcu_k[:, u_rng_i]    # (nc_step × mᵢ), zero-row if nc_step=0
            @views G[_Gu_idx(ws, i, k)] .=
                gu_ik .+ Bᵢₖ' * λᵢₖ .+ Jcuᵢ' * pen_k

            # ── G^i_{x(k)} row  (skipped for k=1 since x₁=x₀ is fixed) ────
            # = ∂J^i/∂x(k) + Aₖᵀ λᵢₖ + Jcxᵀ pen_k
            # TRANSPOSE: Ak is (n×n), Ak' times λᵢₖ (n,) → (n,) ✓
            if k >= 2
                xk_row = _Gx_idx(ws, i, k-1)   # x(k) lives at index k-1 in G
                @views G[xk_row] .+= gx_ik .+ Ak' * λᵢₖ .+ Jcx_k' * pen_k
            end

            # ── G^i_{x(k+1)} row: accumulate −λᵢₖ ──────────────────────────
            # This is ∂/∂x(k+1)[(μᵛ)ᵀ D_k] = −μᵛₖ = −λᵢₖ
            @views G[_Gx_idx(ws, i, k)] .+= -λᵢₖ

            # ── Terminal cost gradient at x(N+1) (k=N only) ─────────────────
            if k == N
                gx_term = if obj_i.terminal_cost isa LQTerminalCost
                    if length(game.metadata.state_dims) == np
                        # Separable LQTerminalCost: Qf is (nᵢ×nᵢ), needs private slice
                        sdim_i = game.metadata.state_dims[i]
                        soff_i = game.metadata.state_offsets[i]
                        s_rng_i = soff_i+1 : soff_i+sdim_i
                        gx_local = terminal_cost_gradient(obj_i.terminal_cost, x_k1[s_rng_i], nothing)
                        gx_full  = zeros(S, n)
                        gx_full[s_rng_i] .= gx_local
                        gx_full
                    else
                        # Shared LQTerminalCost: Qf is (n×n), pass full x
                        terminal_cost_gradient(obj_i.terminal_cost, x_k1, nothing)
                    end
                else
                    # NonlinearTerminalCost / CompositeTerminalCostTerm:
                    # cost_term_gradient uses player_slice internally → full x
                    terminal_cost_gradient(obj_i.terminal_cost, x_k1, nothing)
                end
                @views G[_Gx_idx(ws, i, N)] .+= gx_term
            end
        end

        # ── Dynamics residual row (one per step, not per player) ──────────
        @views G[_Gdyn_idx(ws, k)] .= D_k
    end
end

# ============================================================================
# _build_jacobian!  via ForwardDiff
# ============================================================================

"""
    _build_jacobian!(buf, game, ws)

Compute H = ∂G/∂y via ForwardDiff.jacobian through _G_from_y.

Note: G and y have different orderings (see file header), so H is NOT block-
diagonal. ForwardDiff handles this correctly.
"""
function _build_jacobian!(buf::ALGAMESBuffers{T},
                           game::GameProblem{T},
                           ws  ::ALGAMESWorkspace{T}) where {T}
    y0 = _pack_y(ws)
    buf.H .= ForwardDiff.jacobian(y -> _G_from_y(game, ws, y), y0)
end

"""
    _G_from_y(game, ws_ref, y) -> Vector{S}

Evaluate G at packed y (may contain ForwardDiff dual numbers).
Used by _build_jacobian! and _line_search.
"""
function _G_from_y(game  ::GameProblem{T},
                   ws_ref::ALGAMESWorkspace{T},
                   y     ::AbstractVector{S}) where {T, S}
    ws_tmp = _unpack_y(ws_ref, y)
    G_tmp  = zeros(S, _G_dim(ws_ref))
    _G_into!(G_tmp, game, ws_tmp)
    return G_tmp
end

# ============================================================================
# Newton step and regularisation
# ============================================================================

"""
    _newton_step(H, G) -> δy

Solve H·δy = −G via dense LU.  H is already regularised.
Falls back to −G (steepest descent) if LU fails.
"""
function _newton_step(H::Matrix{T}, G::Vector{T}) where {T}
    F = lu(H; check = false)
    issuccess(F) && return -(F \ G)
    return -G
end

"""
    _regularize!(H, ws, reg)

Add reg·I to the PRIMAL rows of H only (state x and control u blocks).
The dynamics-multiplier rows (Λ_dyn) are NOT regularised.

This matches the reference Algames.jl regularize_residual_jacobian! which
adds reg only to the (x,x) and (u,u) diagonal blocks — not to the (λ,λ)
block. Regularising the dynamics-multiplier rows corrupts the dynamics
residual, which is already linear in the multipliers.
"""
function _regularize!(H::Matrix{T}, ws::ALGAMESWorkspace, reg::Float64) where {T}
    r = T(reg)
    # Primal rows: first n*N (state) + sum(cdim)*N (control) entries of y
    n_primal = (ws.n + sum(ws.control_dims)) * ws.N
    @inbounds for i in 1:n_primal
        H[i, i] += r
    end
end

# ============================================================================
# Line search  (Algorithm 1)
# ============================================================================

"""
    _line_search(game, ws, G_norm_1, δy, solver) -> α

Armijo backtracking: find α ∈ (0,1] such that
  ‖G_trial(y + α·δy)‖₁ ≤ (1 − α·β) · ‖G(y)‖₁

Trial residuals are evaluated with ρ clamped to ρ_init (the fixed
`ρ_trial = 1.0` of the reference). This keeps the line search
well-conditioned as the AL penalty grows in the outer loop — matching
the reference Algames.jl line_search which uses prob.pen.ρ_trial (fixed)
rather than the current prob.pen.ρ (growing).
"""
function _line_search(
    game    ::GameProblem{T},
    ws      ::ALGAMESWorkspace{T},
    G_norm_1::T,
    δy      ::Vector{T},
    solver  ::ALGAMES
) where {T}
    β  = T(solver.ls_beta)
    τ  = T(solver.ls_tau)
    α  = one(T)

    # Build a copy of ws with ρ clamped to ρ_init for trial evaluations.
    # This matches the reference ρ_trial = 1.0 — the line search merit
    # function uses a fixed small penalty so α doesn't collapse as ρ grows.
    ρ_trial = fill(T(solver.ρ_init), length(ws.ρ))

    for _ in 1:solver.ls_iter
        y_trial = _pack_y(ws) .+ α .* δy
        G_trial = _G_from_y_with_rho(game, ws, y_trial, ρ_trial)
        if norm(G_trial, 1) ≤ (1 - α * β) * G_norm_1
            return α
        end
        α *= τ
    end
    return α
end

"""
    _G_from_y_with_rho(game, ws_ref, y, ρ_override) -> Vector{S}

Like _G_from_y but uses ρ_override instead of ws_ref.ρ.
Used by line search to evaluate trial residuals at fixed ρ_trial.
"""
function _G_from_y_with_rho(
    game     ::GameProblem{T},
    ws_ref   ::ALGAMESWorkspace{T},
    y        ::AbstractVector{T},
    ρ_trial  ::Vector{T}
) where {T}
    ws_tmp      = _unpack_y(ws_ref, y)
    ws_tmp.ρ   .= ρ_trial
    G_tmp       = zeros(T, _G_dim(ws_ref))
    _G_into!(G_tmp, game, ws_tmp)
    return G_tmp
end

# ============================================================================
# Convergence helpers
# ============================================================================

"""
    _stationarity_norm(buf, ws) -> T

Mean ‖·‖₁ of the stationarity (x and u) blocks of G.
"""
function _stationarity_norm(buf::ALGAMESBuffers{T}, ws::ALGAMESWorkspace{T}) where {T}
    n_stat = (ws.np * ws.n + sum(ws.control_dims)) * ws.N
    n_stat == 0 && return zero(T)
    return norm(buf.G[1:n_stat], 1) / n_stat
end

"""
    _dynamics_norm(buf, ws) -> T

Mean ‖·‖₁ of the dynamics residual blocks of G.
"""
function _dynamics_norm(buf::ALGAMESBuffers{T}, ws::ALGAMESWorkspace{T}) where {T}
    n_stat = (ws.np * ws.n + sum(ws.control_dims)) * ws.N
    n_dyn  = ws.n * ws.N
    n_dyn == 0 && return zero(T)
    return norm(buf.G[n_stat+1 : n_stat+n_dyn], 1) / n_dyn
end

"""
    _constraint_violation(game, ws) -> T

Maximum primal constraint violation across all steps and constraints.
max(max(C,0)) for inequalities, max(|C|) for equalities.
"""
function _constraint_violation(game::GameProblem{T},
                                ws  ::ALGAMESWorkspace{T}) where {T}
    ws.nc_step == 0 && return zero(T)
    U_flat = _concat_controls(ws)
    vmax   = zero(T)
    for k in 1:ws.N
        x_k = ws.X[:, k]
        u_k = U_flat[:, k]
        for c in Iterators.flatten((game.private_constraints,
                                    game.shared_constraints))
            cv = evaluate_constraint(c, x_k, u_k, nothing, k)
            v  = is_inequality(c) ? maximum(max.(cv, zero(T))) : maximum(abs.(cv))
            vmax = max(vmax, v)
        end
    end
    return vmax
end

# ============================================================================
# Dynamics Jacobian helper
# ============================================================================

"""
    _dyn_jac_discrete(dyn, x_k, u_k, n, k, dt) -> (Ak, Bk)

Returns ∂x_next/∂x (n×n) and ∂x_next/∂u (n×m_total) for the DISCRETE step.

For `LinearDynamics` (already discrete): exact matrices from `get_A`/`get_B_concatenated`.
For continuous dynamics (`SeparableDynamics`, `CoupledNonlinearDynamics`): differentiate
through `_rk4_step` via ForwardDiff — the same function used for the residual.

This is the TRANSPOSE AUDIT point: only `Ak'` and `Bᵢₖ'` appear in the
residual. This function returns the UN-transposed matrices. All callers transpose.
"""
function _dyn_jac_discrete(dyn::LinearDynamics{T}, x_k, u_k, n::Int, k::Int, dt) where {T}
    return get_A(dyn, k), get_B_concatenated(dyn, k)
end

function _dyn_jac_discrete(dyn::DynamicsSpec{T}, x_k, u_k, n::Int, k::Int, dt) where {T}
    J = ForwardDiff.jacobian(
        z -> _discrete_step(dyn, z[1:n], z[n+1:end], nothing, k, dt),
        vcat(x_k, u_k)
    )
    return J[:, 1:n], J[:, n+1:end]
end

# ============================================================================
# Constraint helpers
# ============================================================================

"""
    _con_eval(game, ws, x_k, u_k, k) -> (C_k, Jcx_k, Jcu_k)

Evaluate all constraints + Jacobians at (x_k, u_k, k).
Returns empty arrays of correct shape when nc_step == 0.
"""
function _con_eval(game::GameProblem{T},
                   ws  ::ALGAMESWorkspace{S},
                   x_k, u_k, k::Int) where {T, S}
    nc = ws.nc_step
    n  = ws.n
    m  = sum(ws.control_dims)

    if nc == 0
        return S[], zeros(S, 0, n), zeros(S, 0, m)
    end

    C_list   = Vector{S}[]
    Jcx_list = Matrix{S}[]
    Jcu_list = Matrix{S}[]

    for c in Iterators.flatten((game.private_constraints, game.shared_constraints))
        cv = evaluate_constraint(c, x_k, u_k, nothing, k)

        # Use ForwardDiff through evaluate_constraint for the Jacobian rather than
        # calling constraint_jacobian(c, ...) directly.  Some constraints (e.g.
        # ProximityConstraint) have analytical Jacobians that pre-allocate Float64
        # output arrays and are not ForwardDiff-compatible when x_k/u_k carry
        # dual numbers (as they do inside _build_jacobian!).  Differentiating
        # through evaluate_constraint is always correct and type-compatible.
        z = vcat(x_k, u_k)
        J = ForwardDiff.jacobian(
            z_var -> evaluate_constraint(c, z_var[1:n], z_var[n+1:end], nothing, k),
            z
        )
        Jx = J[:, 1:n]
        Ju = J[:, n+1:end]

        push!(C_list, cv);   push!(Jcx_list, Jx);   push!(Jcu_list, Ju)
    end

    return vcat(C_list...), vcat(Jcx_list...), vcat(Jcu_list...)
end

"""
    _al_pen(C_k, ws, k, game) -> pen_k

Compute active-set AL penalty Iᵨ·(λ + ρ·C) for step k (Eq. 5).

  Inactive inequality (C_j < 0 AND λ_j ≈ 0): pen_j = 0
  Otherwise:                                   pen_j = ρ_j·(λ_j + ρ_j·C_j)
"""
function _al_pen(C_k ::AbstractVector{S},
                 ws  ::ALGAMESWorkspace{S},
                 k   ::Int,
                 game::GameProblem{T}) where {T, S}
    ws.nc_step == 0 && return S[]

    λ_k = view(ws.λ, (k-1)*ws.nc_step+1 : k*ws.nc_step)
    ρ_k = view(ws.ρ, (k-1)*ws.nc_step+1 : k*ws.nc_step)

    pen = Vector{S}(undef, ws.nc_step)
    off = 0
    n   = ws.n
    m   = sum(ws.control_dims)

    for c in Iterators.flatten((game.private_constraints, game.shared_constraints))
        nc = _nc_con(c, n, m)
        for j in 1:nc
            g  = off + j
            λj = λ_k[g];  ρj = ρ_k[g];  Cj = C_k[g]
            if is_inequality(c) && Cj < zero(S) && abs(λj) < 1e-14
                pen[g] = zero(S)   # inactive: Eq. (5)
            else
                pen[g] = ρj * (λj + ρj * Cj)
            end
        end
        off += nc
    end
    return pen
end

"""
    _nc_con(c, n, m) -> Int

Output dimension of constraint c evaluated at a zero point.
Uses Float64 to avoid type issues during ForwardDiff passes.
"""
function _nc_con(c::AbstractConstraint, n::Int, m::Int)
    return length(evaluate_constraint(c, zeros(n), zeros(m), nothing, 1))
end

# ============================================================================
# Misc
# ============================================================================

"""
    _concat_controls(ws) -> Matrix  (m_total × N)

Stack all players' controls column-by-column.
"""
_concat_controls(ws::ALGAMESWorkspace) = vcat(ws.U...)

"""
    _rollout_flat(dyn, x0, U, N, dt) -> Matrix  (n × N+1)

Forward simulate `dyn` from `x0` under per-player controls `U`.
- `LinearDynamics`: `evaluate_dynamics` is already a discrete map; use directly.
- `SeparableDynamics`/`CoupledNonlinearDynamics`: `evaluate_dynamics` returns the
  continuous RHS ẋ; integrate with RK4 over timestep `dt`.
"""
function _rollout_flat(dyn::DynamicsSpec{T}, x0::Vector{T},
                        U::Vector{Matrix{T}}, N::Int, dt::T = T(0.1)) where {T}
    n = length(x0)
    X = Matrix{T}(undef, n, N+1)
    X[:, 1] .= x0
    U_flat   = vcat(U...)
    for k in 1:N
        X[:, k+1] .= _discrete_step(dyn, X[:, k], U_flat[:, k], nothing, k, dt)
    end
    return X
end

"""
    _discrete_step(dyn, x, u, p, k, dt) -> x_next

One discrete step of the dynamics, dispatching on type:
- `LinearDynamics`: exact discrete map `x_next = Ax + Bu`, `dt` unused.
- Continuous dynamics (`SeparableDynamics`, `CoupledNonlinearDynamics`):
  RK4 integration over timestep `dt` with zero-order hold on u.
  ForwardDiff-compatible: x and u may be dual numbers.
"""
_discrete_step(dyn::LinearDynamics, x, u, p, k::Int, dt) =
    evaluate_dynamics(dyn, x, u, p, k)

function _discrete_step(dyn::DynamicsSpec, x, u, p, k::Int, dt)
    t = (k - 1) * dt
    k1 = evaluate_dynamics(dyn, x,              u, p, t)
    k2 = evaluate_dynamics(dyn, x .+ dt/2 .* k1, u, p, t + dt/2)
    k3 = evaluate_dynamics(dyn, x .+ dt/2 .* k2, u, p, t + dt/2)
    k4 = evaluate_dynamics(dyn, x .+ dt    .* k3, u, p, t + dt)
    return x .+ (dt/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

"""
    _count_nc_step(game) -> Int

Number of constraint output dimensions per timestep (assumed constant).
"""
function _count_nc_step(game::GameProblem{T}) where {T}
    isempty(game.private_constraints) && isempty(game.shared_constraints) && return 0
    n = total_state_dim(game.dynamics)
    m = total_control_dim(game.dynamics)
    return sum(_nc_con(c, n, m)
               for c in Iterators.flatten((game.private_constraints,
                                           game.shared_constraints));
               init = 0)
end