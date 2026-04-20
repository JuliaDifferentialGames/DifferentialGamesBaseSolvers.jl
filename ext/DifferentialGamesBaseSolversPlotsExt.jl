# ============================================================================
# ext/DifferentialGamesBaseSolversPlotsExt.jl
#
# Plots.jl extension for DifferentialGamesBaseSolvers.
# Activated automatically when both packages are loaded:
#
#   using DifferentialGamesBaseSolvers
#   using Plots
#
# Provides:
#   plot(sol::GNEPSolution; kwargs...)
#       Two-panel figure: 2D trajectories (left) + convergence residuals (right).
#
#   plot_trajectories(sol; kwargs...)
#       Only the 2D trajectory panel. Useful for composing your own layouts.
#
#   plot_convergence(sol; kwargs...)
#       Only the convergence residuals panel.
#
# Both single-panel functions return a Plots.Plot and accept any Plots keyword
# argument, passed through to the underlying plot call.
#
# ── State layout assumption ──────────────────────────────────────────────────
# The 2D trajectory plot reads x-position and y-position from the joint state
# trajectory stored in sol.state_trajectory (n × N+1 matrix).
#
# For SeparableDynamics games (PDGNEProblem), each player's private state is
# a contiguous block: player i occupies rows [(i-1)*nᵢ+1 : i*nᵢ] and the
# first two components of that block are (x, y).
#
# For shared-state games (LQGameProblem), the x/y indices must be supplied
# explicitly via the `xy_indices` keyword argument.
#
# The plot uses game.metadata to detect the layout automatically when possible.
# ============================================================================

module DifferentialGamesBaseSolversPlotsExt

using DifferentialGamesBaseSolvers
using DifferentialGamesBase
using Plots

# Palette: one colour per player, cycling if needed
const _PLAYER_COLOURS = [:royalblue, :crimson, :forestgreen,
                          :darkorange, :purple, :teal]

_player_colour(i::Int) = _PLAYER_COLOURS[mod1(i, length(_PLAYER_COLOURS))]

# ============================================================================
# Public API — overload Plots.plot
# ============================================================================

"""
    plot(sol::GNEPSolution; layout=:both, xy_indices=nothing, kwargs...)

Plot an ALGAMES solution. Returns a `Plots.Plot`.

## Keyword arguments
- `layout`      : `:both` (default) — two-panel [trajectories | convergence];
                  `:trajectories` — trajectory panel only;
                  `:convergence`  — convergence panel only.
- `xy_indices`  : For shared-state games, a `Vector{Tuple{Int,Int}}` giving
                  the `(x_idx, y_idx)` joint-state indices for each player.
                  For PDGNEPs, detected automatically from metadata.
- Any remaining kwargs are forwarded to `Plots.plot`.
"""
function Plots.plot(sol::GNEPSolution; layout=:both, xy_indices=nothing, kwargs...)
    if layout == :trajectories
        return plot_trajectories(sol; xy_indices, kwargs...)
    elseif layout == :convergence
        return plot_convergence(sol; kwargs...)
    else
        p1 = plot_trajectories(sol; xy_indices)
        p2 = plot_convergence(sol)
        return plot(p1, p2; layout=(1,2), size=(900,400), kwargs...)
    end
end

# ============================================================================
# plot_trajectories
# ============================================================================

"""
    plot_trajectories(sol::GNEPSolution; xy_indices=nothing, kwargs...) -> Plot

2D trajectory plot. Each player's path is drawn as a coloured line with
scatter markers at the initial position (○) and final position (★).

For PDGNEPs (SeparableDynamics), x/y indices are auto-detected from
`game.metadata`: player i occupies `state_offsets[i]+1 : state_offsets[i]+state_dims[i]`
and positions are the first two components of that slice.

For shared-state games, supply `xy_indices::Vector{Tuple{Int,Int}}` where
entry `i` gives `(x_joint_idx, y_joint_idx)` for player i.
"""
function plot_trajectories(sol::GNEPSolution;
                            xy_indices = nothing,
                            kwargs...)
    game = sol.game
    np   = game.n_players
    X    = sol.state_trajectory   # (n × N+1)

    # Resolve position indices in joint state
    xi, yi = _resolve_xy_indices(game, np, xy_indices)

    p = plot(;
        aspect_ratio = :equal,
        xlabel       = "x",
        ylabel       = "y",
        title        = "Trajectories",
        legend       = :outertopright,
        kwargs...
    )

    for i in 1:np
        col  = _player_colour(i)
        xs   = X[xi[i], :]
        ys   = X[yi[i], :]

        # Path
        plot!(p, xs, ys;
            color     = col,
            linewidth = 2,
            label     = "Player $i"
        )
        # Start marker
        scatter!(p, [xs[1]], [ys[1]];
            color       = col,
            markershape = :circle,
            markersize  = 6,
            label       = false
        )
        # End marker
        scatter!(p, [xs[end]], [ys[end]];
            color       = col,
            markershape = :star5,
            markersize  = 8,
            label       = false
        )
    end

    return p
end

# ============================================================================
# plot_convergence
# ============================================================================

"""
    plot_convergence(sol::GNEPSolution; kwargs...) -> Plot

Semi-log plot of the per-outer-iteration convergence residuals stored in
`sol.solver_info[:history]`. Shows stationarity (opt), dynamics (dyn), and
constraint violation (con) as a function of outer AL iteration.

Returns a plain `Plots.Plot` if convergence history is available, or a text
annotation plot if it is not (e.g. solution was loaded from a file without
history).
"""
function plot_convergence(sol::GNEPSolution; kwargs...)
    info = sol.solver_info
    if !haskey(info, :history)
        return plot(; title="Convergence (no history stored)",
                      annotations=(0.5, 0.5, "Run solve() to generate history"),
                      kwargs...)
    end

    hist = info[:history]
    opt  = Float64.(hist[:opt])
    dyn  = Float64.(hist[:dyn])
    con  = Float64.(hist[:con])
    iters = 1:length(opt)

    p = plot(;
        yscale  = :log10,
        xlabel  = "Outer iteration",
        ylabel  = "Residual (log₁₀ scale)",
        title   = "Convergence",
        legend  = :topright,
        kwargs...
    )

    plot!(p, iters, max.(opt, 1e-16); label="‖G_stat‖",  color=:royalblue,  linewidth=2)
    plot!(p, iters, max.(dyn, 1e-16); label="‖dyn‖",     color=:crimson,    linewidth=2)
    plot!(p, iters, max.(con, 1e-16); label="max(C)",     color=:forestgreen, linewidth=2)

    # Mark convergence with a vertical dashed line if solver converged
    if sol.converged && sol.iterations < length(opt)
        vline!(p, [sol.iterations];
            linestyle = :dash,
            color     = :gray,
            label     = "converged"
        )
    end

    return p
end

# ============================================================================
# Internals
# ============================================================================

"""
    _resolve_xy_indices(game, np, xy_indices) -> (xi, yi)

Returns two length-np vectors of joint-state indices for x and y positions.

- If `xy_indices` is provided as `Vector{Tuple{Int,Int}}`, use it directly.
- If the game uses `SeparableDynamics` (PDGNEProblem), detect from metadata:
  player i occupies `state_offsets[i]+1 : state_offsets[i]+state_dims[i]`
  and positions are the first two components (index 1 and 2 of private state).
- Otherwise fall back to `(1, 2)` for all players with a warning.
"""
function _resolve_xy_indices(game, np, xy_indices)
    if xy_indices !== nothing
        @assert length(xy_indices) == np "xy_indices must have one entry per player"
        xi = [xy_indices[i][1] for i in 1:np]
        yi = [xy_indices[i][2] for i in 1:np]
        return xi, yi
    end

    meta = game.metadata
    if length(meta.state_dims) == np
        # PDGNEProblem: contiguous private-state blocks, positions at +1 and +2
        xi = [meta.state_offsets[i] + 1 for i in 1:np]
        yi = [meta.state_offsets[i] + 2 for i in 1:np]
        return xi, yi
    end

    # Shared-state fallback: assume positions at 1 and 2 for all players.
    # The user should pass xy_indices for anything more specific.
    if np > 1
        @warn "plot_trajectories: shared-state game with np=$np players. " *
              "Assuming x=x[1], y=x[2] for all players. " *
              "Pass xy_indices=[(xi1,yi1), ...] to override."
    end
    xi = fill(1, np)
    yi = fill(2, np)
    return xi, yi
end

end # module