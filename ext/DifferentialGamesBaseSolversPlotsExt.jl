module DifferentialGamesBaseSolversPlotsExt

using DifferentialGamesBaseSolvers
using DifferentialGamesBase
using Plots
using Printf

# ============================================================================
# animate_solution — dispatch on player count / state structure
# ============================================================================

"""
    animate_solution(sol::GNEPSolution; kwargs...) -> Animation

Animate a `GNEPSolution` produced by iLQGames. Dispatches automatically:
- 2-player, 4-state game  → `animate_figure8`   (shared unicycle, iLQGames.jl style)
- 3-player, 12-state game → `animate_intersection` (3-vehicle intersection)

All keyword arguments are forwarded to the specific animation function.
"""
function DifferentialGamesBaseSolvers.animate_solution(
    sol::GNEPSolution;
    fps::Int             = 20,
    subsample::Int       = 1,
    trail_len::Int       = -1,   # -1 → full horizon for figure-8, 25 for intersection
    arrow_scale::Float64 = 0.3,
    figsize::Tuple       = (500, 500),
    kwargs...
)
    np = sol.game.n_players
    n  = size(sol.state_trajectory, 1)
    N  = size(sol.state_trajectory, 2) - 1

    if np == 2 && n == 4
        tl = trail_len < 0 ? N : trail_len
        return animate_figure8(sol; fps, subsample, trail_len=tl, figsize, kwargs...)
    elseif np == 3 && n == 12
        tl = trail_len < 0 ? 25 : trail_len
        return animate_intersection(sol; fps, subsample, trail_len=tl,
                                    arrow_scale, figsize, kwargs...)
    else
        error("animate_solution: no animation defined for np=$np, n=$n. " *
              "Call animate_figure8 or animate_intersection directly.")
    end
end

# ============================================================================
# Shared helpers
# ============================================================================

# Draw a heading arrow at (px, py) pointing in direction φ.
# Returns nothing (mutates the current plot via plot!).
function _draw_arrow!(px, py, φ, scale, color; lw=2)
    dx = scale * cos(φ)
    dy = scale * sin(φ)
    plot!([px, px + dx], [py, py + dy];
          arrow=true, color=color, lw=lw, label=false)
end

# Draw a filled circle (vehicle body) at (px, py).
function _draw_vehicle!(px, py, r, color)
    θs = range(0, 2π, length=20)
    xs = px .+ r .* cos.(θs)
    ys = py .+ r .* sin.(θs)
    plot!(xs, ys; seriestype=:shape, fillcolor=color, linecolor=color,
          alpha=0.7, label=false)
end

# Return the slice of indices into the trail that fall within [1, N+1].
function _trail_range(k, trail_len, N)
    lo = max(1, k - trail_len)
    return lo:k
end

# ============================================================================
# Figure-8 animation
# ============================================================================

"""
    animate_figure8(sol; kwargs...) -> Animation

Animate the two-player shared unicycle figure-8 game, matching the layout
of the `plot_traj` visualisation from iLQGames.jl:

  Left column  : u₁ (steering, red) and u₂ (acceleration, green) vs time step
  Right panel  : x-y trajectory with a growing trail and current-position dot

# Keyword arguments
- `fps::Int = 20`           : frames per second
- `subsample::Int = 2`      : animate every nth step
- `trail_len::Int = 200`    : timesteps of history shown on the xy panel
- `figsize::Tuple = (700,400)` : figure pixel dimensions
"""
function animate_figure8(
    sol::GNEPSolution;
    fps::Int             = 20,
    subsample::Int       = 2,
    trail_len::Int       = 200,
    figsize::Tuple       = (700, 400),
    kwargs...
)
    X   = sol.state_trajectory          # (4 × N+1)
    U1  = sol.trajectories[1].controls  # (1 × N) steering ω
    U2  = sol.trajectories[2].controls  # (1 × N) acceleration a
    N   = size(X, 2) - 1
    dt  = sol.game.time_horizon.dt
    ts  = 0:N-1   # time step indices for controls

    px_all = X[1, :]; py_all = X[2, :]
    u1_all = vec(U1);  u2_all = vec(U2)

    # Fixed axis limits across frames
    u1_lim = (minimum(u1_all) - 0.2, maximum(u1_all) + 0.2)
    u2_lim = (min(0.0, minimum(u2_all)) - 0.05, maximum(u2_all) + 0.05)
    pad    = 0.5
    xy_lim = (minimum(px_all) - pad, maximum(px_all) + pad)
    yy_lim = (minimum(py_all) - pad, maximum(py_all) + pad)

    frames = 1:subsample:(N+1)

    anim = @animate for k in frames
        lo = max(1, k - trail_len)
        tr = lo:k

        l = @layout [grid(2,1){0.38w}  a{0.62w}]
        plot(; layout=l, size=figsize, background_color=:white,
               left_margin=4Plots.mm, bottom_margin=4Plots.mm)

        # ── Top-left: u₁ (steering) ──────────────────────────────────────────
        # Show full control history in grey, highlight up to current step in red
        plot!(ts, u1_all;
              subplot=1, color=:lightgrey, lw=1, label=false,
              xlims=(0, N), ylims=u1_lim, grid=true,
              ylabel="u₁", xlabel="time step", title="")
        if k > 1
            plot!(ts[1:min(k-1,N)], u1_all[1:min(k-1,N)];
                  subplot=1, color=:red, lw=1.5, label=false)
        end
        vline!([k-1]; subplot=1, color=:black, lw=0.8, ls=:dash, label=false)

        # ── Bottom-left: u₂ (acceleration) ───────────────────────────────────
        plot!(ts, u2_all;
              subplot=2, color=:lightgrey, lw=1, label=false,
              xlims=(0, N), ylims=u2_lim, grid=true,
              ylabel="u₂", xlabel="time step", title="")
        if k > 1
            plot!(ts[1:min(k-1,N)], u2_all[1:min(k-1,N)];
                  subplot=2, color=:green, lw=1.5, label=false)
        end
        vline!([k-1]; subplot=2, color=:black, lw=0.8, ls=:dash, label=false)

        # ── Right: x-y trajectory ─────────────────────────────────────────────
        plot!(; subplot=3,
              xlims=xy_lim, ylims=yy_lim,
              aspect_ratio=:equal, grid=true,
              xlabel="pₓ [m]", ylabel="pᵧ [m]",
              title="", legend=false)

        # Full planned path (faint black, always visible for context)
        plot!(X[1, :], X[2, :];
              subplot=3, color=:black, lw=0.6, alpha=0.18, label=false)

        # Growing trail
        if length(tr) > 1
            plot!(X[1, tr], X[2, tr];
                  subplot=3, color=:black, lw=1.4, alpha=0.7, label=false)
        end

        # Current position — filled diamond (matching reference)
        scatter!([X[1, k]], [X[2, k]];
                 subplot=3, marker=:diamond, ms=7,
                 color=:black, markerstrokewidth=0, label=false)
    end

    return anim
end

# ============================================================================
# Three-player intersection animation
# ============================================================================

"""
    animate_intersection(sol; kwargs...) -> Animation

Animate the three-player unicycle intersection game.

State layout:
  x[1:4]  = (px₁, py₁, φ₁, v₁)  — Player 1 (south → north), red
  x[5:8]  = (px₂, py₂, φ₂, v₂)  — Player 2 (west  → east),  blue
  x[9:12] = (px₃, py₃, φ₃, v₃)  — Player 3 (north → south), green

The animation shows:
- Each vehicle as a filled circle with a heading arrow
- Trail of past positions per vehicle
- Goal position markers (★)
- Road layout (grey lane markings)
- Proximity radius circles when vehicles are close
"""
function animate_intersection(
    sol::GNEPSolution;
    fps::Int             = 20,
    subsample::Int       = 2,
    trail_len::Int       = 25,
    arrow_scale::Float64 = 0.3,
    figsize::Tuple       = (650, 650),
    vehicle_r::Float64   = 0.15,
    road_hw::Float64     = 0.8,      # half-width of road markings
    prox_r::Float64      = 1.5,      # proximity warning radius
    goals::Vector        = [[0.0, 4.0], [4.0, 0.0], [0.0, -4.0]],
    kwargs...
)
    X   = sol.state_trajectory       # (12 × N+1)
    N   = size(X, 2) - 1
    dt  = sol.game.time_horizon.dt

    # Per-player state accessors
    px(i, k) = X[4*(i-1)+1, k]
    py(i, k) = X[4*(i-1)+2, k]
    φ(i, k)  = X[4*(i-1)+3, k]
    v(i, k)  = X[4*(i-1)+4, k]

    colors = [:crimson, :royalblue, :forestgreen]
    names  = ["P1 (S→N)", "P2 (W→E)", "P3 (N→S)"]

    # Fixed axis limits around the intersection
    ax_lim = 5.5
    xlims  = (-ax_lim, ax_lim)
    ylims  = (-ax_lim, ax_lim)

    frames = 1:subsample:(N+1)

    anim = @animate for k in frames
        t_now = (k-1) * dt

        plot(; size=figsize, xlims=xlims, ylims=ylims,
               aspect_ratio=:equal, grid=false, legend=:topright,
               title=@sprintf("3-Player Intersection  t=%.1fs", t_now),
               xlabel="x (m)", ylabel="y (m)",
               background_color=:white)

        # ── Road markings ────────────────────────────────────────────────────
        # Vertical road (Players 1 & 3)
        plot!([-road_hw, -road_hw], [-ax_lim, ax_lim];
              color=:lightgrey, lw=1, label=false)
        plot!([ road_hw,  road_hw], [-ax_lim, ax_lim];
              color=:lightgrey, lw=1, label=false)
        # Horizontal road (Player 2)
        plot!([-ax_lim, ax_lim], [-road_hw, -road_hw];
              color=:lightgrey, lw=1, label=false)
        plot!([-ax_lim, ax_lim], [ road_hw,  road_hw];
              color=:lightgrey, lw=1, label=false)
        # Centre dashed lines
        plot!([0.0, 0.0], [-ax_lim, ax_lim];
              color=:grey, lw=0.8, ls=:dash, label=false)
        plot!([-ax_lim, ax_lim], [0.0, 0.0];
              color=:grey, lw=0.8, ls=:dash, label=false)

        # ── Goal markers ─────────────────────────────────────────────────────
        for i in 1:3
            scatter!([goals[i][1]], [goals[i][2]];
                     marker=:star5, ms=12, color=colors[i],
                     markerstrokewidth=0, alpha=0.5, label=false)
        end

        # ── Proximity warning circles (when vehicles are close) ───────────
        for i in 1:3, j in (i+1):3
            d = norm([px(i,k) - px(j,k), py(i,k) - py(j,k)])
            if d < prox_r
                # Draw a faint warning circle centred between the two vehicles
                mx = (px(i,k) + px(j,k)) / 2
                my = (py(i,k) + py(j,k)) / 2
                θs = range(0, 2π, length=40)
                r  = d / 2
                plot!(mx .+ r .* cos.(θs), my .+ r .* sin.(θs);
                      color=:orange, lw=1, alpha=0.4, label=false)
            end
        end

        # ── Per-vehicle trail, body, arrow ────────────────────────────────
        for i in 1:3
            c  = colors[i]
            tr = _trail_range(k, trail_len, N)

            # Trail
            if length(tr) > 1
                plot!([px(i, kk) for kk in tr],
                      [py(i, kk) for kk in tr];
                      color=c, lw=1.5, alpha=0.45, label=false)
            end

            # Vehicle body
            _draw_vehicle!(px(i, k), py(i, k), vehicle_r, c)

            # Heading arrow
            _draw_arrow!(px(i, k), py(i, k), φ(i, k), arrow_scale, c; lw=2)

            # Legend entry (only at first plotted position for cleanliness)
            scatter!([px(i, k)], [py(i, k)];
                     color=c, ms=0, label=@sprintf("%s  v=%.1f", names[i], v(i, k)))
        end
    end

    return anim
end

# ============================================================================
# save_animation — convenience wrapper
# ============================================================================

"""
    save_animation(anim::Animation, path::String; fps=20)

Save an `Animation` to `path`. Extension determines format:
- `.gif`  : animated GIF (default, broadly compatible)
- `.mp4`  : MPEG-4 video (requires ffmpeg on PATH)

# Example
```julia
using Plots
using DifferentialGamesBaseSolvers
sol  = solve_figure8()
anim = animate_solution(sol)
save_animation(anim, "figure8.gif")
```
"""
function DifferentialGamesBaseSolvers.save_animation(
    anim::Plots.Animation,
    path::String;
    fps::Int = 20
)
    if endswith(path, ".gif")
        gif(anim, path; fps=fps)
    elseif endswith(path, ".mp4")
        mp4(anim, path; fps=fps)
    else
        error("save_animation: unsupported extension in '$path'. Use .gif or .mp4.")
    end
    @info "Animation saved" path=path fps=fps frames=length(anim.frames)
    return path
end

end # module DifferentialGamesBaseSolversPlotsExt