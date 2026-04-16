module DifferentialGamesBaseSolvers

using LinearAlgebra
using DifferentialGamesBase
using Printf

# Import everything needed for method extension and type use.
# GameSolver and solve are defined in DGB/solve.jl — import to extend.
import DifferentialGamesBase:
    GameSolver,
    solve,
    _solve,
    WarmstartData,
    solver_capabilities

# ============================================================================
# Solvers
# ============================================================================

# LQ Games
include("solvers/FNELQ/src/fnelq.jl")

# iLQGames
include("solvers/iLQGames/src/ilqgames.jl")

# ALGAMES
include("solvers/ALGAMES/src/algames.jl")
include("solvers/ALGAMES/src/utils.jl")

# ============================================================================
# Visualisation stubs (implementations live in the Plots weak-dep extension)
#
# Defining the functions here — as stubs that error with a helpful message —
# means the names are always exported from this module regardless of whether
# Plots.jl is loaded. The extension DifferentialGamesBaseSolversPlotsExt
# overloads these methods when `using Plots` triggers it.
# ============================================================================

"""
    animate_solution(sol::GNEPSolution; kwargs...) -> Animation

Animate an iLQGames solution. Requires Plots.jl to be loaded:

```julia
using Plots
using DifferentialGamesBaseSolvers
sol  = solve_figure8(verbose=true)
anim = animate_solution(sol)
gif(anim, "figure8.gif"; fps=20)
```

Dispatches on player count and state dimension:
- 2-player, 4-state  → figure-8 single unicycle animation
- 3-player, 12-state → three-vehicle intersection animation

See `save_animation` to write the result to disk.
"""
function animate_solution(sol::GNEPSolution; kwargs...)
    error(
        "animate_solution requires Plots.jl. Load it first:\n\n" *
        "    using Plots\n\n" *
        "then call animate_solution again."
    )
end

"""
    save_animation(anim, path; fps=20)

Save an animation to disk. Requires Plots.jl. Supported formats:
- `.gif`  — animated GIF (no external dependencies)
- `.mp4`  — MPEG-4 video (requires ffmpeg on PATH)
"""
function save_animation(anim, path::String; fps::Int=20)
    error(
        "save_animation requires Plots.jl. Load it first:\n\n" *
        "    using Plots\n\n" *
        "then call save_animation again."
    )
end

# ============================================================================
# Exports
# ============================================================================

export FNELQ
export iLQGames
export ALGAMES
export animate_solution
export save_animation

end # module