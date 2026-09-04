module DifferentialGamesBaseSolvers

using LinearAlgebra
using DifferentialGamesBase
using ForwardDiff
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

# L-IBR
include("solvers/LIBR/src/libr.jl")

# Yi & Pavel — Distributed Primal-Dual GNE
include("solvers/YIPAVEL/src/yipavel.jl")

# Inga et al. — Inverse LQ Differential Games
include("solvers/INVLQ/src/invlq.jl")

# Receding Horizon Nash wrapper (wraps any inner GameSolver)
include("solvers/RHN/src/rhn.jl")

# Dynamic Tube MPC solver (extends RHN with robust tube optimization)
include("solvers/DTRHN/src/dtmpc.jl")

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
export LIBR
export YiPavel
export InverseLQGames
export RecedingHorizonNash
export DynamicTubeMPC
export DynamicTubeMPCProblem
export DynamicTubeMPCSolution
export DTMPCParameters
export TubeChannel
# DTMPC builders (augmented PD-GNEP)
export build_augmented_player_dynamics_discrete
export build_augmented_cost
export build_tightened_collision_constraint
export build_tightened_actuator_constraint
export build_tightened_rot_actuator_constraint
export build_tightened_detection_constraint
export build_initial_containment_constraint
export build_augmented_game
# DTMPC layout / accessors
export augmented_state_dims
export agent_physical_block
export extract_agent_augmented
export block_offsets
export n_tube_channels_agent
export tube_channels
export uses_multi_tube
export ρ_of
export tube_delta_hat
export dtmpc_free_initial_mask
export augmented_control_dims
export dtmpc_bandwidth_trajectory
# DTMPC utility functions
export lsh_boundary_layer_control
export lsh_rotational_boundary_layer_control
export tighten_coupled_distance
export check_coupled_constraint_satisfaction
export create_independent_uncertainty_models
export create_coupled_uncertainty_models
export default_drag_uncertainty
export default_coupled_drag_uncertainty
export validate_coupled_constraints
export get_coupled_constraint_violations
export all_coupled_constraints_satisfied
# DTMPC accessor functions
export dtmpc_state_trajectory
export dtmpc_control_trajectory
export dtmpc_tube_trajectory
export dtmpc_total_cost
export dtmpc_n_players
export dtmpc_state_dim
export dtmpc_agent_state_dim
export dtmpc_simulation_steps
export dtmpc_coupled_pairs
export animate_solution
export save_animation

end # module