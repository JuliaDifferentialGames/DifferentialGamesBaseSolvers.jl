module DifferentialGamesBaseSolvers

using LinearAlgebra
using DifferentialGamesBase
using ForwardDiff
using Printf
using Random
using SparseArrays

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

# Stochastic LQ Games (additive noise; certainty-equivalent to FNELQ)
include("solvers/SFNELQ/src/sfnelq.jl")

# Chance-constrained iterative LQ stochastic games (Zhong et al., RA-L 2023)
include("solvers/CCiLQGames/src/ccilqgames.jl")

# iLQGames
include("solvers/iLQGames/src/ilqgames.jl")

# ALGAMES
include("solvers/ALGAMES/src/algames.jl")
include("solvers/ALGAMES/src/utils.jl")

# DGSQP — Dynamic Game SQP (Zhu & Borrelli, ICRA 2023)
include("solvers/DGSQP/src/dgsqp.jl")

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
    animate_solution(sol::GNEPSolution; fps=20, subsample=1, xy_indices=nothing, kwargs...) -> Plots.Animation

Animate the 2D trajectories of any `GNEPSolution`: each player's path grows frame by frame
and its current position is marked. Requires Plots.jl to be loaded:

```julia
using Plots
anim = animate_solution(sol)
gif(anim, "game.gif"; fps=20)
```

Positions are read from the joint state: for `PDGNEProblem` solutions the first two components
of each player's state block are used automatically; for shared-state games pass
`xy_indices = [(x_idx, y_idx), …]`, one tuple per player. `subsample` keeps every
`subsample`-th time node.

See [`save_animation`](@ref) to write the result to disk.
"""
function animate_solution(sol; kwargs...)
    error(
        "animate_solution requires Plots.jl and a GNEPSolution. Load Plots first:\n\n" *
        "    using Plots\n\n" *
        "then call animate_solution again."
    )
end

"""
    save_animation(anim, path; fps=20)

Save a `Plots.Animation` to disk. Requires Plots.jl. Supported formats:
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
# FALCON stubs (implementations live in the FALCON weak-dep extension)
#
# FALCON (Outland and Arya) is a separate package that solves non-convex PD-GNEPs by sequential
# convexification and an augmented Lagrangian. It is a soft dependency: `using FALCON` activates
# DifferentialGamesBaseSolversFALCONExt, which gives these constructors their methods. The
# lowercase names avoid clashing with FALCON's own exports (`FALCONSolver`).
# ============================================================================

"""
    falcon(; kwargs...) -> FALCON.FALCONSolver

Construct a [FALCON](https://github.com/BennetOutland/FALCON.jl) solver for a `PDGNEProblem`. Requires
the FALCON package, which is a soft dependency:

```julia
using FALCON                       # activates the extension
solver = falcon(; ϵ_opt = 2e-2, ϵ_feas = 1e-4)
sol    = solve(game, solver; warmstart = ws)
```

Keyword arguments are those of `FALCON.FALCONSolver`.
"""
function falcon end

function __init__()
    # The stubs have no methods until FALCON is loaded; say so instead of a bare MethodError.
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        if exc.f === falcon
            print(io, "\nfalcon needs the FALCON package: run `using FALCON` first ",
                  "(install it with `Pkg.add(url = \"https://github.com/BennetOutland/FALCON.jl\")`).")
        end
    end
end

# ============================================================================
# Exports
# ============================================================================

export FNELQ
export SFNELQ
export CCiLQGames
export additive_noise_cost
export gaussian_trajectory_moments
export get_covariance_trajectory, get_expected_costs, get_mean_trajectory
export iLQGames
export ALGAMES
export DGSQP
export falcon
export dgsqp_paper_options
export dgsqp_racing_options
export dgsqp_tuned_options
export OpenLoopGame
export SeparableRolloutGame
export dgsqp_solve
export DGSQPResult
export best_response_gap
export control_ranges
export split_controls
export join_controls
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