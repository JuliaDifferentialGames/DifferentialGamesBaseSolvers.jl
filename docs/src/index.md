# DifferentialGamesBaseSolvers.jl

**DifferentialGamesBaseSolvers.jl** provides numerical solvers for computing Nash equilibria in differential games. All solvers implement the `solve(game, solver; kwargs...)` interface defined in [DifferentialGamesBase.jl](https://github.com/JuliaDifferentialGames/DifferentialGamesBase.jl).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaDifferentialGames/DifferentialGamesBaseSolvers.jl")
```

Most users should install the top-level package instead:

```julia
Pkg.add(url="https://github.com/JuliaDifferentialGames/DifferentialGames.jl")
```

## Available Solvers

### FNELQ — Feedback Nash Equilibrium for LQ Games

Exact backward-Riccati solver for discrete-time, finite-horizon, linear-quadratic games. Returns time-varying feedback gain matrices `Kᵢ(t)` and per-player costs.

```julia
sol = solve(game, FNELQ())
K_gains = sol.solver_info[:feedback_gains]   # per-player gain matrices
```

**Requirements:** `LQGameProblem` or `LTVLQGameProblem`

### iLQGames — Iterative Linear-Quadratic Games

Approximate feedback Nash equilibria for nonlinear games by iterating:
1. Linearize dynamics and quadraticize costs around the current trajectory
2. Solve the resulting LQ game with FNELQ
3. Roll out the new strategy with line-search backtracking

```julia
sol = solve(game, iLQGames(; max_iter=100, ε_conv=1e-6, verbose=false))
```

**Requirements:** `GameProblem` with ForwardDiff-compatible dynamics and costs

**Reference:** Fridovich-Keil et al., "Efficient Iterative Linear-Quadratic Approximations for Nonlinear Multi-Player General-Sum Differential Games" (ICRA 2020)

### ALGAMES — Augmented Lagrangian Games

Handles constrained games using an augmented Lagrangian outer loop. Supports both private and shared inequality/equality constraints.

```julia
sol = solve(game, ALGAMES(; max_outer=20, ρ_init=1.0))
```

**Requirements:** `GameProblem` with constraints

**Reference:** Le Cleac'h et al., "ALGAMES: A Fast Augmented Lagrangian Solver for Constrained Dynamic Games" (RSS 2020)

## Visualization (Optional)

Load `Plots.jl` to enable the animation extension:

```julia
using Plots
using DifferentialGamesBaseSolvers

sol  = solve(game, iLQGames())
anim = animate_solution(sol)
gif(anim, "game.gif"; fps=20)
```

## API Reference

```@docs
FNELQ
iLQGames
ALGAMES
animate_solution
save_animation
```

## Adding New Solvers

Use the `ExampleSolver` template in `src/solvers/ExampleSolver/`. The minimum interface is:

```julia
struct MySolver <: GameSolver end

function solver_capabilities(::Type{MySolver})
    return [:feedback_nash, :unconstrained]
end

function DifferentialGamesBase._solve(
    game::GameProblem{T},
    solver::MySolver,
    warmstart::Union{Nothing, WarmstartData},
    verbose::Bool
) where {T}
    # ... your implementation ...
    return GNEPSolution(game, trajectories; converged=true)
end
```
