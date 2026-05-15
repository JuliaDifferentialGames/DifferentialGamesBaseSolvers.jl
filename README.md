# DifferentialGamesBaseSolvers.jl

[![CI](https://github.com/JuliaDifferentialGames/DifferentialGamesBaseSolvers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaDifferentialGames/DifferentialGamesBaseSolvers.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-DifferentialGames.jl-blue.svg)](https://JuliaDifferentialGames.github.io/DifferentialGames.jl/stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

Numerical solvers for the [DifferentialGames.jl](https://github.com/JuliaDifferentialGames/DifferentialGames.jl) ecosystem. Implements the `solve(game, solver)` interface defined in [DifferentialGamesBase.jl](https://github.com/JuliaDifferentialGames/DifferentialGamesBase.jl).

> ⚠️ **Work in progress** — API may change before v1.0.0.

## Available Solvers

| Solver | Type | Constraints | Notes |
|--------|------|-------------|-------|
| `FNELQ` | Feedback Nash | ✗ | Exact for LQ games; O(N·n³) |
| `iLQGames` | Feedback Nash | ✗ | Nonlinear games; iterative LQ approximation |
| `ALGAMES` | Open-loop Nash | ✓ | Augmented Lagrangian; private + shared constraints |

## Installation

Most users should install the umbrella package:

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaDifferentialGames/DifferentialGames.jl")
```

To use this package directly:

```julia
Pkg.add(url="https://github.com/JuliaDifferentialGames/DifferentialGamesBaseSolvers.jl")
```

## Quick Start

```julia
using DifferentialGamesBase, DifferentialGamesBaseSolvers, LinearAlgebra

# LQ game
n = 4
game = LQGameProblem(
    0.9 * I(n),
    [I(n)[:, 1:1], I(n)[:, 3:3]],
    [diagm(ones(n)), diagm(ones(n))],
    [fill(0.1, 1, 1), fill(0.1, 1, 1)],
    [diagm(ones(n)), diagm(ones(n))],
    ones(n), 2.0; dt=0.1
)

# Exact LQ feedback Nash equilibrium
sol = solve(game, FNELQ())
println("Player 1 cost: ", get_cost(sol, 1))
```

## Solver Details

### FNELQ
Solves discrete-time finite-horizon LQ games exactly via backward Riccati recursion. Use this as the inner loop for nonlinear games (it's what `iLQGames` calls internally).

### iLQGames
Iterative LQ approximation for nonlinear games. At each iteration:
1. Linearize dynamics and quadraticize costs around the current trajectory
2. Solve the resulting LQ game with FNELQ
3. Roll out the updated strategy with Armijo line search

Based on Fridovich-Keil et al., ICRA 2020.

### ALGAMES
Augmented Lagrangian solver for constrained games. Handles box constraints on state/control, proximity constraints, and general shared nonlinear constraints.

Based on Le Cleac'h et al., RSS 2020.

## Visualization (Optional)

```julia
using Plots
using DifferentialGamesBaseSolvers

sol  = solve(game, iLQGames())
anim = animate_solution(sol)
gif(anim, "game.gif"; fps=20)
```

## Adding New Solvers

See `src/solvers/ExampleSolver/` for a documented template. The interface is:

```julia
struct MySolver <: GameSolver end

solver_capabilities(::Type{MySolver}) = [:feedback_nash, :unconstrained]

function DifferentialGamesBase._solve(
    game::GameProblem{T}, solver::MySolver,
    warmstart, verbose::Bool
) where {T}
    # ... implementation ...
    return GNEPSolution(game, trajectories; converged=true)
end
```

## License

MIT License — see LICENSE file for details.

## Disclosure of Generative AI Usage

Generative AI (Claude Sonnet 4.5/4.6) was used in the creation of this library as a programming aid including guided code generation, assistance with performance optimization, and documentation. All code and documentation has been reviewed by the author(s) for accuracy.
