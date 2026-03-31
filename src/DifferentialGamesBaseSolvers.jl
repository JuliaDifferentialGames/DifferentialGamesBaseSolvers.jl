module DifferentialGamesBaseSolvers

# ============================================================================
# Dependencies
# ============================================================================

# Core DifferentialGamesBase types — all game problem, dynamics, cost, and
# solution types are re-exported from here; solvers never import them directly.
using DifferentialGamesBase

# Extend the solve/_solve dispatch from DifferentialGamesBase
import DifferentialGamesBase: solve, _solve

# AD and results — required by quadraticization.jl (DiffResults single-pass)
using ForwardDiff
using DiffResults

# Standard library
using LinearAlgebra
using SparseArrays

# ============================================================================
# FNELQ — Feedback Nash Equilibrium via LQ Dynamic Programming
# ============================================================================

include("solvers/FNELQ/src/fnelq.jl")

# ============================================================================
# PANGOLIN — Penalty-Augmented Nash Game via Operating-point Linearization
#            and Iterative Nash (iLQGames with AL constraints)
#
# Load order is strict — each file depends on types from the previous ones:
#   quadraticization  → OperatingPoint, LTVLQApproximation, lq_approximation!
#   al_augmentation   → ALOptions, ALAugmentedObjective, ALSolverState, dual update
#   PANGOLIN          → FeedbackStrategy, solve_lq_game!, rollout!, _solve
# ============================================================================

include("solvers/PANGOLIN/src/quadraticization.jl")
include("solvers/PANGOLIN/src/al_augmentation.jl")
include("solvers/PANGOLIN/src/PANGOLIN.jl")

# ============================================================================
# Exports
# ============================================================================

# ── FNELQ ────────────────────────────────────────────────────────────────────
export FNELQ

# ── PANGOLIN solver ──────────────────────────────────────────────────────────
export PANGOLIN
export PANGOLINCache

# ── Quadraticization layer (Phase 2) ─────────────────────────────────────────
export OperatingPoint
export LTVLQApproximation
export GeneralStageCost
export GeneralTerminalCost
export quadraticize
export linearize
export lq_approximation!

# ── AL augmentation layer (Phase 3) ──────────────────────────────────────────
export ALOptions
export ALAugmentedObjective
export ALSolverState
export augmented_stage_cost
export update_multipliers!
export update_shared_multipliers!
export maybe_update_penalty!
export constraint_violation
export reset_multipliers!
export reset_al_state!
export constraint_dim

# ── PANGOLIN internals (useful for testing and extension) ────────────────────
export FeedbackStrategy
export scale_feedforward!
export solve_lq_game!
export rollout!
export backtrack_scale!
export has_converged

# ── Solver protocol (re-exported for convenience) ────────────────────────────
export solve, _solve

end