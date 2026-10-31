# ============================================================================
# ext/DifferentialGamesBaseSolversFALCONExt.jl
#
# FALCON weak-dependency extension. Activated when FALCON is loaded alongside this package:
#
#   using DifferentialGamesBaseSolvers
#   using FALCON
#
# Gives the stub `falcon` its method. FALCONSolver is a regular
# DifferentialGamesBase `GameSolver`, so `solve(game, falcon(...))` needs nothing else.
# ============================================================================
module DifferentialGamesBaseSolversFALCONExt

using DifferentialGamesBaseSolvers
using FALCON

DifferentialGamesBaseSolvers.falcon(; kwargs...) = FALCON.FALCONSolver(; kwargs...)

end
