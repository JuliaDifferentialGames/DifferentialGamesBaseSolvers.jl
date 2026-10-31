# ============================================================================
# runtests.jl — Main test runner for DTRHN solver
#
# NOT wired into the package's `test/runtests.jl` / `Pkg.test()`: the
# theory-faithful integration tests require FALCON (methods/FALCON), which is
# an external, UNREGISTERED package (local dev path only) and therefore
# cannot be a declared dependency of this package (would break `Pkg.test()`'s
# isolated sandbox env for every other contributor). Run this suite manually
# via an environment that has both packages available, e.g. the shared
# ~/.julia/environments/v1.12 dev environment used throughout this project:
#
#   julia --project=@v1.12 -e '
#       using DifferentialGamesBase, DifferentialGamesBaseSolvers, FALCON
#       include("src/solvers/DTRHN/test/runtests.jl")'
# ============================================================================

println("Running DTRHN solver tests...")

# Core implementation tests (A1–B7, integration with FALCON)
include("dtmpc_tests.jl")

# Theory verification: empirical validation of the Theory Spine v0.2's observable claims
include("theory_validation.jl")

println("DTRHN tests completed.")