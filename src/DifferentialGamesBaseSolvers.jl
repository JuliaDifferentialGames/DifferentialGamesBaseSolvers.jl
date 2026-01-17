module DifferentialGamesBaseSolvers

# Core DifferentialGamesBase types
using DifferentialGamesBase

# Import solve to extend it
import DifferentialGamesBase: solve
import DifferentialGamesBase: _solve

# Standard library
using LinearAlgebra
using SparseArrays

# Includes 
include("solvers/FNELQ/src/fnelq.jl")
include("solvers/FNELQ/src/utils.jl")

export FNELQ, _solve

end
