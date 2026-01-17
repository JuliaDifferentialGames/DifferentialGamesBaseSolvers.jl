"""
    extract_lq_matrices(game::GameProblem{T})

Extract LQ matrices from GameProblem structure.
Returns (A, B, Q, R, Qf, M, q, r) where vectors index by player.
"""
function extract_lq_matrices(game::GameProblem{T}) where {T}
    # Extract dynamics
    dyn = game.dynamics
    @assert dyn isa LinearDynamics "Requires linear dynamics"
    
    A = dyn.A
    B_vec = dyn.B
    
    n_players = game.n_players
    n = size(A, 1)
    
    # Extract cost matrices from objectives
    Q = Matrix{T}[]
    R = Matrix{T}[]
    Qf = Matrix{T}[]
    M = Matrix{T}[]
    q = Vector{T}[]
    r = Vector{T}[]
    
    for i in 1:n_players
        obj = game.objectives[i]
        
        stage = obj.stage_cost
        @assert stage isa LQStageCost "Requires LQ stage costs"
        
        terminal = obj.terminal_cost
        @assert terminal isa LQTerminalCost "Requires LQ terminal costs"
        
        push!(Q, stage.Q)
        push!(R, stage.R)
        push!(Qf, terminal.Qf)
        push!(M, stage.M)
        push!(q, stage.q)
        push!(r, stage.r)
    end
    
    return A, B_vec, Q, R, Qf, M, q, r
end