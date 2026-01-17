using DifferentialGamesBase

"""
    FNELQ <: GameSolver

Solves N-player discrete-time LQ games for feedback Nash equilibrium via 
dynamic programming.

Based on Başar & Olsder (1999), Theorem 6.17 (discrete-time case). Computes 
feedback Nash equilibrium by solving coupled matrix equations backward in time.

# Example
```julia
game = LQGameProblem(A, B, Q, R, Qf, x0, 10.0; dt=0.1)
solver = FNELQ()
sol = solve(game, solver; verbose=true)

# Extract feedback gains
P_gains = sol.solver_info[:feedback_gains]  # Proportional gains P(k)
α_gains = sol.solver_info[:affine_gains]    # Affine terms α(k)
```

# References
Başar, T., & Olsder, G. J. (1999). Dynamic Noncooperative Game Theory (2nd ed.). 
SIAM. Section 6.3, Theorem 6.17 (Equation 6.17a).
"""
struct FNELQ <: GameSolver
    check_singularity::Bool
    rcond_threshold::Float64
    
    function FNELQ(;
        check_singularity::Bool = true,
        rcond_threshold::Float64 = 1e-10
    )
        new(check_singularity, rcond_threshold)
    end
end

solver_capabilities(::Type{FNELQ}) = [
    :LQGame, 
    :FeedbackPolicies,
    :UnconstrainedGame,
    :DiscreteTime
]

"""
    _solve(game::GameProblem{T}, solver::FNELQ, warmstart, verbose::Bool)

Solve discrete-time unconstrained LQ game for feedback Nash equilibrium.
"""
function _solve(
    game::GameProblem{T},
    solver::FNELQ,
    warmstart::Union{Nothing, WarmstartData},
    verbose::Bool
) where {T}
    
    # Extract LQ structure
    A, B_vec, Q, R, Qf, M, q, r = extract_lq_matrices(game)
    n = size(A, 1)
    n_players = game.n_players
    
    # Validate discrete-time
    @assert game.time_horizon isa DiscreteTime "FNELQ requires discrete-time formulation"
    
    # Validate unconstrained
    @assert is_unconstrained(game) "FNELQ only handles unconstrained games"
    
    # Validate no cross terms (M) or linear terms for now
    if any(norm(M[i]) > eps(T) for i in 1:n_players)
        @warn "FNELQ: Cross terms M ≠ 0 not yet implemented, ignoring"
    end
    has_affine_terms = any(norm(q[i]) > eps(T) || norm(r[i]) > eps(T) for i in 1:n_players)
    
    # Control dimensions and indexing
    control_dims = [size(B_vec[i], 2) for i in 1:n_players]
    total_controls = sum(control_dims)
    control_ranges = [sum(control_dims[1:i-1])+1:sum(control_dims[1:i]) for i in 1:n_players]
    
    # Concatenate B matrices
    B = hcat(B_vec...)  # [n × total_controls]
    
    # Time discretization
    N = Int(round(game.time_horizon.tf / game.time_horizon.dt))
    dt = game.time_horizon.dt
    
    verbose && @info "FNELQ: Solving $(n_players)-player discrete LQ game" n_states=n time_steps=N dt=dt control_dims=control_dims
    
    t_start = time()
    
    # Initialize cost-to-go (dynamic programming value function)
    Z = [copy(Qf[i]) for i in 1:n_players]  # Quadratic term: Vᵢ(x) = xᵀZᵢx + ζᵢᵀx
    ζ = has_affine_terms ? [zeros(T, n) for i in 1:n_players] : nothing
    
    # Storage for feedback strategies
    P_history = [Matrix{T}[] for _ in 1:n_players]  # P[i][k] is player i's gain at time k
    α_history = has_affine_terms ? [Vector{T}[] for _ in 1:n_players] : nothing
    
    # Storage for cost-to-go matrices (for diagnostics)
    Z_history = [Matrix{T}[] for _ in 1:n_players]
    for i in 1:n_players
        push!(Z_history[i], copy(Z[i]))
    end
    
    # Condition number history
    rcond_history = Float64[]
    
    # Backward pass: dynamic programming
    for k in N:-1:1
        # Build coupled S matrix and Y matrix
        S = zeros(T, total_controls, total_controls)
        YP = zeros(T, total_controls, n)
        Yα = has_affine_terms ? zeros(T, total_controls) : nothing
        
        for i in 1:n_players
            rng_i = control_ranges[i]
            Bi = B_vec[i]
            Ri = R[i]
            Zi = Z[i]
            
            # BᵢᵀZᵢBⱼ terms
            BiZi = Bi' * Zi
            
            for j in 1:n_players
                rng_j = control_ranges[j]
                Bj = B_vec[j]
                
                # Diagonal blocks: Rᵢ + BᵢᵀZᵢBᵢ
                # Off-diagonal blocks: BᵢᵀZᵢBⱼ
                S[rng_i, rng_j] = (i == j ? Ri : zeros(T, length(rng_i), length(rng_j))) + BiZi * Bj
            end
            
            # Right-hand side for proportional gains: YP[uⁱ, :] = BᵢᵀZᵢA
            YP[rng_i, :] = BiZi * A
            
            # Right-hand side for affine terms: Yα[uⁱ] = Bᵢᵀζᵢ + rᵢ
            if has_affine_terms
                Yα[rng_i] = Bi' * ζ[i] + r[i]
            end
        end
        
        # Check conditioning
        if solver.check_singularity
            S_rcond = 1.0 / cond(S)
            push!(rcond_history, S_rcond)
            
            if S_rcond < solver.rcond_threshold
                @warn "FNELQ: S matrix poorly conditioned at k=$k" rcond=S_rcond
            end
        end
        
        # Solve S * P = YP for feedback gains
        P = S \ YP  # [total_controls × n]
        
        # Extract per-player gains
        P_k = [P[control_ranges[i], :] for i in 1:n_players]
        
        # Store (prepending since going backward)
        for i in 1:n_players
            pushfirst!(P_history[i], copy(P_k[i]))
        end
        
        # Handle affine terms if present
        α_k = nothing
        β = nothing
        if has_affine_terms
            # Solve S * α = Yα where Yα[uⁱ] = Bᵢᵀζᵢ + rᵢ
            α = S \ Yα
            α_k = [α[control_ranges[i]] for i in 1:n_players]
            
            for i in 1:n_players
                pushfirst!(α_history[i], copy(α_k[i]))
            end
            
            β = -B * α
        end
        
        # Update cost-to-go for next iteration (k-1)
        F = A - B * P
        
        for i in 1:n_players
            Pi = P_k[i]
            Ri = R[i]
            
            # Quadratic term: Zᵢ(k) = FᵀZᵢ(k+1)F + Qᵢ + PᵢᵀRᵢPᵢ
            PRi = Pi' * Ri
            Z[i] = F' * Z[i] * F + Q[i] + PRi * Pi
            
            # Ensure symmetry
            Z[i] = (Z[i] + Z[i]') / 2
            
            pushfirst!(Z_history[i], copy(Z[i]))
            
            # Affine term update if present (Başar & Olsder Eq. 6.17c)
            # ζᵢ(k) = Fᵀ(ζᵢ(k+1) + Zᵢ(k+1)β) + qᵢ + PᵢᵀRᵢαᵢ - Pᵢᵀrᵢ
            if has_affine_terms
                αi = α_k[i]
                ζ[i] = F' * (ζ[i] + Z[i] * β) + q[i] + PRi * αi - Pi' * r[i]
            end
        end
        
        verbose && (k % 10 == 0 || k == N) && @info "  Timestep $(N-k+1)/$(N)" rcond=get(rcond_history, length(rcond_history), NaN)
    end
    
    backward_time = time() - t_start
    verbose && @info "FNELQ: Backward pass complete" time=backward_time
    
    # Forward pass: simulate with feedback policies
    t_forward_start = time()
    
    x_traj = zeros(T, n, N+1)
    x_traj[:, 1] = game.initial_state
    
    u_trajs = [zeros(T, control_dims[i], N+1) for i in 1:n_players]
    
    for k in 1:N
        x_k = x_traj[:, k]
        
        # Compute controls: uᵢ(k) = -Pᵢ(k)·x(k) - αᵢ(k)
        u_k = Vector{Vector{T}}(undef, n_players)
        for i in 1:n_players
            u_k[i] = -P_history[i][k] * x_k
            if has_affine_terms
                u_k[i] -= α_history[i][k]
            end
            u_trajs[i][:, k] = u_k[i]
        end
        
        # Dynamics: x(k+1) = A·x(k) + Σᵢ Bᵢ·uᵢ(k)
        u_total = sum(B_vec[i] * u_k[i] for i in 1:n_players)
        x_traj[:, k+1] = A * x_k + u_total
    end
    
    forward_time = time() - t_forward_start
    total_time = time() - t_start
    
    verbose && @info "FNELQ: Forward pass complete" time=forward_time total_time=total_time
    
    # Build trajectories
    t_vec = range(T(0), game.time_horizon.tf, length=N+1)
    trajectories = Trajectory{T}[]
    for i in 1:n_players    
        u_traj_i = u_trajs[i][:, 1:N]  # [mᵢ × N]
        traj = Trajectory(
            i,
            x_traj,              # [n × (N+1)] - shared state for all players
            u_traj_i,          # [mᵢ × N] - player i's controls
            collect(t_vec),      # [(N+1)] time vector
            game.time_horizon.tf # terminal time
        )
        push!(trajectories, traj)
    end
    
    # Compute costs
    costs = zeros(T, n_players)
    for i in 1:n_players
        stage_cost_sum = zero(T)
        for k in 1:N
            xk = x_traj[:, k]
            uk = u_trajs[i][:, k]
            stage_cost_sum += xk' * Q[i] * xk + uk' * R[i] * uk
            if has_affine_terms
                stage_cost_sum += 2 * q[i]' * xk + 2 * r[i]' * uk
            end
        end
        xf = x_traj[:, N+1]
        terminal_cost = xf' * Qf[i] * xf
        costs[i] = stage_cost_sum + terminal_cost
    end
    
    # Build solver info
    solver_info = Dict{Symbol, Any}(
        :feedback_gains => P_history,
        :cost_to_go_matrices => Z_history,
        :costs => costs,
        :backward_pass_time => backward_time,
        :forward_pass_time => forward_time,
        :rcond_history => rcond_history
    )
    
    if has_affine_terms
        solver_info[:affine_gains] = α_history
    end
    
    # Always converges in one backward pass (no iteration)
    converged = true
    if solver.check_singularity && any(rcond_history .< solver.rcond_threshold)
        converged = false
        @warn "FNELQ: Solution may be unreliable due to ill-conditioned S matrices"
    end
    
    return GameSolution(
        game,
        trajectories;
        equilibrium_type = :FeedbackNash,
        converged = converged,
        iterations = 1,  # Single backward pass
        solve_time = total_time,
        solver_info = solver_info
    )
end