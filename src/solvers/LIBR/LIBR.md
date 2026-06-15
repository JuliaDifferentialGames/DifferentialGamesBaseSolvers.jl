# LIBR
**Reference:** Miller, K. & Mitra, S. (2022). Multi-agent motion planning using differential games with lexicographic preferences. IEEE CDC 2022, pp. 5751-5757.
**Type:** Discrete-Time, Convex, Lexicographic Nash Equilibrium
**Description:** Lexicographic Iterative Best Response solver for LexicographicGameProblem. Implements Algorithm 1 of Miller & Mitra (2022) where each agent sequentially updates its strategy with two gradient-descent phases: first minimizing collision cost, then minimizing personal cost subject to collision cost constraint.

---

## Problem Class
- **Game Type:** Convex Game with Lexicographic Preferences
- **Players:** N-player (general)
- **Information Structure:** Open-loop (trajectory-based)
- **Dynamics:** Separable (each player has private dynamics)
- **Cost/Objective:** Lexicographic: Jᵢ(z) = (Jᵢᶜᵒˡ, Jᵢᵖᵉʳ) ordered lexicographically

---

## Mathematical Formulation

### Dynamics (Separable)
```
# Each player i has private dynamics
xᵢ_{k+1} = fᵢ(xᵢ_k, uᵢ_k, t_k)
```

### Lexicographic Cost for Player i
```
Jᵢ(z) = (Jᵢᶜᵒˡ(z), Jᵢᵖᵉʳ(z))  ∈ ℝ²
```
where z = (z₁, ..., z_N) with zᵢ = (Xᵢ, Uᵢ) and:
- Jᵢᶜᵒˡ: collision cost (higher priority)
- Jᵢᵖᵉʳ: personal cost (lower priority)

### Lexicographic Order
```
Jᵢ(z) < Jᵢ(z')  iff  Jᵢᶜᵒˡ(z) < Jᵢᶜᵒˡ(z')  OR  (Jᵢᶜᵒˡ(z) = Jᵢᶜᵒˡ(z') AND Jᵢᵖᵉʳ(z) < Jᵢᵖᵉʳ(z'))
```

### Assumptions
- Each player's dynamics fᵢ are convex in (xᵢ, uᵢ)
- Cost functions Jᵢᶜᵒˡ and Jᵢᵖᵉʳ are convex
- Private feasible sets Ωᵢ are convex
- Shared coupling constraint is affine: A x ≤ b

---

## Algorithm

### Core Idea
Each outer IBR iteration updates every agent sequentially with two gradient-descent phases. Phase 1 minimizes collision cost with all others fixed. Phase 2 minimizes personal cost subject to the collision cost not exceeding the Phase 1 optimum plus a small slack.

### Steps
1. **Initialization:**
   ```
   U[i] = 0 ∈ ℝ^{mᵢ × N} for all i
   X[i] = rollout(fᵢ, x_i0, U[i]) for all i
   ```

2. **Main Loop (iter = 1, 2, ..., max_iter):**
   - **Step 2a:** Store previous controls: prev_U = copy(U)
   - **Step 2b:** For each player i = 1, ..., N:
     - **Phase 1:** Minimize Jᵢᶜᵒˡ with z_{-i} fixed
       ```
       U[i] ← gradient_descent(Jᵢᶜᵒˡ(·, z_{-i}), U[i])
       J*_col = Jᵢᶜᵒˡ(z_i, z_{-i})
       X[i] ← rollout(fᵢ, x_i0, U[i])
       ```
     - **Phase 2:** Minimize Jᵢᵖᵉʳ subject to Jᵢᶜᵒˡ ≤ J*_col + slack
       ```
       U[i] ← gradient_descent(Jᵢᵖᵉʳ + ρ·max(0, Jᵢᶜᵒˡ−J*_col)², U[i])
       X[i] ← rollout(fᵢ, x_i0, U[i])
       ```
   - **Step 2c:** Check convergence: max_i ‖U[i] − prev_U[i]‖_∞ < tol

3. **Output:**
   ```
   z* = (X*, U*)  # Lexicographic Nash Equilibrium
   ```

### Key Equations
- Phase 1 update: U[i] = P_{Ωᵢ}[U[i] − τᵢ ∇_{U[i]} Jᵢᶜᵒˡ]
- Phase 2 update: U[i] = P_{Ωᵢ}[U[i] − τᵢ (∇_{U[i]} Jᵢᵖᵉʳ + ρ ∇_{U[i]} max(0, Jᵢᶜᵒˡ−J*_col)²)]
- Euler integration: x_{k+1} = x_k + dt · fᵢ(x_k, u_k, nothing, (k−1)·dt)

---

## Implementation Details

### Inputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `game` | LexicographicGameProblem | Game with lexicographic costs | - |
| `dynamics` | Vector{AbstractDynamics} | Private dynamics for each player | - |
| `collision_costs` | Vector{Function} | Jᵢᶜᵒˡ functions | Each: (X_i, U_i, z_{-i}) → ℝ |
| `personal_costs` | Vector{Function} | Jᵢᵖᵉʳ functions | Each: (X_i, U_i, z_{-i}) → ℝ |

### Outputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `X[i]` | Matrix | State trajectory for player i | nᵢ × (N+1) |
| `U[i]` | Matrix | Control trajectory for player i | mᵢ × N |
| `costs` | Vector{Tuple} | (Jᵢᶜᵒˡ, Jᵢᵖᵉʳ) for each player | Tuple{Float64, Float64} |

### Fields/Parameters
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| **IBR loop** |
| `max_iter` | Int | 100 | Maximum outer IBR iterations |
| `tol` | Float64 | 1e-4 | Convergence: max change in any player's control |
| **Per-phase gradient descent** |
| `inner_iter` | Int | 200 | Gradient steps per phase per player per IBR iteration |
| `step_size` | Float64 | 0.01 | Initial step size α₀ for backtracking line search |
| `ls_beta` | Float64 | 0.5 | Step contraction factor for line search |
| `ls_max_iter` | Int | 15 | Maximum backtracking steps |
| `ls_armijo` | Float64 | 1e-4 | Armijo sufficient-decrease constant c |
| **Phase 2 penalty** |
| `col_penalty` | Float64 | 1e3 | Quadratic penalty weight ρ for Jᵢᶜᵒˡ > J*_col |
| `col_slack` | Float64 | 1e-6 | Absolute slack: Phase 2 target is Jᵢᶜᵒˡ ≤ J*_col + δ |

---

## Theoretical Notes
- **Existence/Uniqueness:** Lexicographic Nash Equilibrium (L-NE) guaranteed to exist as a pure strategy by Proposition 2 of Miller & Mitra (2022)
- **Convergence:** Monotonic improvement in lexicographic order; convergence to L-NE under standard assumptions
- **Complexity:** O(max_iter · N · inner_iter · (nᵢ + mᵢ)) per player
- **Numerical Stability:** Line search ensures robustness; projection maintains feasibility

---

## Dependencies
- `LinearAlgebra`
- `DifferentialGamesBase`
- `ForwardDiff` (for gradient computation)

---

## Example Usage
```julia
using DifferentialGames

# Define lexicographic game
n_players = 2
n_states = [2, 2]  # Each player has 2D state
m_controls = [1, 1]  # Each player has 1D control

# Private dynamics for each player
dynamics = [
    NonlinearDynamics((x, u, t) -> [x[1] + 0.1*u[1]; x[2]]),  # Player 1
    NonlinearDynamics((x, u, t) -> [x[1]; x[2] + 0.1*u[1]])   # Player 2
]

# Lexicographic costs: (collision_cost, personal_cost)
collision_costs = [
    (X1, U1, z2) -> 100*norm(X1[:,end] - z2.X[:,end])^2,  # Player 1
    (X2, U2, z1) -> 100*norm(X2[:,end] - z1.X[:,end])^2   # Player 2
]
personal_costs = [
    (X1, U1, z2) -> sum(U1.^2),  # Player 1
    (X2, U2, z1) -> sum(U2.^2)   # Player 2
]

# Create lexicographic game problem
game = LexicographicGameProblem(dynamics, collision_costs, personal_costs,
                                 initial_states=[[1.0; 0.0], [0.0; 1.0]], horizon=50, dt=0.1)

# Solve with LIBR
solver = LIBR(max_iter=100, tol=1e-4, inner_iter=100, col_penalty=1e3)
sol = solve(game, solver)

# Extract trajectories
X1_opt = sol.trajectories[1].states  # Player 1 state trajectory
U1_opt = sol.trajectories[1].controls  # Player 1 control trajectory

# Get final costs
final_costs = sol.costs  # Vector of (Jᵢᶜᵒˡ, Jᵢᵖᵉʳ) tuples
```

---

## References
- Miller, K. & Mitra, S. (2022). Multi-agent motion planning using differential games with lexicographic preferences. IEEE CDC 2022, pp. 5751-5757.

---

## Validation
- **Test Cases:** Validated on multi-agent motion planning problems
- **Benchmarks:** Compared with standard Nash equilibrium solutions
- **Numerical Tests:** Includes tests for convergence, lexicographic ordering, and constraint satisfaction

---

## Performance Tips
- Start with smaller step_size for highly nonlinear problems
- Increase inner_iter for more accurate gradient descent
- Adjust col_penalty to balance collision avoidance and personal cost
- Use col_slack to control how strictly collision cost constraint is enforced
- For problems with known good initial guesses, provide initial trajectories

---

## Limitations
- Requires convex dynamics and costs
- Separable dynamics only (no shared state)
- Open-loop strategies (no feedback)
- No handling of hard constraints (only penalty-based)
- Lexicographic approach may not capture all Pareto-optimal solutions