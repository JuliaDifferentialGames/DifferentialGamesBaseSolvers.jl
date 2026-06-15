# ALGAMES
**Reference:** Le Cleac'h, Schwager, Manchester (2021). ALGAMES: A Fast Augmented Lagrangian Solver for Constrained Dynamic Games. https://arxiv.org/abs/2104.08452
**Type:** Discrete-Time, Nonlinear, Constrained, Open-Loop Nash Equilibrium
**Description:** Augmented Lagrangian solver for open-loop (Generalized) Nash equilibria. Implements Algorithm 3 of Le Cleac'h et al. (2021) with an outer Augmented Lagrangian loop and inner Newton iterations. Converges to a Normalized Nash Equilibrium (NNE) for games with shared inequality constraints.

---

## Problem Class
- **Game Type:** Nonlinear Dynamic Games with Constraints
- **Players:** N-player (general)
- **Information Structure:** Open-loop (trajectory-based)
- **Dynamics:** Nonlinear (may be continuous or discrete)
- **Cost/Objective:** General nonlinear costs
- **Constraints:** Shared inequality constraints on states and controls

---

## Mathematical Formulation

### Dynamics
```
# Discrete-time:
x_{k+1} = f(x_k, u_{1,k}, ..., u_{N,k}, k)

# Continuous-time (auto-discretized):
ẋ(t) = f(x(t), u_1(t), ..., u_N(t), t)
```

### Cost for Player i
```
J_i = \sum_{k=0}^{N-1} L_i(x_k, u_{1,k}, ..., u_{N,k}, k) + Lf_i(x_N)
```

### Constraints
```
# Shared inequality constraints
h_j(x_k, u_{1,k}, ..., u_{N,k}, k) ≤ 0  for all j, k
```

### Assumptions
- Dynamics are twice continuously differentiable
- Cost functions are twice continuously differentiable
- Constraint functions are twice continuously differentiable
- Constraint qualification conditions hold

---

## Algorithm

### Core Idea
Uses an Augmented Lagrangian method with Newton inner iterations to solve for open-loop Nash equilibria with shared constraints. The outer loop handles the AL updates while the inner loop solves the stationarity conditions using Newton's method with line search.

### Steps
1. **Initialization:**
   - Initialize state and control trajectories X, U
   - Initialize Lagrange multipliers λ = 0
   - Initialize penalty parameter ρ = ρ_init

2. **Outer AL Loop (l = 1, ..., outer_iter):**
   - **Step 2a:** Newton Inner Loop (k = 1, ..., inner_iter):
     - Assemble residual G = [stationarity conditions; dynamics residual]
     - Assemble Jacobian H = ∂G/∂y
     - Solve δy = −H⁻¹G with regularization
     - Perform Armijo line search with backtracking
     - Update y = y + η·δy
     - Check inner convergence: ‖G‖ < tol_opt and ‖D‖ < tol_dyn
   - **Step 2b:** Dual Ascent Update:
     ```
     λ^(l+1) = λ^(l) + ρ^(l) · c(X^(l+1), U^(l+1))
     ```
   - **Step 2c:** Penalty Update:
     ```
     ρ^(l+1) = min(γ · ρ^(l), ρ_max)
     ```
   - **Step 2d:** Check Outer Convergence: primal + feasibility + cost stability

3. **Output:**
   ```
   (X*, U*)  # Open-loop Nash Equilibrium
   λ*        # Lagrange multipliers
   ```

### Key Equations
- Augmented Lagrangian: L_AL = L + λᵀ c + (ρ/2) ‖c‖²
- Stationarity conditions: ∇_{u_i} L_AL_i = 0 for all i
- Dynamics residual: D = [x_{k+1} - f(x_k, u_k, k)]
- Constraint violation: C = max_j c_j(x, u)
- Newton step: δy = −(H + reg·I)⁻¹ G
- Dual update: λ^(l+1) = λ^(l) + ρ^(l) · c(x^(l+1), u^(l+1))
- Penalty update: ρ^(l+1) = min(γ · ρ^(l), ρ_max)

---

## Implementation Details

### Inputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `game` | GameProblem | Constrained nonlinear game problem | - |
| `dynamics` | AbstractDynamics | Nonlinear dynamics function | - |
| `objectives` | Vector{Objective} | Player cost functions | - |
| `constraints` | Vector{Constraint} | Shared inequality constraints | - |

### Outputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `X` | Matrix | State trajectory | n × (N+1) |
| `U` | Vector{Matrix} | Control trajectories | U[i]: m_i × N |
| `λ` | Vector | Lagrange multipliers | Depends on constraints |
| `costs` | Vector | Final cost values | Scalar per player |

### Fields/Parameters
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| **Outer AL loop** |
| `outer_iter` | Int | 50 | Maximum AL iterations |
| `ρ_init` | Float64 | 1.0 | Initial penalty weight ρ⁽⁰⁾ |
| `ρ_increase` | Float64 | 10.0 | Geometric multiplier γ > 1 |
| `ρ_max` | Float64 | 1e8 | Cap on ρ |
| **Newton inner loop** |
| `inner_iter` | Int | 10 | Maximum Newton steps per outer iteration |
| `reg` | Float64 | 1e-3 | Tikhonov regularization on H |
| **Regularisation and line search** |
| `ls_iter` | Int | 20 | Maximum backtracks in Armijo line search |
| `ls_beta` | Float64 | 0.1 | Sufficient-decrease fraction β |
| `ls_tau` | Float64 | 0.5 | Step contraction τ |
| **Convergence tolerances** |
| `tol_opt` | Float64 | 1e-4 | ‖G_stationarity‖₁ / len |
| `tol_dyn` | Float64 | 1e-4 | ‖D‖₁ / len |
| `tol_con` | Float64 | 1e-3 | max constraint violation |
| **Warm-start** |
| `reset_duals` | Bool | true | If false, load λ, ρ, μ from WarmstartData |

---

## Theoretical Notes
- **Solution Concept:** Open-loop Nash equilibrium
- **Equilibrium Type:** With shared inequality constraints, converges to a Normalized Nash Equilibrium (NNE, Rosen 1967) because identical dual-ascent updates enforce equal multipliers on shared constraints (Eqs. 20-23 and Section 6.4 of the paper)
- **Convergence:** Guaranteed under standard AL conditions and constraint qualification
- **Complexity:** O(outer_iter · inner_iter · N · d²) where d is total decision dimension
- **Numerical Stability:** Regularization ensures Jacobian invertibility; line search ensures robustness

---

## Dependencies
- `LinearAlgebra`
- `DifferentialGamesBase`
- `ForwardDiff` (for automatic differentiation of residuals and Jacobians)

---

## Example Usage
```julia
using DifferentialGames

# Define constrained nonlinear game
dynamics = NonlinearDynamics((x, u, t) -> [x[1] + 0.1*u[1]; x[2] + 0.1*u[2]])
objectives = [
    Objective((x, u, t) -> 100x[1]^2 + u[1]^2),  # Player 1
    Objective((x, u, t) -> 100x[2]^2 + u[2]^2)   # Player 2
]

# Define shared constraints
constraints = [
    InequalityConstraint((x, u, t) -> x[1] + x[2] - 1.0),  # x1 + x2 ≤ 1
    ControlLimit((u, t) -> u[1], -0.5, 0.5),  # |u1| ≤ 0.5
    ControlLimit((u, t) -> u[2], -0.5, 0.5)   # |u2| ≤ 0.5
]

game = ConstrainedGameProblem(dynamics, objectives, constraints,
                               initial_state=[1.0; -1.0], horizon=50, dt=0.1)

# Solve with ALGAMES
solver = ALGAMES(outer_iter=50, ρ_init=1.0, ρ_increase=10.0,
                 inner_iter=10, reg=1e-3,
                 tol_opt=1e-4, tol_dyn=1e-4, tol_con=1e-3)
sol = solve(game, solver)

# Extract solution
X_opt = sol.state_trajectory  # n × (N+1)
U_opt = [sol.trajectories[i].controls for i in 1:2]  # Control trajectories
λ_opt = sol.multipliers  # Lagrange multipliers
costs_opt = sol.costs  # Final cost values

# Check constraint satisfaction
constraint_violations = sol.solver_info[:constraint_violations]
```

---

## References
- Le Cleac'h, Schwager, Manchester (2021). ALGAMES: A Fast Augmented Lagrangian Solver for Constrained Dynamic Games. https://arxiv.org/abs/2104.08452
- Rosen, J. B. (1967). Normalized Nash equilibria and non-zero-sum games. SIAM Journal on Control, 5(3), 363-374.

---

## Validation
- **Test Cases:** Validated on various constrained nonlinear games
- **Benchmarks:** Compared with unconstrained solutions and other constraint-handling methods
- **Numerical Tests:** Includes tests for convergence, constraint satisfaction, and dual updates

---

## Performance Tips
- Start with moderate penalty parameters (ρ_init) and increase as needed
- Use regularization (reg) for ill-conditioned problems
- Adjust convergence tolerances based on required accuracy
- For problems with many constraints, consider constraint prioritization
- Use warm-starting when solving similar problems sequentially
- Limit inner_iter for real-time applications

---

## Limitations
- Open-loop solutions only (no feedback)
- Requires smooth (twice differentiable) dynamics, costs, and constraints
- Local method: convergence depends on initialization
- May require careful tuning of penalty parameters
- No guarantee of global optimality for non-convex problems
- Computation time scales with problem size and number of constraints