# iLQGames
**Reference:** Fridovich-Keil, D., Ratner, E., Peters, L., Dragan, A. D., & Tomlin, C. J. (2020). Efficient Iterative Linear-Quadratic Approximations for Nonlinear Multi-Player General-Sum Differential Games. ICRA 2020. arXiv:1909.04694
**Type:** Discrete-Time, Nonlinear, Non-Zero-Sum, Feedback Nash Equilibrium
**Description:** Iterative Linear-Quadratic Games solver for nonlinear multi-player general-sum differential games. Finds a feedback Nash equilibrium of successive LQ approximations by linearising the dynamics and quadraticising the costs around the current nominal trajectory, then solving the resulting LQ subgame exactly via FNELQ.

---

## Problem Class
- **Game Type:** Nonlinear Differential Games with LQ approximations
- **Players:** N-player (general)
- **Information Structure:** Feedback (closed-loop)
- **Dynamics:** Nonlinear (automatically discretized if continuous)
- **Cost/Objective:** General nonlinear costs (quadraticised around nominal trajectory)

---

## Mathematical Formulation

### Dynamics
```
# Continuous-time (auto-discretized):
ẋ(t) = f(x(t), u_1(t), ..., u_N(t), t)

# Discrete-time (after discretization):
x_{k+1} = f_d(x_k, u_{1,k}, ..., u_{N,k}, k)
```

### Cost for Player i
```
J_i = \sum_{k=0}^{N-1} L_i(x_k, u_{1,k}, ..., u_{N,k}, k) + Lf_i(x_N)
```
where L_i and Lf_i are general nonlinear cost functions that are quadraticised around the nominal trajectory.

### Assumptions
- Dynamics are twice continuously differentiable
- Cost functions are twice continuously differentiable
- Initial nominal trajectory is feasible (or can be initialized)

---

## Algorithm

### Core Idea
Iteratively linearizes the dynamics and quadraticizes the costs around the current nominal trajectory, then solves the resulting LQ subgame exactly. This provides a feedback Nash equilibrium for the local LQ approximation at the fixed point.

### Steps
1. **Initialization:**
   - Initialize nominal state trajectory X and control trajectory U
   - Set initial regularization parameter μ = μ_init

2. **Main Loop (l = 1, 2, ..., max_iter):**
   - **Step 2a:** Linearize dynamics around (X, U) to get LTV system
   - **Step 2b:** Quadraticize costs around (X, U) to get LQ cost approximations
   - **Step 2c:** Assemble LQ subgame with regularization μ
   - **Step 2d:** Solve LQ subgame exactly using FNELQ solver
   - **Step 2e:** Extract feedback strategy from FNELQ solution
   - **Step 2f:** Perform line search to find step size η
   - **Step 2g:** Update nominal trajectory: X_new, U_new = rollout with η
   - **Step 2h:** Check convergence: max(‖X_new - X‖_∞, ‖U_new - U‖_∞) < ε_conv
   - **Step 2i:** Update regularization: adjust μ based on convergence

3. **Feedback Strategy:**
   ```
   u_{i,k} = -P_i(k) x_k - α_i(k)  (from last FNELQ solve)
   ```

### Key Equations
- Linearization: `A(k) = ∂f/∂x`, `B_i(k) = ∂f/∂u_i` at (X[k], U[k])
- Quadraticization: `Q_i(k) = ∂²L_i/∂x²`, `R_i(k) = ∂²L_i/∂u_i²`, etc.
- Regularization: μ added to diagonal of S matrix in FNELQ
- Line search: Geometric backtracking with factor β

---

## Implementation Details

### Inputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `game` | GameProblem | Nonlinear game problem | - |
| `dynamics` | AbstractDynamics | Nonlinear dynamics function | - |
| `objectives` | Vector{Objective} | Player cost functions | - |

### Outputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `strategy` | FeedbackStrategy | Feedback gains and affine terms | P[i][k]: m_i×n, α[i][k]: m_i |
| `trajectories` | Vector{Trajectory} | State and control trajectories | X: n×(N+1), U_i: m_i×N |
| `solver_info` | Dict | Additional solver information | Includes open_loop_strategy |

### Fields/Parameters
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `max_iter` | Int | 200 | Maximum outer iterations |
| `ε_conv` | Float64 | 0.05 | Trajectory-change convergence threshold |
| `β` | Float64 | 0.5 | Line search backtrack factor |
| `η_min` | Float64 | β^20 | Minimum line search step |
| `max_state_step` | Float64 | 1.0 | Maximum ‖X_new - X_curr‖_∞ accepted in line search |
| `μ_init` | Float64 | 1.0 | Initial S-regularization |
| `μ_max` | Float64 | 1e6 | Maximum regularization |
| `μ_scale` | Float64 | 10.0 | Regularization growth factor on ill-conditioning |
| `μ_decay` | Float64 | 0.5 | Regularization decay factor on accepted step |
| `discretization` | AbstractDiscretizationMethod | ZOHDiscretization() | Method for auto-discretising continuous dynamics |

---

## Theoretical Notes
- **Existence/Uniqueness:** Convergence to feedback Nash equilibrium of LQ approximation at fixed point (§IV-B, Fridovich-Keil et al. 2020)
- **Convergence:** Guaranteed under Assumptions 1-3 of Fridovich-Keil et al. (2020) when step-sizes satisfy Theorem 2
- **Complexity:** O(max_iter · N · (Σm_i)^3) where N is horizon, Σm_i is total control dimension
- **Numerical Stability:** Regularization ensures S matrix is well-conditioned; line search ensures robustness

---

## Dependencies
- `LinearAlgebra`
- `DifferentialGamesBase`
- `ForwardDiff` (for automatic differentiation)
- FNELQ solver (used internally)

---

## Example Usage
```julia
using DifferentialGames

# Define nonlinear game
dynamics = NonlinearDynamics((x, u, t) -> x + 0.1 * [u[1] + u[2]; u[1] - u[2]])
objectives = [
    Objective((x, u, t) -> 100x[1]^2 + u[1]^2),  # Player 1
    Objective((x, u, t) -> 100x[2]^2 + u[2]^2)   # Player 2
]
game = GameProblem(dynamics, objectives, initial_state=[1.0; -1.0], horizon=50, dt=0.1)

# Solve with iLQGames
solver = iLQGames(max_iter=100, ε_conv=1e-4, μ_init=1.0)
sol = solve(game, solver)

# Extract feedback strategy
P_1 = sol.strategy.P[1]  # Feedback gains for player 1
α_1 = sol.strategy.α[1]  # Affine terms for player 1

# Get converged nominal trajectory
X_opt = sol.trajectories[1].states
U_opt = sol.trajectories[1].controls
```

---

## References
- Fridovich-Keil, D., Ratner, E., Peters, L., Dragan, A. D., & Tomlin, C. J. (2020). Efficient Iterative Linear-Quadratic Approximations for Nonlinear Multi-Player General-Sum Differential Games. ICRA 2020. arXiv:1909.04694

---

## Validation
- **Test Cases:** Validated on various nonlinear games including double integrator, simple pursuit
- **Benchmarks:** Compared with analytical solutions for special cases
- **Numerical Tests:** Includes tests for convergence, regularization, and line search

---

## Performance Tips
- Start with higher regularization (μ_init) for highly nonlinear problems
- Reduce ε_conv for higher accuracy (at cost of more iterations)
- Use ZOHDiscretization for continuous-time systems
- Pre-allocate memory for large problems by providing initial trajectories
- For problems with known good initial guesses, use warm-starting

---

## Limitations
- Requires smooth (twice differentiable) dynamics and costs
- Local method: convergence depends on initialization
- May get stuck in local minima for highly non-convex problems
- Regularization may affect solution accuracy for very small μ
- No constraint handling (unconstrained game only)