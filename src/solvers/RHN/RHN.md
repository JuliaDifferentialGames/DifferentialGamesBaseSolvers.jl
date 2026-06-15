# RecedingHorizonNash
**Reference:** Mattingley, Wang, Boyd (2011) — Receding Horizon Control; Mahajan et al. (2026) — Conic-TinyMPC; Laine et al. (2023) — GFNE
**Type:** Discrete-Time, Model Predictive Control, Feedback Nash Equilibrium
**Description:** Receding-horizon wrapper for any Nash equilibrium solver. At each simulation step, the sub-game is solved over a fixed prediction window, only the first optimal control is applied, and the previous solution is shifted to warm-start the next solve. This provides a closed-loop Nash strategy via the receding-horizon loop.

---

## Problem Class
- **Game Type:** General (wraps any supported game type)
- **Players:** N-player (general)
- **Information Structure:** Feedback (closed-loop via MPC)
- **Dynamics:** Any supported by inner solver
- **Cost/Objective:** Any supported by inner solver
- **Horizon:** Finite prediction window with receding horizon

---

## Mathematical Formulation

### Receding Horizon Loop
```
for t = 1, ..., N_sim:
  1. sub_game = with_initial_state(game_window, X[:,t])
  2. sol = solve(sub_game, inner_solver; warmstart=ws)
  3. u_t = apply_strategy(sol.strategy, X[:,t], 1)
  4. X[:,t+1] = dynamics(X[:,t], u_t, t)
  5. ws = shift_warmstart(sol)
```

### Warm-Start Shift
For FeedbackStrategy:
```
# Shift gains, feedforward, and nominal trajectory left by 1 step
P_shifted[i][k] = P_old[i][k+1] for k = 1,...,N-1
P_shifted[i][N] = P_old[i][N]  # Repeat last element

α_shifted[i][k] = α_old[i][k+1] for k = 1,...,N-1
α_shifted[i][N] = α_old[i][N]  # Repeat last element
```

For OpenLoopStrategy:
```
U_shifted[i][k] = U_old[i][k+1] for k = 1,...,N-1
U_shifted[i][N] = U_old[i][N]  # Repeat last element
```

---

## Algorithm

### Core Idea
Provides a closed-loop control strategy by repeatedly solving finite-horizon optimal control problems over a receding window. The warm-start shifting ensures that each sub-problem solve starts from a good initial guess, significantly reducing computation time for nonlinear solvers.

### Steps
1. **Initialization:**
   - Initialize current state X[:,1] = x0
   - Initialize warm-start data (optional)

2. **Simulation Loop (t = 1, ..., N_sim):**
   - **Step 2a:** Create sub-game with current state as initial condition
   - **Step 2b:** Solve sub-game using inner solver with warm-start
   - **Step 2c:** Extract first control action from solution
   - **Step 2d:** Apply control and propagate state
   - **Step 2e:** Shift solution to create warm-start for next iteration

3. **Output:**
   - Full state and control trajectories over simulation horizon
   - Final warm-start data for potential continuation

### Key Equations
- Sub-game: game_t = RecedingHorizonNashProblem(game_template, X[:,t], prediction_horizon)
- Control application: u_t = sol.strategy.P[i][1] * X[:,t] + sol.strategy.α[i][1]
- State propagation: X[:,t+1] = f(X[:,t], u_t, t)
- Warm-start shift: ws_t = shift_warmstart(sol_t)

---

## Implementation Details

### Inputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `inner_solver` | GameSolver | Any GameSolver that handles GameProblem | - |
| `game_template` | GameProblem | Template game problem with full horizon | - |
| `x0` | Vector | Initial state | n |
| `N_sim` | Int | Simulation horizon (number of steps) | Scalar |

### Outputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `X` | Matrix | State trajectory over simulation | n × (N_sim+1) |
| `U` | Vector{Matrix} | Control trajectories for each player | U[i]: m_i × N_sim |
| `strategy` | FeedbackStrategy | Current feedback strategy | P[i][k]: m_i×n, α[i][k]: m_i |
| `warmstart` | WarmstartData | Final warm-start for continuation | - |

### Fields/Parameters
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `inner_solver` | GameSolver | Required | Any GameSolver (FNELQ, iLQGames, ALGAMES, etc.) |
| `warm_start` | Bool | `true` | Whether to shift previous solution as warm-start |
| `verbose_inner` | Bool | `false` | Whether to pass verbose=true to inner solves |

---

## Theoretical Notes
- **Solution Concept:** Closed-loop Nash strategy via receding horizon MPC
- **Equilibrium Type:** Inherits from inner solver (FeedbackNash, OpenLoopNash, etc.)
- **Convergence:** Depends on inner solver; warm-starting typically improves convergence
- **Complexity:** O(N_sim × inner_complexity) where inner_complexity depends on inner solver
- **Numerical Stability:** Warm-start shifting maintains feasibility and provides good initialization

---

## Dependencies
- `LinearAlgebra`
- `DifferentialGamesBase`
- Any inner solver used (FNELQ, iLQGames, ALGAMES, etc.)

---

## Example Usage
```julia
using DifferentialGames

# Define the full-horizon game template
dynamics = NonlinearDynamics((x, u, t) -> [x[1] + 0.1*u[1]; x[2] + 0.1*u[2]])
objectives = [
    Objective((x, u, t) -> 100x[1]^2 + u[1]^2),  # Player 1
    Objective((x, u, t) -> 100x[2]^2 + u[2]^2)   # Player 2
]
game_template = GameProblem(dynamics, objectives, initial_state=[0.0; 0.0], horizon=20, dt=0.1)

# Create receding horizon problem
N_sim = 100  # Simulation steps
x0 = [1.0; -1.0]  # Initial state

# Choose inner solver
inner_solver = iLQGames(max_iter=50, ε_conv=1e-4)

# Create RHN solver
solver = RecedingHorizonNash(inner_solver; warm_start=true, verbose_inner=false)

# Create RHN problem
rhn_problem = RecedingHorizonNashProblem(game_template, x0, N_sim)

# Solve
sol = solve(rhn_problem, solver)

# Extract results
X_traj = sol.state_trajectory  # n × (N_sim+1)
U_traj = [sol.trajectories[i].controls for i in 1:2]  # Control trajectories

# Get final feedback strategy (for current window)
P_final = sol.strategy.P  # Feedback gains
α_final = sol.strategy.α  # Affine terms

# Continue simulation with new initial state
x_new = [0.5; -0.5]
rhn_problem_new = RecedingHorizonNashProblem(game_template, x_new, N_sim)
sol_new = solve(rhn_problem_new, solver; warmstart=sol.warmstart)
```

---

## References
- Mattingley, J., Wang, A., & Boyd, S. (2011). Receding horizon control. Annual Review of Control, Robotics, and Autonomous Systems.
- Mahajan, A., et al. (2026). Conic-TinyMPC: A minimalist MPC solver for real-time applications.
- Laine, T., et al. (2023). Generalized Feedback Nash Equilibrium for dynamic games.

---

## Validation
- **Test Cases:** Validated on various game types with different inner solvers
- **Benchmarks:** Compared with full-horizon solutions for consistency
- **Numerical Tests:** Includes tests for warm-start shifting and closed-loop stability

---

## Performance Tips
- Use warm_start=true for nonlinear inner solvers (iLQGames, ALGAMES) to significantly reduce computation time
- Set verbose_inner=false for production use to reduce overhead
- Choose inner solver based on problem characteristics:
  - FNELQ: Fast for LQ games
  - iLQGames: For nonlinear games without constraints
  - ALGAMES: For open-loop constrained games
- For real-time applications, limit max_iter in inner solver

---

## Limitations
- Closed-loop performance depends on inner solver quality
- Warm-start may not be effective for problems with rapidly changing dynamics
- No stability guarantees (depends on inner solver and problem)
- Computation time scales with prediction horizon
- Memory usage can be high for large prediction horizons