# YiPavel
**Reference:** Yi, P. & Pavel, L. (2017). A distributed primal-dual algorithm for computation of generalized Nash equilibria via operator splitting. IEEE CDC 2017, pp. 3841-3846.
**Type:** Discrete-Time, Convex, Generalized Nash Equilibrium
**Description:** Distributed primal-dual operator-splitting GNE solver for ConvexGameProblem. Implements Algorithm 1 of Yi & Pavel (2017) for computing the variational GNE (v-NE) of an N-player convex game with a globally shared affine coupling constraint. Uses simultaneous updates of primal decisions, auxiliary consensus variables, and local copies of the coupling multiplier.

---

## Problem Class
- **Game Type:** Convex Game with Shared Coupling Constraints
- **Players:** N-player (general)
- **Information Structure:** Distributed (fully-connected communication graph)
- **Dynamics:** Separable (each player has private decision variables)
- **Cost/Objective:** Convex private costs with coupling
- **Constraints:** Shared affine coupling constraint Ax ≤ b

---

## Mathematical Formulation

### Decision Variables
```
# Each player i controls
xᵢ ∈ Ωᵢ ⊂ ℝ^{dᵢ}  (private decision variables)
```

### Cost for Player i
```
fᵢ(x₁, ..., x_N)  (private cost, convex in xᵢ, may couple across players)
```

### Coupling Constraint
```
Σᵢ Aᵢ xᵢ ≤ b  (shared affine constraint)
```

### Assumptions
- Each player's feasible set Ωᵢ is convex and closed
- Each cost function fᵢ is convex in xᵢ
- Coupling constraint is affine
- Communication graph is fully-connected
- Step-sizes satisfy Theorem 2 (Lemma 3) of Yi & Pavel (2017)

---

## Algorithm

### Core Idea
Distributed primal-dual operator splitting algorithm where each player simultaneously updates its primal decision, auxiliary consensus variable, and local copy of the coupling multiplier. The algorithm enforces consensus on the coupling multiplier through Laplacian-based dynamics while maintaining primal feasibility.

### Steps
1. **Initialization:**
   ```
   xᵢ = 0 ∈ ℝ^{dᵢ} for all i
   zᵢ = 0 ∈ ℝ^{m} for all i  (auxiliary consensus variables)
   λᵢ = 0 ∈ ℝ^{m} for all i  (local copies of coupling multiplier)
   ```

2. **Main Loop (k = 0, 1, 2, ..., max_iter):**
   - **Primal Update (simultaneous for all i):**
     ```
     xᵢ^{k+1} = P_{Ωᵢ}[ xᵢ^k − τᵢ (∇_{xᵢ} fᵢ(x^k) − Aᵢᵀ λᵢ^k) ]
     ```
   - **Auxiliary Update (consensus on λ):**
     ```
     zᵢ^{k+1} = zᵢ^k + νᵢ Σ_{j≠i} (λᵢ^k − λⱼ^k)
     ```
   - **Dual Update (simultaneous for all i):**
     ```
     λᵢ^{k+1} = P_{ℝ₊ᵐ}[ λᵢ^k − σᵢ (
                   Aᵢ (2xᵢ^{k+1} − xᵢ^k) − bᵢ
                   + Σ_{j≠i} [2(zᵢ^{k+1}−zⱼ^{k+1}) − (zᵢ^k−zⱼ^k)]
                   + Σ_{j≠i} (λᵢ^k − λⱼ^k)) ]
     ```

3. **Convergence Check:**
   ```
   Δ = max_i ‖xᵢ^{k+1} − xᵢ^k‖_∞ < tol
   ```

### Key Equations
- Primal update: Projected gradient descent on local Lagrangian
- Auxiliary update: PI-like dynamics accumulating Laplacian of λ
- Dual update: Projected ascent with extrapolated constraint residual and consensus correction
- Consensus: Laplacian L = D - A where D is degree matrix, A is adjacency matrix

---

## Implementation Details

### Inputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `game` | ConvexGameProblem | N-player convex game | - |
| `cost_fns` | Vector{Function} | fᵢ functions | Each: (x₁,...,x_N) → ℝ |
| `coupling_A` | Matrix | Coupling constraint matrix | m × Σdᵢ |
| `coupling_b` | Vector | Coupling constraint RHS | m |
| `lb` | Vector{Vector} | Lower bounds for each player | lb[i]: dᵢ |
| `ub` | Vector{Vector} | Upper bounds for each player | ub[i]: dᵢ |

### Outputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `x[i]` | Vector | Decision variables for player i | dᵢ |
| `λ[i]` | Vector | Local coupling multiplier for player i | m |
| `z[i]` | Vector | Auxiliary consensus variable for player i | m |
| `costs` | Vector | Final cost values for each player | Scalar per player |

### Fields/Parameters
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `max_iter` | Int | 5000 | Maximum number of iterations |
| `tol` | Float64 | 1e-6 | Convergence tolerance ‖Δx‖_∞ < tol |
| `τ` | Float64 | 0.05 | Primal step size τᵢ (same for all players) |
| `ν` | Float64 | 0.02 | Auxiliary (z) step size νᵢ (same for all players) |
| `σ` | Float64 | 0.05 | Dual step size σᵢ (same for all players) |

---

## Theoretical Notes
- **Solution Concept:** Variational Generalized Nash Equilibrium (v-GNE)
- **Existence/Uniqueness:** Convergence to v-GNE guaranteed under Assumptions 1-3 of Yi & Pavel (2017) when step-sizes satisfy Theorem 2 (Lemma 3)
- **Convergence:** Linear convergence rate under standard assumptions
- **Complexity:** O(max_iter · N · d) where d is total decision dimension
- **Numerical Stability:** Projection ensures feasibility; step-sizes can be tuned for robustness

---

## Dependencies
- `LinearAlgebra`
- `DifferentialGamesBase`
- `ForwardDiff` (for gradient computation)

---

## Example Usage
```julia
using DifferentialGames

# Define convex game
n_players = 3
n_vars = [2, 2, 2]  # Each player has 2 decision variables

# Cost functions (convex, may couple across players)
cost_fns = [
    (x1, x2, x3) -> 100x1[1]^2 + x1[2]^2 + 0.1*(x1[1]-x2[1])^2 + 0.1*(x1[2]-x3[2])^2,
    (x1, x2, x3) -> 80x2[1]^2 + x2[2]^2 + 0.1*(x2[1]-x1[1])^2 + 0.1*(x2[2]-x3[1])^2,
    (x1, x2, x3) -> 60x3[1]^2 + x3[2]^2 + 0.1*(x3[1]-x1[2])^2 + 0.1*(x3[2]-x2[2])^2
]

# Coupling constraint: x1[1] + x2[1] + x3[1] + x1[2] + x2[2] + x3[2] ≤ 1
coupling_A = [1 0 1 0 1 0; 0 1 0 1 0 1]  # 2 constraints
coupling_b = [1.0; 1.0]

# Bounds
lb = [[-1.0; -1.0], [-1.0; -1.0], [-1.0; -1.0]]
ub = [[1.0; 1.0], [1.0; 1.0], [1.0; 1.0]]

# Create convex game problem
game = ConvexGameProblem(cost_fns, n_vars, lb, ub)

# Solve with YiPavel
solver = YiPavel(max_iter=5000, tol=1e-6, τ=0.05, ν=0.02, σ=0.05)
sol = solve(game, solver;
            coupling_A=coupling_A,
            coupling_b=coupling_b,
            coupling_leq=true)

# Extract solution
x_opt = [sol.x[i] for i in 1:3]  # Decision variables for each player
λ_opt = sol.multipliers  # Coupling multipliers
costs_opt = sol.costs  # Final cost values

# Check constraint satisfaction
constraint_value = coupling_A * vcat(x_opt...) - coupling_b
```

---

## References
- Yi, P. & Pavel, L. (2017). A distributed primal-dual algorithm for computation of generalized Nash equilibria via operator splitting. IEEE CDC 2017, pp. 3841-3846.

---

## Validation
- **Test Cases:** Validated on various convex GNE problems with shared constraints
- **Benchmarks:** Compared with centralized GNE solvers
- **Numerical Tests:** Includes tests for convergence, consensus, and constraint satisfaction

---

## Performance Tips
- Start with smaller step-sizes (τ, ν, σ) for problems with strong coupling
- Increase max_iter for higher accuracy
- Adjust step-sizes based on problem scaling (larger for well-scaled problems)
- For problems with many players, consider hierarchical or sparse communication
- Use warm-starting when solving similar problems sequentially

---

## Limitations
- Requires convex costs and feasible sets
- Shared coupling constraint must be affine
- Fully-connected communication graph (no sparse communication support)
- Same step-sizes for all players (no per-player tuning)
- No handling of nonlinear constraints
- Distributed but requires synchronization (not fully asynchronous)