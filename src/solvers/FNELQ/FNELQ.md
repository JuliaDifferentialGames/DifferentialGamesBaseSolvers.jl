# FNELQ
**Reference:** Başar, T., & Olsder, G. J. (1999). *Dynamic Noncooperative Game Theory* (2nd ed.). SIAM. §6.17
**Type:** Discrete-Time, Non-Zero-Sum, Feedback Nash Equilibrium
**Description:** Computes the exact, closed-form feedback Nash equilibrium for finite-horizon discrete-time Linear-Quadratic (LQ) games by solving the coupled backward Riccati recursion in a single pass. No iteration required; the solution is exact up to floating-point precision.

---

## Problem Class
- **Game Type:** Linear-Quadratic (LQ)
- **Players:** N-player (general)
- **Information Structure:** Feedback (closed-loop)
- **Dynamics:** Linear (time-invariant or time-varying)
- **Cost/Objective:** Quadratic stage and terminal costs with optional affine terms

---

## Mathematical Formulation

### Dynamics
```
x_{k+1} = A(k) x_k + \sum_{i=1}^N B_i(k) u_{i,k} + w_k
```

### Cost for Player i
```
J_i = \sum_{k=0}^{N-1} [x_k^T Q_i(k) x_k + u_{i,k}^T R_i(k) u_{i,k} + 2 q_i(k)^T x_k + 2 r_i(k)^T u_{i,k} + 2 u_{i,k}^T M_i(k) \sum_{j\neq i} u_{j,k}] + x_N^T Qf_i x_N + 2 qf_i^T x_N
```

### Assumptions
- A(k), B_i(k), Q_i(k), R_i(k), M_i(k) are time-varying matrices (LTV case) or time-invariant (LTI case)
- Q_i(k) are positive semi-definite
- R_i(k) are positive definite
- S matrix (coupled Riccati gain system) is non-singular at all time steps

---

## Algorithm

### Core Idea
Solves coupled Riccati equations backward in time for feedback Nash equilibrium. The solution is obtained by solving a single linear system at each time step, making it computationally efficient with O(N) complexity where N is the time horizon.

### Steps
1. **Initialization:**
   ```
   Z_i(N) = Qf_i, ζ_i(N) = qf_i for all players i
   ```

2. **Backward Recursion (k = N-1, ..., 0):**
   - **Step 2a:** Assemble the coupled gain matrix S and right-hand side matrices YP, Yα
   - **Step 2b:** Solve S · P_full = YP for feedback gains P_i(k)
   - **Step 2c:** If affine terms present, solve S · α_full = Yα for affine terms α_i(k)
   - **Step 2d:** Update cost-to-go matrices Z_i(k) and ζ_i(k)

3. **Feedback Strategy:**
   ```
   u_{i,k} = -P_i(k) x_k - α_i(k)
   ```

### Key Equations
- Coupled gain system: `S[ri, rj] = (i == j ? R_i : 0) + B_i^T Z_i B_j`
- Feedback gain: `P_full = S \ YP`
- Affine term: `α_full = S \ Yα`
- Cost-to-go update: `Z_i(k) = F^T Z_i(k+1) F + Q_i + P_i^T R_i P_i + cross terms`
- Closed-loop matrix: `F = A - B P_full`

---

## Implementation Details

### Inputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `A` | Matrix/Vector{Matrix} | State dynamics matrix | n×n (or n×n×N for LTV) |
| `B_i` | Matrix/Vector{Matrix} | Control input matrix for player i | n×m_i (or n×m_i×N for LTV) |
| `Q_i` | Matrix/Vector{Matrix} | State cost matrix for player i | n×n (or n×n×N for LTV) |
| `R_i` | Matrix/Vector{Matrix} | Control cost matrix for player i | m_i×m_i (or m_i×m_i×N for LTV) |
| `M_i` | Matrix/Vector{Matrix} | Cross-control cost matrix | m_i×m_j (or m_i×m_j×N for LTV) |
| `Qf_i` | Matrix | Terminal state cost matrix | n×n |
| `qf_i` | Vector | Terminal affine term | n |

### Outputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `P_i(k)` | Matrix | Feedback gain for player i at time k | m_i×n |
| `α_i(k)` | Vector | Affine term for player i at time k | m_i |
| `Z_i(k)` | Matrix | Cost-to-go matrix at time k | n×n |
| `ζ_i(k)` | Vector | Affine cost-to-go term at time k | n |

### Fields/Parameters
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `check_singularity` | Bool | `true` | Check S matrix condition number |
| `rcond_threshold` | Float64 | `1e-10` | Warn if `rcond(S)` below this threshold |
| `regularization` | Float64 | `0.0` | Tikhonov regularization added to S matrix |

---

## Theoretical Notes
- **Existence/Uniqueness:** Guaranteed if S is non-singular at all time steps (Theorem 6.17, Başar & Olsder 1999)
- **Convergence:** Backward recursion converges in exactly N steps (finite horizon)
- **Complexity:** O(N · (Σm_i)^3) per time step for N players with total control dimension Σm_i
- **Numerical Stability:** More stable than continuous-time Riccati integration for discrete-time systems

---

## Dependencies
- `LinearAlgebra`
- `DifferentialGamesBase`

---

## Example Usage
```julia
using DifferentialGames

# Define game parameters
A = rand(4, 4)
B = [rand(4, 2) for _ in 1:2]  # 2 players, 2 controls each
Q = [rand(4, 4) for _ in 1:2]
R = [rand(2, 2) for _ in 1:2]
Qf = [rand(4, 4) for _ in 1:2]

# Create LQ game problem
game = LQGameProblem(A, B, Q, R, Qf, initial_state=rand(4), horizon=10)

# Solve with FNELQ
solver = FNELQ(check_singularity=true, rcond_threshold=1e-10)
sol = solve(game, solver)

# Extract feedback gains for player 1 at time k=0
P_1_0 = sol.strategy.P[1][1]  # Player 1, time step 0
α_1_0 = sol.strategy.α[1][1]  # Affine term for player 1, time step 0
```

---

## References
- Başar, T., & Olsder, G. J. (1999). *Dynamic Noncooperative Game Theory* (2nd ed.). SIAM.
- Engwerda, J. C. (2005). LQ dynamic optimization and differential games. John Wiley & Sons.

---

## Validation
- **Test Cases:** Validated against analytical solutions for 2-player LQ games
- **Benchmarks:** Compared with open-loop solutions for consistency
- **Numerical Tests:** Includes tests for time-invariant and time-varying systems, with and without affine terms

---

## Performance Tips
- Pre-allocate matrices for large N to reduce memory overhead
- Use sparse matrices if A or B are sparse (though FNELQ currently uses dense algebra)
- For ill-conditioned problems, increase `regularization` parameter
- For large-scale problems, consider using the LTV version with pre-allocated arrays

---

## Limitations
- Assumes linear dynamics and quadratic costs
- Not suitable for games with non-convex costs or nonlinear dynamics
- Requires S matrix to be non-singular (may fail for certain problem instances)
- No constraint handling (unconstrained game only)