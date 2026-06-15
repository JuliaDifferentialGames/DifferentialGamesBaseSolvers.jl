# InverseLQGames
**Reference:** Inga, J., Bischoff, E., Molloy, T.L., Flad, M., Hohmann, S. (2019). Solution sets for inverse non-cooperative linear-quadratic differential games. IEEE Control Systems Letters, 3(4), 871-876. DOI: 10.1109/LCSYS.2019.2919271
**Type:** Continuous-Time, Inverse Problem, Linear-Quadratic, Infinite Horizon
**Description:** Algebraic solver for inverse infinite-horizon LQ differential games. Computes the canonical parameter set characterizing all cost-function matrices consistent with an observed Nash equilibrium. Given feedback gains K* = (K₁*,...,Kₙ*), finds ALL (Qᵢ, Rᵢ₁,...,RᵢN) consistent with Nash.

---

## Problem Class
- **Game Type:** Inverse Linear-Quadratic (LQ) Differential Game
- **Players:** N-player (general)
- **Information Structure:** Feedback (closed-loop)
- **Dynamics:** Linear time-invariant
- **Cost/Objective:** Quadratic costs to be inferred
- **Horizon:** Infinite

---

## Mathematical Formulation

### Dynamics
```
ẋ(t) = A x(t) + Σᵢ Bᵢ uᵢ(t)
```

### Control (Feedback Nash Equilibrium)
```
uᵢ(t) = −Kᵢ x(t)   (linear feedback Nash equilibrium)
```

### Cost for Player i (to be inferred)
```
Jᵢ = ½ ∫₀^∞ (xᵀ Qᵢ x + Σⱼ uⱼᵀ Rᵢⱼ uⱼ) dt
```

### Assumptions
- System is stabilizable (closed-loop matrix F is Hurwitz)
- Observed feedback gains K* constitute a Nash equilibrium
- Rᵢᵢ ≻ 0 (positive definite control costs)

---

## Algorithm

### Core Idea
Characterizes all cost matrices (Qᵢ, Rᵢ₁,...,RᵢN) that are consistent with an observed Nash equilibrium by computing the null space of specially constructed matrices Mᵢ. The solution set Θ = ∩ᵢ ker(Mᵢ) with Rᵢᵢ ≻ 0 gives all possible cost parameters.

### Steps (Exact K* Mode)
1. **Compute Closed-Loop Matrix:**
   ```
   F = A − Σᵢ Bᵢ Kᵢ*  (eq. 3)
   ```

2. **Form Kronecker Sum:**
   ```
   F⊕ = Fᵀ ⊗ Iₙ + Iₙ ⊗ Fᵀ  (eq. 9)
   ```

3. **Compute Sᵢ Matrices:**
   ```
   Sᵢ = (Iₙ ⊗ Bᵢᵀ) F⊕⁻¹  (eq. 10)
   ```

4. **Form Mᵢ Matrices:**
   ```
   Mᵢ = [Sᵢ | Sᵢ K₁^⊗ | ⋯ | (Sᵢ Kᵢ^⊗ + Kᵢᵀ ⊗ I_{pᵢ}) | ⋯ | Sᵢ Kₙ^⊗]  (eq. 14)
   ```
   where Kᵢ^⊗ = Kᵢᵀ ⊗ Kᵢᵀ

5. **Compute Null Space:**
   ```
   θᵢ ∈ ker(Mᵢ)  (eq. 13)
   Θ = ∪ᵢ ker(Mᵢ) with Rᵢᵢ ≻ 0  (eq. 19)
   ```

### Steps (Trajectory Mode)
0. **Estimate K̂ᵢ from Observed Trajectories:**
   ```
   K̂ᵢ = argmin Σₖ ‖Kᵢ x^[k] + uᵢ^[k]‖₂²  (eq. 22)
   K̂ᵢ = −Uᵢ Xᵀ (X Xᵀ + λI)⁻¹, λ = ridge
   ```
Then proceed as above with K̂ᵢ in place of Kᵢ*.

---

## Implementation Details

### Inputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `A` | Matrix | State dynamics matrix | n×n |
| `B_i` | Vector{Matrix} | Control input matrices | n×m_i for each player i |
| `K_i` | Vector{Matrix} | Observed feedback gains | m_i×n for each player i |
| `X` | Matrix | State trajectory samples | n×K_samples |
| `U_i` | Vector{Matrix} | Control trajectory samples | m_i×K_samples for each player i |

### Outputs
| Variable | Type | Description | Shape/Dimensions |
|----------|------|-------------|------------------|
| `Q_i` | Vector{Matrix} | Inferred state cost matrices | n×n for each player i |
| `R_ij` | Vector{Vector{Matrix}} | Inferred control cost matrices | m_j×m_j for each player i,j |
| `null_space` | Vector{Matrix} | Null space basis for each player | Depends on problem |
| `converged` | Bool | Whether solution satisfies tolerance | Scalar |

### Fields/Parameters
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `tol` | Float64 | 1e-6 | Residual threshold; converged is true when all ‖Mᵢθᵢ‖ ≤ tol |
| `K_ridge` | Float64 | 0.0 | Tikhonov regularization λ for K estimation; set > 0 for noisy trajectory data |
| `svd_tol` | Float64 | 1e-8 | Relative singular-value threshold for null-space computation |

---

## Theoretical Notes
- **Existence/Uniqueness:** Solution set Θ = ∩ᵢ ker(Mᵢ) with Rᵢᵢ ≻ 0 characterizes all cost matrices consistent with observed Nash equilibrium (Theorem 2, Inga et al. 2019)
- **Convergence:** Exact algebraic solution (no iteration required for exact K* mode)
- **Complexity:** O(n³ + Σ(pᵢ³)) for n states and pᵢ controls per player
- **Numerical Stability:** SVD-based null space computation with relative tolerance; ridge regularization for trajectory estimation

---

## Dependencies
- `LinearAlgebra`
- `DifferentialGamesBase`

---

## Example Usage
```julia
using DifferentialGames

# Exact K* mode: given observed feedback gains
A = [0 1; -1 0]  # 2D system
B = [[1 0], [0 1]]  # 2 players, each with 1D control
K_star = [[1 2], [3 1]]  # Observed feedback gains

# Create inverse LQ game problem
problem = InverseLQGame(A, B, K_star)

# Solve
solver = InverseLQGames(tol=1e-6, svd_tol=1e-8)
sol = solve(problem, solver)

# Extract inferred cost matrices
Q_inferred = sol.Q  # Q_i for each player
R_inferred = sol.R  # R_ij for each player i, control j

# Trajectory mode: estimate K from observed data
X_samples = rand(2, 100)  # 100 state samples
U_samples = [rand(1, 100), rand(1, 100)]  # Control samples for 2 players

problem_trajectory = InverseLQGame(A, B, X_samples, U_samples)
sol_trajectory = solve(problem_trajectory, InverseLQGames(K_ridge=0.1))
```

---

## References
- Inga, J., Bischoff, E., Molloy, T.L., Flad, M., Hohmann, S. (2019). Solution sets for inverse non-cooperative linear-quadratic differential games. IEEE Control Systems Letters, 3(4), 871-876. DOI: 10.1109/LCSYS.2019.2919271

---

## Validation
- **Test Cases:** Validated on known LQ games with analytical solutions
- **Benchmarks:** Verified that inferred costs produce the observed equilibrium
- **Numerical Tests:** Includes tests for both exact K* and trajectory estimation modes

---

## Performance Tips
- Use exact K* mode when feedback gains are known precisely
- Increase K_ridge for noisy trajectory data
- Adjust svd_tol for problems with nearly singular Mᵢ matrices
- For large systems, consider sparse matrix representations

---

## Limitations
- Infinite horizon only (no finite horizon support)
- Linear dynamics only
- Requires observed Nash equilibrium (either as K* or trajectory data)
- Solution set may be empty if observed gains are not a Nash equilibrium
- No handling of state or control constraints