# SFNELQ
**Reference:** Başar, T., & Olsder, G. J. (1999). *Dynamic Noncooperative Game Theory* (2nd ed.). SIAM. §6 (feedback Nash for discrete-time stochastic LQ games).
**Type:** Discrete-Time, Non-Zero-Sum, Stochastic Feedback Nash Equilibrium
**Description:** Computes the risk-neutral feedback Nash equilibrium of a discrete-time stochastic LQ game whose process noise is additive and Gaussian. Delegates gain computation to `FNELQ` under certainty equivalence and reports the expected cost correction attributable to the noise.

---

## Problem Class
- **Game Type:** Linear-Quadratic (LQ) with additive Gaussian process noise
- **Players:** N-player (general)
- **Information Structure:** Feedback (closed-loop)
- **Dynamics:** `StochasticGameProblem` wrapping `LinearDynamics`, with an `ItoSDEDynamics`/`ConstantDiffusion` additive-noise term (built via `ProcessNoiseModel`)
- **Cost/Objective:** Quadratic stage and terminal costs with optional affine terms

---

## Mathematical Formulation

### Dynamics
```
x_{k+1} = A(k) x_k + \sum_{i=1}^N B_i(k) u_{i,k} + w_k,   w_k ~ N(0, W)
```

### Cost for Player i
Identical to FNELQ's — see [FNELQ](@ref) — but evaluated in EXPECTATION over `w_k`:
```
J_i = E[ \sum_{k=0}^{N-1} ( ... FNELQ's stage cost ... ) + ... terminal cost ... ]
```

### Certainty Equivalence
Because the noise is additive (not state- or control-multiplicative), it enters only the
value function's constant term, never its quadratic (`P`) or linear (`p`/feedforward)
parts. Consequently the equilibrium feedback gains are IDENTICAL to the deterministic
`FNELQ` gains computed on the same nominal game; the noise only shifts each player's
expected cost upward by a computable constant.

---

## Algorithm

### Core Idea
1. Solve the nominal (noise-free) game with `FNELQ`, obtaining the feedback gains and the
   cost-to-go sequence `Z[i][k]` (FNELQ's `solver_info[:cost_to_go_matrices]`).
2. Add the additive-noise correction to each player's expected cost:
   ```
   noise_cost[i] = Σ_{k=1}^{N} ½ tr(Z[i][k+1] · W)
   ```
3. Return the same feedback strategy as FNELQ, with trajectory costs and `solver_info`
   updated to reflect the noise correction.

### Key Equations
- Noise correction: `noise_cost[i] = Σ_k ½ tr(P_{k+1}^i W)`
- Gains: identical to FNELQ's `K[i][k]` — see [FNELQ](@ref) for the full Riccati recursion

---

## Implementation Details

### Inputs
Same as `FNELQ`, plus a `StochasticGameProblem` wrapping the nominal `GameProblem` with an
additive-noise `ItoSDEDynamics` (built via `ProcessNoiseModel(Σ)`).

### Outputs
| Variable | Type | Description |
|----------|------|-------------|
| `strategy.gains` | `Vector{Vector{Matrix}}` | Identical to the nominal FNELQ gains |
| `solver_info[:noise_cost]` | `Vector{Float64}` | Per-player expected cost from additive noise |
| `solver_info[:certainty_equivalent]` | `Bool` | Always `true` in this round's scope |
| `solver_info[:expected_costs]` | `Vector{Float64}` | Deterministic cost + `noise_cost` |

### Fields/Parameters
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `check_singularity` | Bool | `true` | Forwarded to the inner `FNELQ` solve |
| `rcond_threshold` | Float64 | `1e-10` | Forwarded to the inner `FNELQ` solve |
| `regularization` | Float64 | `0.0` | Forwarded to the inner `FNELQ` solve |

---

## Theoretical Notes
- **Scope of this port:** state-/control-multiplicative noise (where certainty equivalence
  does *not* hold) is not supported — `DifferentialGamesBase.jl`'s `StochasticGameProblem`
  only models additive noise today. `solve(sgame, SFNELQ())` throws a clear error if
  `is_certainty_equivalent(sgame) == false`. See
  [DifferentialGamesBase.jl's stochastic docs](https://github.com/JuliaDifferentialGames/DifferentialGamesBase.jl)
  for the roadmap.
- **Existence/Uniqueness:** Same as FNELQ (Theorem 6.17, Başar & Olsder 1999) — the noise
  correction never affects feasibility.

---

## Dependencies
- `LinearAlgebra`
- `DifferentialGamesBase`
- `FNELQ` (this package)

---

## Example Usage
```julia
using DifferentialGames

A = [1.0 1.0; 0.0 1.0]
B = [reshape([0.0, 1.0], 2, 1)]
Q = diagm([1.0, 1.0]); R = reshape([0.1], 1, 1); Qf = diagm([10.0, 10.0])
game = LQGameProblem(A, B, [Q], [R], [Qf], [1.0, 0.0], 2.0; dt = 1.0)

Σ = diagm([0.01, 0.01])
sgame = StochasticGameProblem(game, ProcessNoiseModel(Σ))

sol = solve(sgame, SFNELQ())
sol.solver_info[:noise_cost]            # per-player expected cost from noise
sol.solver_info[:certainty_equivalent]  # true
```

---

## References
- Başar, T., & Olsder, G. J. (1999). *Dynamic Noncooperative Game Theory* (2nd ed.). SIAM.

---

## Validation
- **Test Cases:** `test/sfnelq_tests.jl` — gains and mean trajectory identical to
  deterministic `FNELQ` under additive noise; expected cost equals the deterministic cost
  plus the reported `noise_cost`; rejects a multiplicative-noise `StochasticGameProblem`.

---

## Limitations
- Additive noise only — no state-/control-multiplicative noise support this round.
- Requires an LQ nominal game (same requirement as `FNELQ`).
- No constraint handling (unconstrained game only, inherited from `FNELQ`).
