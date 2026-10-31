# Dynamic Tube MPC (DTMPC)

**References:** Lopez, Slotine & How (2019), *Dynamic Tube MPC for Nonlinear Systems*;
Outland & Arya (2026), *FALCON* (PD-GNEP base algorithm, external).
**Type:** Discrete-time, robust, tube-based receding-horizon **game** solver.

A per-agent **Dynamic Tube MPC** wrapper that realizes the design proposal
*"Per-Agent Dynamic-Tube MPC Wrapper."* Each agent carries a **scalar tube-width**
state that is **co-optimized inside the game** (Proposition 1), and coupled
collision constraints are tightened **decision-dependently** (Eq. 7) so that tube
allocation becomes part of the equilibrium (§3).

> The verification status of every theoretical claim is tracked in
> [`THEORY_AUDIT.md`](./THEORY_AUDIT.md).

---

## Problem class

- **Horizon game:** must be a **PD-GNEP** (`SeparableDynamics`): each agent `i`
  has its own physical state `xᵢ ∈ ℝ^{nᵢ}`, control `uᵢ ∈ ℝ^{mᵢ}`, continuous
  field `ẋᵢ = fᵢ(xᵢ, uᵢ, p, t)`, and an **LQ** objective.
- **Inner solver:** any constrained `GameSolver` (validated with `ALGAMES`).
- **Uncertainty:** per-agent bound `Δᵢ(xᵢ, x₋ᵢ) → scalar`; relative-sensing
  uncertainty is supported via a frozen-neighbor reference (obligation T-1).

A shared-state game (e.g. `LQGameProblem`) is **rejected** — DTMPC requires
per-agent state blocks.

---

## Augmented formulation (the core idea)

At each receding-horizon step DTMPC builds an **augmented PD-GNEP** and hands it to
the inner solver. Per agent:

```
augmented state    x̃ᵢ = [ xᵢ (nᵢ) ; sᵢ (1) ; αᵢ (1) ]
augmented control  ũᵢ = [ uᵢ (mᵢ) ; vᵢ (1) ]
```

with the augmented **continuous** field (Lopez Eq. 13; the solver integrates it):

```
ẋᵢ = fᵢ(xᵢ, uᵢ)
ṡᵢ = −αᵢ·sᵢ + Δᵢ(xᵢ) + Dᵢ + ηᵢ        # tube width
α̇ᵢ = vᵢ                               # bandwidth
```

augmented quadratic cost (built as a block-diagonal `LQStageCost`/`LQTerminalCost`):

```
ℓ̃ᵢ = ℓᵢ(xᵢ, uᵢ) + wΦ·sᵢ² + wα·αᵢ² + wv·vᵢ²
```

and the decision-dependent tightened collision constraint (shared, Eq. 7):

```
c_ij(x̃) = d_min − ‖rᵢ − rⱼ‖ + sᵢ + sⱼ ≤ 0
```

At **execution**, the true state is rolled forward under the ancillary
boundary-layer sliding control `uᵢ = ūᵢ + κᵢ(xᵢ − x̄ᵢ)` (Eq. 2), driven by the
solved scalar `sᵢ`.

---

## API

```julia
using DifferentialGamesBase, DifferentialGamesBaseSolvers

# 1. A PD-GNEP horizon game: 2 agents, each a 2-D double integrator
di(xi, ui, p, t) = [xi[2], ui[1]]
obj(i) = PlayerObjective(i, DiagonalLQStageCost([1.0, 0.1], [0.5]),
                            DiagonalLQTerminalCost([5.0, 1.0]))
players = [PlayerSpec(i, 2, 1, [Float64(i), 0.0], di, obj(i)) for i in 1:2]
game = PDGNEProblem(players, 0.8, 0.1)          # tf, dt

# 2. Per-agent scalar-tube parameters
params = DTMPCParameters(Float64, 2;
    s₀ = [0.1, 0.1], α₀ = [1.0, 1.0],
    D  = [0.05, 0.05], η = [0.01, 0.01],
    coupled_pairs = [(1, 2)], collision_min_distance = 0.5, pos_dim = 1)

# 3. Per-agent uncertainty models  Δᵢ(xᵢ, x₋ᵢ) -> scalar
um = create_independent_uncertainty_models(x -> 0.05, 2)

# 4. Problem + solve (abstract over the inner solver)
prob = DynamicTubeMPCProblem(game, [1.0,0.0, 2.0,0.0], 6, params, um)
sol  = solve(prob, DynamicTubeMPC(ALGAMES()))
```

### Solution accessors

| Call | Returns |
|---|---|
| `get_state_trajectory(sol)` / `dtmpc_state_trajectory(sol)` | joint physical state `(n × N_sim+1)` |
| `get_control_trajectory(sol, i)` | agent `i` nominal physical controls |
| `get_total_cost(sol, i)` | agent `i` accumulated cost |
| `dtmpc_tube_trajectory(sol, i)` | agent `i` scalar tube-width trajectory `sᵢ` |
| `dtmpc_bandwidth_trajectory(sol, i)` | agent `i` scalar bandwidth trajectory `αᵢ` |
| `all_coupled_constraints_satisfied(sol)` | post-hoc collision-constraint monitor |

### Key builders (reusable / testable)

`build_augmented_player_dynamics`, `build_augmented_cost`,
`build_tightened_collision_constraint`, `build_augmented_game`,
`extract_agent_augmented`, `per_agent_boundary_layer_sliding_control`,
`create_independent_uncertainty_models`, `create_coupled_uncertainty_models`.

---

## Requirements & caveats

- Original costs must be **time-invariant LQ** and `w_v > 0` (keeps the augmented
  control weight `R̃` positive definite).
- The ancillary gain `boundary_layer_gain[i]` is an optional `mᵢ × nᵢ` matrix
  (designed offline); `nothing` ⇒ no ancillary correction.
- The closed-loop solve uses a deterministic (default zero) disturbance for
  reproducibility; pass `disturbance = t -> d` to inject closed-loop disturbances.
- **Inner-solver limitation:** ALGAMES enforces shared inequality constraints only
  approximately, so the exact coupled tube-allocation equilibrium (§3) is not fully
  reproduced; the mechanism is correctly wired and its first-order effect verified
  (see `THEORY_AUDIT.md`, D4/D6). The solver is abstracted to allow swapping in
  FALCON.

## Not claimed (proposal §7)

Feedback / Markov-perfect Nash optimality, time-consistency of the open-loop NE,
and robustness to a *strategic* opponent. The tube hedges against bounded
disturbance (*nature*), not a rational adversary.

## Tests

- `test/dtmpc_tests.jl` — builders (A1–A5), augmented-game assembly (B1),
  end-to-end ALGAMES integration (2- & 3-agent), warm-start parity, accessors.
- `test/theory_validation.jl` — D1 equation fidelity, D2 separability, D3 tube
  invariance (Monte-Carlo), D4 robust constraint satisfaction, D5 margin scaling,
  D6 authority-cost allocation.

Run: `include("src/solvers/DTRHN/test/runtests.jl")` (96 tests).
