# DTMPC Theory Audit

Verification status of the design proposal *"Per-Agent Dynamic-Tube MPC Wrapper
Around FALCON"* against the implementation in `src/solvers/DTRHN/`.

Each item is classified:

- **✅ realized** — implemented and exercised by the test suite.
- **🔬 validated-empirically** — the theorem's *observable consequence* is checked
  numerically (the analytical proof is not, and is not claimed here).
- **📝 analytical-open** — a proof obligation that code cannot discharge; tracked, not done.
- **⚠️ partial / inner-solver-limited** — wired correctly but only approximately
  realized with the current inner solver (ALGAMES); flagged for a stronger solver.
- **🚫 not-claimed** — explicitly out of scope (proposal §7/§9).

Test references are to `test/dtmpc_tests.jl` (core) and `test/theory_validation.jl` (D-series).

---

## 1. Core structural result — Proposition 1 (tube-augmented game is still a PD-GNEP)

| Claim | Status | Evidence |
|---|---|---|
| Tube width is a per-agent **scalar** augmented state `sᵢ` (+ `αᵢ`, control `vᵢ`) | ✅ realized | `build_augmented_player_dynamics`; A1/A3 tests |
| Augmented dynamics stay **block-diagonal / separable** ⇒ still PD-GNEP | ✅ realized · 🔬 D2 | `build_augmented_game` ⇒ `is_pd_gnep == true`, `has_separable_dynamics == true` (B1 test, D2) |
| Tube co-optimized **inside** the game (not integrated externally) | ✅ realized | augmented game solved each step; `α`/`v` co-optimized — integration test asserts `α` changes (i.e. `v ≠ 0`) |
| Strong convexity / inner-solver convergence **inherited** | 🔬 partial | inherited from the inner solver (ALGAMES) on the augmented PD-GNEP; we validate convergence empirically, not the constant |

The augmented per-agent **continuous** field (Lopez–Slotine–How Eq. 13) is

```
ẋᵢ = fᵢ(xᵢ,uᵢ) ;   ṡᵢ = −αᵢ sᵢ + Δᵢ(xᵢ) + Dᵢ + ηᵢ ;   α̇ᵢ = vᵢ
```

(`SeparableDynamics` is integrated by the solver, so the tube term is a derivative,
not a discrete Euler map — this differs from the original proposal sketch, which
wrote `πᵢ` as a discrete recursion; the continuous form is the faithful one for
this codebase.) Equation fidelity is asserted in **D1**.

The augmented cost is built as a block-diagonal **LQ** cost
(`Q̃ = blkdiag(Q, diag(wΦ, wα))`, `R̃ = blkdiag(R, wv)`) rather than a generic
nonlinear cost, so the inner solver gets exact analytical derivatives. This also
sidesteps a pre-existing bug in `DifferentialGamesBase`
(`NonlinearStageCost`'s AD-gradient fallback calls an undefined
`automatic_differentiation_gradient`). Requires LTI-LQ original costs and `wv > 0`.

## 2. Decision-dependent constraint tightening — §3, Eq. 7

| Claim | Status | Evidence |
|---|---|---|
| Collision constraint tightened by tube **decision variables**: `d_min − ‖rᵢ−rⱼ‖ + sᵢ + sⱼ ≤ 0` | ✅ realized | `build_tightened_collision_constraint` → `SharedNonlinearInequality` over the joint augmented state; A5 test |
| Tightening is **enforced in-solver** (not a post-hoc check) | ✅ realized · 🔬 D4 | constraint is added to the augmented game's `shared_constraints`; D4 shows tightening-on yields strictly larger separation than tightening-off |
| Tube **allocation is part of the equilibrium** (margin to the cheaper agent) | ⚠️ partial · 🔬 D6 | the cost-asymmetry form is reproducible (cheaper `w_Φ` ⇒ larger tube at horizon end; D6). The *full* coupled-allocation-through-the-shared-constraint is only approximately realized — see inner-solver caveat below. |

**Inner-solver caveat (ALGAMES).** ALGAMES's augmented-Lagrangian enforcement of
shared inequality constraints is *approximate* in the regimes probed here: it can
report `converged` with a residual collision-constraint violation, and it does not
reliably drive the *tube states* to satisfy a constraint that can only be met by
tube shrinkage (the tube is a state reachable only indirectly through `v → α → s`).
Consequently the sharp §3 allocation experiment (positions frozen, constraint
satisfiable only by trading tube widths) does **not** converge cleanly with
ALGAMES. The mechanism is correctly *wired* (A5 builds the right decision-dependent
constraint; the inner solver consumes it; D4 confirms the margin acts); a stronger
PD-GNEP solver (the external **FALCON**, the proposal's intended base) is the path
to the exact allocation. The solver is abstracted (`DynamicTubeMPC{<:GameSolver}`),
so swapping it in requires no DTMPC changes.

## 3. Separability tripwire — §4 / obligation **T-1**

| Item | Status | Evidence |
|---|---|---|
| Relative-sensing `Δᵢ(xᵢ, x₋ᵢ)` breaks separability unless neighbor is frozen | ✅ realized | `build_augmented_game` freezes the neighbor reference at the current outer-iterate nominal; coupled uncertainty models supported via `create_coupled_uncertainty_models` |
| Frozen-neighbor tube is a **conservative over-approximation** | 📝 analytical-open (**T-1**) | the *construction* preserves separability (D2 checks block-diagonality); the conservativeness bound is a standard sensitivity argument not proved here |

## 4. Robustness vs. suboptimality — Proposition 2 / obligation **R-1**

| Item | Status | Evidence |
|---|---|---|
| Executed profile is an `ε`-OLNE with `ε = O(diam Dᵢ)` | 🔬 D5 (observable) · 📝 R-1 (analytical-open) | D5 shows the tube margin scales monotonically with the disturbance bound `D`; the Lipschitz/gap-function bound with a **decision-dependent** margin through the AL loop (the genuinely new wrinkle vs. E2) is not proved |

## 5. Recursive feasibility & stability — obligations **F-1**, **F-2**

| Item | Status |
|---|---|
| Per-agent tube invariance / ISS (boundary-layer) | 🔬 D3 (observable): Monte-Carlo `‖xᵢ − x̄ᵢ‖ ≤ sᵢ` under sampled disturbances with a designed ancillary gain |
| **F-1** joint robust terminal ingredients for the game | 📝 analytical-open — not implemented (no robust terminal set synthesis) |
| **F-2** disturbed gap-function decrease ⇒ practical tracking | 📝 analytical-open — not implemented |

The ancillary boundary-layer controller (`per_agent_boundary_layer_sliding_control`,
Eq. 2 / Lopez Eq. 9) is applied at *execution* time and is driven by the augmented
solution's scalar `sᵢ`. It requires a user-supplied per-agent gain `Kᵢ ∈ ℝ^{mᵢ×nᵢ}`;
with `nothing` the nominal control is applied unchanged (no correction). The
closed-loop RH solve uses a deterministic (default zero) disturbance so it is
reproducible; disturbance realizations for invariance testing are injected by the
D3 harness, not the main loop.

## 6. Equilibrium-semantics caveats — §7 (do **not** oversell)

| Claim | Status |
|---|---|
| Receding-horizon Nash *trajectory*, **not** a feedback / Markov-perfect NE | 🚫 not-claimed |
| Open-loop NE are time-inconsistent (re-commit each step) | 🚫 not-claimed |
| Tube robustness ≠ strategic robustness (tube is for *nature*, not a rational opponent) | 🚫 not-claimed — strategic robustness is the separate E5 robust-game variant |

## 7. Proof-obligation ledger (proposal §9)

| ID | Statement | Status here |
|---|---|---|
| T-1 | Frozen-neighbor tube is a conservative over-approximation | construction ✅ realized; bound 📝 analytical-open |
| R-1 | `ε`-OLNE bound with decision-dependent margin through the AL loop | observable 🔬 D5; proof 📝 analytical-open |
| F-1 | Joint robust terminal set under tightened coupling (game) | 📝 analytical-open |
| F-2 | Disturbed gap-function decrease ⇒ practical tracking | 📝 analytical-open |

## 8. Citation-verification list (status)

The following are cited in the proposal with `[verify]` tags. They are **not**
independently re-verified in this code work; they remain documentation TODOs:
Lopez–Slotine–How (2019) venue/pages; Slotine & Li (1991); Mayne–Seron–Raković
(2005); Raković et al. (2005); Richards & How (2007); Conte et al. (2016);
Köhler–Zeilinger–Grüne (2023) & Grüne et al. (2010); Başar & Olsder (1999) §6;
Fiacco (1983) / Robinson (1980). FALCON (Outland & Arya 2026) is the intended
PD-GNEP base solver and is **external** to this repository (`methods/FALCON/`).

---

## Summary

The **load-bearing structural claim (Proposition 1)** is fully realized: the
per-agent scalar tube is an augmented PD-GNEP state co-optimized inside the game,
and the tube/bandwidth genuinely respond to the optimization (D1–D2, integration
tests). The **decision-dependent tightening (§3, Eq. 7)** is correctly built and
enforced in-solver, with its first-order effect confirmed (D4) and the
authority-cost allocation demonstrated in its reproducible form (D6). The
**observable consequences** of the robustness theorems are validated numerically
(D3 invariance, D5 margin scaling). The remaining items (T-1 bound, R-1, F-1, F-2)
are analytical obligations left open, and the equilibrium-semantics limits of §7
are not overclaimed. The one place reality falls short of the proposal is the
**exact coupled tube-allocation equilibrium**, limited by ALGAMES's approximate
shared-constraint handling — the motivation for keeping the inner solver
abstracted and reserving the exact result for FALCON.
