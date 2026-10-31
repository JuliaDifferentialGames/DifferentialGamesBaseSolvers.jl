# DGSQP
**References:** Zhu & Borrelli (2023), *A Sequential Quadratic Programming Approach to the Solution of Open-Loop Generalized Nash Equilibria*, ICRA, https://arxiv.org/abs/2203.16478. Zhu & Borrelli (2024), *… for Autonomous Racing*, https://arxiv.org/abs/2404.00186. Reference implementation: https://github.com/zhu-edward/DGSQP.
**Type:** Discrete-Time, Nonlinear, Constrained, Open-Loop (Generalized) Nash Equilibrium
**Description:** Finds a local open-loop generalized Nash equilibrium (normalized, shared multipliers) by solving one convex QP per iteration on the KKT system of the game, with a merit function and a non-monotone step-acceptance scheme. Two step-acceptance variants: the ICRA 2023 solver (`DGSQP()`, watchdog line search) and the racing paper's solver (`dgsqp_racing_options()`: decaying regularization, d-steps and m-steps).

---

## Problem Class
- **Game type:** `GameProblem` (`LQGameProblem`, `PDGNEProblem`, any private/shared constraints), `OpenLoopGame` (cost and constraints as functions of the control sequences) or `SeparableRolloutGame` (separable dynamics with node-wise coupling and exact structured derivatives)
- **Strategy:** open-loop (`OpenLoopStrategy`, `equilibrium_type = :OpenLoopNash`)

## Usage
```julia
sol = solve(game, DGSQP())                         # ICRA 2023 defaults
sol = solve(game, dgsqp_paper_options())           # the paper's experiment settings
res = dgsqp_solve(openloop_or_rollout_game, w0; solver = dgsqp_racing_options())
```
Full documentation, the options, the replication of both papers and the comparison with FALCON are in the manual page *Solvers → DGSQP*.

## Capabilities
`:LQGame :NonlinearGame :GNEP :ConstrainedGame :UnconstrainedGame :DiscreteTime :SeparableDynamics :SharedConstraints`

## Tests
`test/dgsqp_tests.jl` (unit, analytic, LQ cross-check against ALGAMES, best-response certificates, structured derivatives, paper study) and `test/austin_benchmark_tests.jl` (needs DifferentialGamesBenchmarks). `test/paper_study.jl` reproduces Table I of the ICRA paper (run on a schedule by `.github/workflows/DGSQP.yml`).

## Known limitations
Local method; the "elastic mode" QP of the reference (for infeasible QPs) and BFGS Hessians are not ported; the exact Hessian of a generic `GameProblem` uses nested forward-mode differentiation and is slow for long horizons (use `SeparableRolloutGame`).
