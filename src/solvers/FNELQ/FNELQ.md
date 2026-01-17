# FNELQ
Based on Başar & Olsder (1999), Theorem 6.17 (discrete-time case). Computes 
feedback Nash equilibrium by solving coupled matrix equations backward in time.

## Mathematical Formulation
Dynamics: xₖ₊₁ = A xₖ + Σᵢ Bᵢ uᵢₖ

Cost for player i:
    Jᵢ = Σₖ [xₖᵀQᵢxₖ + uᵢₖᵀRᵢuᵢₖ + 2qᵢᵀxₖ + 2rᵢᵀuᵢₖ] + xₙᵀQfᵢxₙ

## Algorithm (Başar & Olsder, Section 6.3)
Dynamic programming backward recursion:

1. Initialize: Zᵢ(N) = Qfᵢ, ζᵢ(N) = 0

2. For k = N-1, ..., 0:
   a. Solve coupled linear system for feedback gains:
      S·P = Y where:
      - S[uⁱ, uʲ] = Rⁱ·δᵢⱼ + BᵢᵀZᵢBⱼ  (block structure)
      - Y[uⁱ, :] = BᵢᵀZᵢA
   
   b. Compute affine term: α = S⁻¹·Yα where Yα[uⁱ] = Bᵢᵀζᵢ + rᵢ
   
   c. Update cost-to-go:
      F = A - B·P
      β = -B·α
      Zᵢ(k) = FᵀZᵢ(k+1)F + Qᵢ + PᵀRᵢP
      ζᵢ(k) = Fᵀ(ζᵢ(k+1) + Zᵢ(k+1)β) + qᵢ + PᵀRᵢα - Pᵀrᵢ

3. Feedback strategy: uᵢₖ = -Pᵢ(k)xₖ - αᵢ(k)

## Key Difference from Open-Loop
- Solves linear matrix equation (no iteration needed at each time step)
- Coupling through S matrix construction (all players simultaneously)
- Guaranteed unique solution if S is non-singular

## Fields
- `check_singularity::Bool` : Check S matrix condition number (default: true)
- `rcond_threshold::Float64` : Warn if rcond(S) below this (default: 1e-10)

## Theoretical Notes
- Existence requires S non-singular at all time steps (Theorem 6.17)
- S singularity indicates lack of unique Nash equilibrium
- Scales as O(N²·m²) per time step for N players with m controls each
- More numerically stable than continuous-time Riccati integration