# Topology Optimization

**Topology optimization of a cantilever beam in Wolfram Mathematica (SIMP + optimality criteria), with spanwise bending-stiffness evaluation, an FE load test with stress visualization, and STL export of the optimized design.**

A compact, self-checking Mathematica notebook for minimum-compliance topology optimization on a 5:2 cantilever. It follows the classical density-based approach (SIMP penalization, sensitivity filter, optimality-criteria update) and adds post-processing that is useful for a beam-type structure: the cross-sectional bending stiffness $EI(x)$ along the span, a stress plot of the final design, and a watertight 3D mesh for CAD or 3D printing.

Only long-standing core functions are used.

## Problem

Minimize the compliance of a clamped beam under a vertical point load, subject to a volume constraint:

```
min   c(ρ) = Fᵀ u
s.t.  K(ρ) u = F
      Σ ρₑ / N ≤ V*          (volume fraction, default V* = 0.4)
      0.001 ≤ ρₑ ≤ 1
```

| Item | Setting |
|---|---|
| Domain | Rectangle, aspect ratio 5:2, default 100 × 40 elements |
| Elements | Bilinear Q4, plane stress, unit thickness, element size 1 |
| Supports | Left edge fully clamped |
| Load | Vertical point load at mid-height of the right edge (split over two nodes if the number of element rows is odd) |
| Material | E₀ = 1, Eₘᵢₙ = 10⁻⁹, ν = 0.3 (dimensionless) |

## Method

- **SIMP interpolation:** Eₑ = Eₘᵢₙ + ρₑᵖ (E₀ − Eₘᵢₙ), p = 3.
- **Sensitivities:** ∂c/∂ρₑ = −p ρₑ^(p−1) (E₀ − Eₘᵢₙ) uₑᵀ Kₑ uₑ.
- **Sensitivity filter** (mesh-independence): weights Hₑf = max(0, r_min − dist(e, f)), default r_min = 1.6 elements, built once as a sparse matrix.
- **Optimality-criteria update** with bisection on the Lagrange multiplier, move limit 0.2, lower bound 0.001.
- **Stopping criterion:** max |Δρ| < 0.01 after more than 15 iterations, or `maxIter`.


<img src="./pics/convergence.jpg" alt="convergence"  width="80%" />

<img src="./pics/topology.jpg" alt="topology"  />

<img src="./pics/stiffness.jpg" alt="local stiffness" width="80%"  />

<img src="./pics/fem.jpg" alt="FE results"  />


## Reference results

Default settings with a unit load (the compliance scales with P²):

| Quantity | Value |
|---|---|
| Compliance, uniform start ρ = 0.4 (iteration 1) | ≈ 1119 |
| Compliance at convergence (57 iterations) | ≈ 139.8 |
| Final volume fraction | ≈ 0.4000 |
| Section equilibrium M_num / M_ref | mean 0.9996, range 0.9979 to 1.0007 |

With the default `maxIter = 30` the stopping criterion is not yet met (max |Δρ| ≈ 0.034); set `maxIter` to 60 or more to run to convergence.

## Limitations

- 2D, single load case, compliance minimization with a volume constraint only.
- Only a sensitivity filter is used: no projection filter, so intermediate densities remain; no robust formulation; no stress constraint (stresses are evaluated after the optimization only).
- All quantities are dimensionless (E₀ = 1, element size 1).
- The σₓₓ field uses the penalized (SIMP) element modulus.

## References

- O. Sigmund, *A 99 line topology optimization code written in Matlab*, Structural and Multidisciplinary Optimization 21 (2001), 120–127.
- E. Andreassen, A. Clausen, M. Schevenels, B. S. Lazarov, O. Sigmund, *Efficient topology optimization in MATLAB using 88 lines of code*, Structural and Multidisciplinary Optimization 43 (2011), 1–16.
- M. P. Bendsøe, O. Sigmund, *Topology Optimization: Theory, Methods and Applications*, Springer, 2003.
