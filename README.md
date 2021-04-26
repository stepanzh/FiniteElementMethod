# FiniteElementMethod
My practice at Finite Element Method.

# Literature

1. *Mats G. Larson, Fredrik Bengzon.*
The Finite Element Method: Theory, Implementation, and Applications.
Springer, 2013.
ISBN 978-3-642-33286-9

2. *Chen, Z., Huan, G., & Ma, Y.*
Computational methods for multiphase flows in porous media.
Society for Industrial and Applied Mathematics, 2006.
ISBN 0-89871-606-3

# Polynomial Approximation
## L2-projection
### One-dimensional
- ![Script](PolynomialApproximation/l2_projection.jl).
- ![Plot of projection](PolynomialApproximation/plots/l2_projection.pdf).
- ![Plot of error](PolynomialApproximation/plots/l2_projection_error.pdf)

### Two-dimensional
- ![Script](PolynomialApproximation/l2_projection_2d.jl).

# Standard FEM

## Stationary elliptic with constant source (Poisson's equation)

- ![Script](StandartFEM/elliptic_stationary_zero_bc.jl).
- ![Plot](StandartFEM/plots/elliptic_stationary_zero_bc_uniform.pdf).

## General stationary parabolic problem

- ![Script](StandartFEM/parabolic_stationary.jl).
- ![Plot](StandartFEM/plots/parabolic_stationary.pdf).
