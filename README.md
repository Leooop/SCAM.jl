# SCAM 

A [Julia](http://julialang.org) 0-dimensional implementation of the Sub-Criticaly Altered Maxwell model.

#### Author
- Léo Petit, École Normale Supérieure de Paris, France.

#### License

`SCAM` is licensed under the [MIT license](./LICENSE.md).

#### Installation

`SCAM` is not a registered package and can be installed from the package REPL with
```julia
pkg> add https://github.com/Leooop/SCAM.jl.git
```
or similarly with
```julia
julia> using Pkg ; Pkg.add("https://github.com/Leooop/SCAM.jl.git")
```
Requires Julia v1.7 or higher

#### Introduction

This rheological model considers the growth of tensile cracks in a compressive state of stress,based on the wing-crack model of Ashby & Sammis (1991) coupled with a sub-critical crack growth law (Charles, 1958).
The material is assumed elastically incompressible, in view of a use in long-term tectonic simulations.
It uses an empirical linear dependance of the shear modulus on damage $D = \frac{4}{3}\pi N_v (l+\alpha a)^3$, were $N_v$ is the number of cracks per unit volume, $l$ is the length of each tensile crack growing from the tips of closed penny-shaped cracks of radius $a$ and oriented at an angle $\psi=\cos^{-1}{\alpha}$.
The penny-shaped cracks normals are assumed to be be contained in the $\sigma_1$--$\sigma_3$ plane, such that long term behavior post crack coalescence (at $D \sim  1$) can be represented by 2-dimentional Mohr-Coulomb plasticity in the same plane.

This model in the 0-D approximation is able to accurately describe the deformation of compact rocks under various confining pressures, strains rate and under creep conditions (stress kept constant). The dependence of temperature is, for now, not implemented.

This model can be coupled with the unregistered packages [DataFormatter.jl](https://github.com/Leooop/DataFormatter.jl) and [ParametersEstimatator.jl](https://github.com/Leooop/ParametersEstimator.jl) to perform bayesian parameters inversion against triaxial experimental data under constant strain rate or brittle creep conditions.

#### Usage



#### Examples


