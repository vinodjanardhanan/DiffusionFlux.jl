```@meta
CurrentModule = DiffusionFlux
```

# DiffusionFlux
DiffusionFlux is a package for the calculation of diffusion flux. The packge implements different models for the calculation
of diffusion fluxes.

Documentation for [DiffusionFlux](https://github.com/vinodjanardhanan/DiffusionFlux.jl).

## Installation
To install the package, use the following commands in the julia REPL
```julia
julia> using Pkg
julia> Pkg.add("DiffusionFlux")
```


## Models
### Fickian diffusion
```math
j_k = -D_{km} M_k \frac{dc_k}{dx}
```
In the above eqaution $j_k$ is the mass flux (kg/m$^2$s), $D_{km}$ is the diffusion coefficient of species $k$ in the rest of the mixture (m$^2$/s), $M_k$ is the molecular weight (kg/mol), $c_k$ is the concentration of species $k$ (mol/m$^3$), and $x$ is the spatial coordinate.  The flux calculated using the above formulation does not ensure that $\sum j_k = 0$. Therefore the fluxes are corrected according to 
```math
j_k^{\prime} = j_k - Y_k j_{corr}
```
where
```math
j_{corr} = \sum j_k
```

Two different formulations are available for the calculation of species flux through porous media
### Modified Fickian (porous media flow)
```math
j_k = M_k \left(-D^e_{k} \frac{dc_k}{dx} - \frac{B_g c_k}{\mu}\frac{p}{dx}\right)
```
```math
\frac{1}{D^e_{k}} = \frac{1}{D^e_{k,Kn}} + \frac{1}{D^e_{km}}
```
```math
D^e_{k,kn} = \frac{\epsilon}{\tau}D_{k,Kn}
```
```math
D^e_{km} = \frac{\epsilon}{\tau}D_{km}
```
In the above equations, $\epsilon$ is the porosity, $\tau$ is the tortuosity, $B_g$ is the permeability  (m$^2$), $p$ is the pressure (Pa), and $D_{k,Kn}$ is the Knudsen diffusion coefficient of species $k$ (m$^2$/s).
### Dusty Gas Model (porous media flow)

```math
j_k = -M_k\left[ \sum_{l=1}^{N_g} D_{kl}^\mathrm{DGM} \frac{dc_l}{dx} + \left( \sum_{l=1}^{N_g}\frac{D_{kl}^\mathrm{DGM}c_l}{D_{l,Kn}} \right)\frac{B_g}{\mu} \frac{dp}{dx}  \right]
```
```math
D_{kl}^\mathrm{DGM} = H^{-1}
```
```math
h_{kl} = \left[ \frac{1}{D^e_{k,Kn}} + \sum_{j \ne k} \frac{X_j}{D^e_{kj}} \right] \delta_{kl} + (\delta_{kl}-1)\frac{X_k}{D^e_{kl}}
```
In the above equation $X_k$ is the mole fraction of species $k$, $D^e_{kj}$ is the effective binary diffusion coeffcients, and $\delta_{kl}$ is the kronecker delta.
## General interfaces
```@index
```

```@autodocs
Modules = [DiffusionFlux]
```
