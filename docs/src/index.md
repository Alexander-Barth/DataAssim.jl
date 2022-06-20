# DataAssim.jl

Documentation for [DataAssim.jl](https://github.com/Alexander-Barth/DataAssim.jl)

## Simulation driver


The observations `yo` are provided as a vector of vectors (possibly of different length). Similarily their error covariance `R` is a vector of matrices.
The elements of the vectors `yo` and `R` can be constructed as needed by the type `DataAssim.VectorFun`:

```julia
using DataAssim, LinearAlgebra
R = DataAssim.VectorFun(Matrix{Bool},10,n -> (n < 5 : Matrix(I,3,3) : Matrix(I,5,5)))
```

The model ``ℳ`` implementes the following API:

* ``ℳ(n,x)`` apply the model to x to forecast the index state vector. n is the time index.
* ``tgl(ℳn,x,Δx)`` apply the tangent linear model for a perturbation ``Δx`` around ``x`` (foreward differentiation).
* ``adj(ℳn,x,Δx)`` apply the adjoin model (reverse differentiation).

For any ``x``, ``Δx₁``, ``Δx₂`` and ``n``, tangent linear and adjoint model are related by:

```math
Δx₂ ⋅ tgl(ℳ,n,x,Δx₁) = adj(ℳ,n,x,Δx₂) ⋅ Δx₁
```

where ⋅ is the inner product.

The same API should also be implemented for the observation ``𝓗`` where ``𝓗(n,x)`` represents the observed part of the state vector.
Note that for the Kalman Filter method the adjoint (of the model or the observation operator) is not needed.


```@docs
FreeRun
KalmanFilter
fourDVar
```

## Ensemble methods

```@docs
ETKF
EnKF
EnSRF
EAKF
SEIK
ESTKF
serialEnSRF
local_ETKF
local_EnKF
local_EnSRF
local_EAKF
local_SEIK
local_ESTKF
```


## Models

```@docs
AbstractModel
ModelMatrix
ModelFun
LinShallowWater1DModel
Lorenz63Model
```

## Utility functions

```@docs
compact_locfun
```
