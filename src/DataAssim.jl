"""
Data Assimiation module

see https://alexander-barth.github.io/DataAssim.jl/latest/
"""
module DataAssim
using Test
using LinearAlgebra
using Printf
using Plots
using Statistics
#using DIVAnd
#import DIVAnd: pack, unpack
using Optim
using Krylov
using LinearOperators
include("types.jl")
include("models.jl")
include("fourDVar.jl")
include("TwinExperiment.jl")
include("diagnostic.jl")

export FreeRun, fourDVar, TwinExperiment, LinShallowWater1DModel
export AbstractModel, ModelFun
export pack, unpack
export ModelMatrix
export compact_locfun

end
