"""
Data Assimiation module

see https://alexander-barth.github.io/DataAssim.jl/latest/
"""
module DataAssim
using Test
using LinearAlgebra
using Printf
using Statistics
#using DIVAnd
#import DIVAnd: pack, unpack
using Optim

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
