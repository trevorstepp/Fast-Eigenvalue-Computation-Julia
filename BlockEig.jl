module BlockEig

using Random
using LinearAlgebra
#using BenchmarkTools

include("build.jl")
include("algorithm.jl")
include("benchmark.jl")
include("verify.jl")

export build_block_matrix, eig_KxK_diagblocks, time_block_method, time_dense_eig

end