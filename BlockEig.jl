module BlockEig

using Random
using LinearAlgebra
#using BenchmarkTools

include("build.jl")
include("algorithm.jl")
include("verify.jl")

export build_block_matrix, eig_KxK_diagblocks, verify_results

end