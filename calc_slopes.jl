using CSV
using DataFrames
using Statistics
using Plots

df = CSV.read("julia_timings.csv", DataFrame)
n_vals = df.n 
block = df.block_time
dense = df.dense_time

slope_block = cov(log.(n_vals), log.(block)) / var(log.(n_vals))
slope_dense = cov(log.(n_vals), log.(dense)) / var(log.(n_vals))

println("Block slope = $slope_block")
println("Dense slope = $slope_dense")

p = scatter(n_vals, block, xscale=:log10, yscale=:log10, label="Block", 
    marker=:circle, grid = true, minorgrid = true, legend = :topleft)
scatter!(n_vals, dense, xscale=:log10, yscale=:log10, label="Dense", marker=:diamond)

# reference lines
plot!(n_vals, n_vals .^ 1 .* block[1]/n_vals[1], linestyle=:dash, label="O(n)")
plot!(n_vals, n_vals .^ 3 .* dense[1]/n_vals[1]^3, linestyle=:dash, label="O(n^3)")

xlabel!("n")
ylabel!("Time (s)")
title!("Scaling of Block vs Dense Eigensolvers")
savefig(p, "scaling_plot.png")