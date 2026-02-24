using CSV, DataFrames
using Plots

# load saved data
df = CSV.read("timings.csv", DataFrame)

# create plot with block algorithm times 
plot(
    df.n,
    df.block_time,
    yscale = :log10,
    marker = :circle,
    label = "Block approach",
    xlabel = "discretization size (n)",
    ylabel = "run time (s)",
    title = "All eigenpairs runtime (K=3)",
    grid = true,
    minorgrid = true,
    legend = :topleft
)

# add dense algorithm times
plot!(
    df.n,
    df.dense_time,
    marker = :square,
    label = "eigen"
)

savefig("K=3.png")
display(current())