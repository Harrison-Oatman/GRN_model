using Pkg


Pkg.add("StatsBase")
Pkg.add("Distributions")
Pkg.add("ConfParser")
Pkg.add("DataFrames")
Pkg.add("GLM")

println("successfully reached the end")

include("sim.jl")

println("starting main")
main()
