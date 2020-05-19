# GRN_model

This repository contains the files to perform simulations of the evolution of gene regualtory networks. It uses a variant of a model developed by Andreas Wagner in 1994 (Wagner 1994, PNAS). The simulations are written in Julia 1.0.

To run the simulations and perform analyses of the results, run the files in the following order:

  1. Run sim.jl to perform the simulations. 
    $ julia sim.jl -set RANDOM_SEED 100 -set TREATMENT Genes10_Mu110 -set GENES 10
    
  2. Run landscape_analysis.jl to calculate the fitness landscape of evolved networks
    $ julia landscape_analysis.jl -set TREATMENT Genes10_Mu110
  
  3. Run analysis.jl to calculate statistics of the evolved networks, including data generated in the previous step.
    $ julia analysis.jl -set TREATMENT Genes10_Mu110
  
  4. Run combine_analysis.jl to combine the CSV files generated in step 3 into one file.
