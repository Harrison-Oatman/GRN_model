"""
This file runs the simulation. It:
    1) loads the configuration parameters
    2) sets them to any user-inputted values
    3) creates directory for simulation
    4) saves configuration parameters to text file
    5) runs the simulation for a set number of replicates
"""

include("sim_functions.jl")

function main()
    """
    This function only exists to eliminate global variables
    """

    conf = ConfParse("sim.cfg")
    parse_conf!(conf)

    user_input!(conf)

    #Create dictionary to store simulation variables
    sim_dict = Dict{String,Any}()
    sim_dict["random_seed"] = parse(Int64,retrieve(conf, "RANDOM_SEED"))
    sim_dict["genes"] = parse(Int64,retrieve(conf, "GENES"))
    sim_dict["ancestor_connectivity"] = parse(Float64,retrieve(conf, "ANCESTOR_CONNECTIVITY"))
    sim_dict["interaction_type"] = retrieve(conf,"INTERACTION_TYPE")
    sim_dict["neutral_walk_length"] = parse(Int64,retrieve(conf,"WALK_LENGTH"))
    sim_dict["mu_gain"]= parse(Float64,retrieve(conf, "GAIN_MUTATION_RATE"))
    sim_dict["mu_loss"] = parse(Float64,retrieve(conf, "LOSS_MUTATION_RATE"))
    sim_dict["mu_change"] = parse(Float64,retrieve(conf, "CHANGE_MUTATION_RATE"))
    sim_dict["sigma"] = parse(Float64,retrieve(conf, "SIGMA"))
    sim_dict["treatment"] = retrieve(conf,"TREATMENT")

    #Create treatment and replicate folder if necessary
    dir = sim_dict["treatment"]
    if isdir(dir) == false
        mkpath(dir)
    end

    #Record actual config parameters to know which parameters were user-inputted
    save!(conf,sim_dict["treatment"]*"/actual_sim.cfg")

    for i=1:parse(Int64,retrieve(conf,"REPLICATES"))
        sim_dict["random_seed"] = parse(Int64,retrieve(conf, "RANDOM_SEED")) + i
        println("Running replicate "*string(sim_dict["random_seed"]))

        dir = sim_dict["treatment"]*"/replicate_"*string(sim_dict["random_seed"])
        if isdir(dir) == false
            mkpath(dir)
        end

        Random.seed!(sim_dict["random_seed"])
        mu_values = Dict("gain" => sim_dict["mu_gain"],
                         "change" => sim_dict["mu_change"],
                         "loss" => sim_dict["mu_loss"])

        if occursin("Sep",sim_dict["treatment"]) == false
            S_0, ancestor_grn, ancestor_S_final = create_ancestor(sim_dict)
            save(dir,"ancestor",S_0,ancestor_grn,ancestor_S_final)

            genotype = Genotype(ancestor_grn,Mutation[])
            genotype_evolved_S_final = []
            genotype,genotype_evolved_S_final = neutral_walk(genotype, mu_values,
                                                             [S_0],
                                                             [ancestor_S_final],
                                                             sim_dict["interaction_type"],
                                                             dir,"evolved",
                                                             sim_dict["neutral_walk_length"])
        else
            S_0, ancestor_grn, ancestor_S_final = create_sep_ancestors(sim_dict)
            save(dir,"ancestor",S_0,ancestor_grn,ancestor_S_final)

            genotype = Genotype(ancestor_grn,Mutation[])
            genotype_evolved_S_final = []
            genotype,genotype_evolved_S_final = neutral_walk(genotype, mu_values,
                                                             S_0,
                                                             ancestor_S_final,
                                                             sim_dict["interaction_type"],
                                                             dir,"evolved",
                                                             sim_dict["neutral_walk_length"])
        end
    end
end

main()
