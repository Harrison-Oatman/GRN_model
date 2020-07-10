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
    sim_dict["walk_length"] = parse(Int64,retrieve(conf,"WALK_LENGTH"))
    sim_dict["mu_gain"]= parse(Float64,retrieve(conf, "GAIN_MUTATION_RATE"))
    sim_dict["mu_loss"] = parse(Float64,retrieve(conf, "LOSS_MUTATION_RATE"))
    sim_dict["mu_change"] = parse(Float64,retrieve(conf, "CHANGE_MUTATION_RATE"))
    sim_dict["sigma"] = parse(Float64,retrieve(conf, "SIGMA"))
    sim_dict["treatment"] = retrieve(conf,"TREATMENT")

    sim_dict["experiment_type"] = retrieve(conf,"EXPERIMENT_TYPE")
    sim_dict["selection_method"] = retrieve(conf,"SELECTION_METHOD")

    sim_dict["reverse_signal"] = parse(Bool,retrieve(conf,"REVERSE_SIGNAL"))

    sim_dict["parent_folder"] = retrieve(conf,"PARENT_FOLDER")

    #Create treatment and replicate folder if necessary
    dir = sim_dict["parent_folder"] * "/" * sim_dict["treatment"]
    if isdir(dir) == false
        mkpath(dir)
    end

    #Record actual config parameters to know which parameters were user-inputted
    save!(conf,dir*"/actual_sim.cfg")

    for i=1:parse(Int64,retrieve(conf,"REPLICATES"))
        sim_dict["random_seed"] = parse(Int64,retrieve(conf, "RANDOM_SEED")) + i
        println("Running replicate "*string(sim_dict["random_seed"]))
        run_sim(sim_dict)
    end
end

main()
