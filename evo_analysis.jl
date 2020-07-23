

include("evo_analysis_functions.jl")

struct Mutation
    location::Int64
    ancestral_value::Float64
    new_value::Float64
    fitness_diff::Float64
    incoming_interactions::Int64
    outgoing_interactions::Int64
end


function main()
    parent_folder = "results"
    # treatment_list = ["Size10","Size30","Size50","Size70","Size90","Size110"]
    treatment_list = ["Size10","Size30","Size50","Size70","Size90","Size110","Size130","Size150"]
    x_values = [10,30,50,70,90,110,130,150]

    num_mutations = assemble_series(parent_folder,treatment_list,"adaptive_combined_evolved_num_mutations")

    ancestor_conv_time = assemble_series(parent_folder,treatment_list,"hybrid_ancestor_convergence_time")
    evolved_conv_time = assemble_series(parent_folder,treatment_list,"adaptive_combined_evolved_convergence_time")

    avg_steps = assemble_series(parent_folder,treatment_list,"adaptive_combined_evolved_steps")

    avg_final_connectivity = assemble_series(parent_folder,treatment_list,"adaptive_combined_evolved_avg_connectivity")

    mutation_list = get_mutation_lists(parent_folder,treatment_list,"adaptive_combined_evolved_Fixed_Mutations")

    quadrants,incoming,outgoing,fitness_difference = series_from_mutators(mutation_list,x_values)

    produce_plots(num_mutations,x_values,false,"Number of mutations to optimal")
    produce_plots([ancestor_conv_time,evolved_conv_time],x_values,false,"Convergence time",["ancestor","evolved"])
    produce_plots(avg_steps,x_values,true,"Steps to optimal state")

    produce_plots([avg_final_connectivity,incoming,outgoing],x_values,false,
                    "Average connectivity",["All genes","Incoming mutators","Outgoing mutators"])

    produce_plots(quadrants,x_values,false,"Quadrants of mutations")

    produce_plots(fitness_difference,x_values,false,"Fitness difference of mutations")

    to_csv(avg_steps, "analysis", "average_steps.csv",x_values)
    to_csv(ancestor_conv_time, "analysis", "ancestor_conv_time.csv",x_values)
    to_csv(evolved_conv_time, "analysis", "evolved_conv_time.csv",x_values)
    to_csv(num_mutations, "analysis", "num_mutations.csv",x_values)
    to_csv(avg_final_connectivity, "analysis", "avg_final_connectivity.csv",x_values)
    to_csv(quadrants, "analysis", "quadrants.csv",x_values)
    to_csv(incoming, "analysis", "incoming.csv",x_values)
    to_csv(outgoing, "analysis", "outgoing.csv",x_values)
    to_csv(fitness_difference, "analysis", "fitness_difference.csv",x_values)


    # generate_landscape_csv(parent_folder, treatment_list)

    # return [quadrants,incoming,outgoing,fitness_difference]
end

println("starting analysis")
m = main()
