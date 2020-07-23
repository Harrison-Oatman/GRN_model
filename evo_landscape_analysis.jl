include("evo_analysis_functions.jl")

function main()
    parent_folder = "results"
    treatment_list = ["Size10","Size30","Size50","Size70","Size90","Size110","Size130","Size150"]
    # treatment_list = ["Size10"]
    x_values = [10,30,50,70,90,110,130,150]
    # x_values = [10]

    generate_partial_landscape(parent_folder, treatment_list)

    prop_pos = assemble_series(parent_folder,treatment_list,"prop_pos")
    prop_neg = assemble_series(parent_folder,treatment_list,"prop_neg")
    prop_neu = assemble_series(parent_folder,treatment_list,"prop_neu")
    prop_complete = assemble_series(parent_folder,treatment_list,"prop_complete")

    produce_plots([prop_pos,prop_neg,prop_neu,prop_complete],x_values, false,
                 "Proportion of Quadrant 3 Fitness Effects",
                 ["positive (includes complete)", "negative", "neutral", "complete"],
                 false, true)
end

main()
