"""
This script calculates fitness landscapes, for both ineractions and genes, for
each replicate and each experiment per replicate for a user-provided treatment.
"""

include("analyze_functions.jl")

function main()
    """
    This function only exists to eliminate global variables
    """
    conf = ConfParse("analyze.cfg")
    parse_conf!(conf)

    user_input!(conf)

    treatment = retrieve(conf, "TREATMENT")
    replicate_names = [dir for dir in readdir(treatment) if occursin("replicate",dir)]
    for rep in replicate_names

        current_dir = treatment*"/"*rep
        println(current_dir)

        experiment_list = ["hybrid_ancestor","ancestor_1","ancestor_2",
                           "reduced_ancestor_1","reduced_ancestor_2",
                           "neutral_separate_evolved",
                           "neutral_combined_evolved"]
        if occursin("Control",treatment)
            experiment_list = ["ancestor","evolved"]
        end

        for experiment in experiment_list

            S_0 = load(current_dir*"/"*experiment*"_S_0.dat")
            grn = load(current_dir*"/"*experiment*"_GRN.dat")
            S_final = load(current_dir*"/"*experiment*"_S_final.dat")

            for lt in ["interactions","genes"]

                landscape = analyze_fitness_landscape(lt,S_0,grn,S_final)

                f = open(treatment*"/"*rep*"/"*experiment*"_"*"landscape_"*lt*".csv","w")
                write(f,"Row,Column,Value,FinalGeneStates,Fitness")

                genes = length(grn[1,:])
                for row=0:genes

                    col_values = Int64[i for i=1:genes]
                    if row == 0
                        col_values = Int64[0]
                    end

                    for col in col_values
                        if lt == "genes" && row != col
                            #if genes, row == col are the only values generated
                            continue
                        end
                        if row == 0 && col != 0
                            #row = 0 denotes ancestor
                            continue
                        end
                        write(f,"\n")
                        write(f,string(row)*","*string(col)*",")
                        write(f,string(landscape[(row,col)]["Value"])*",")
                        for i=1:length(landscape[(row,col)]["Final_States"])
                            if i != 1
                                write(f,";")
                            end
                            for j=1:length(landscape[(row,col)]["Final_States"][i])
                                if j != 1
                                    write(f,":")
                                end
                                write(f,string(landscape[(row,col)]["Final_States"][i][j]))
                            end
                        end
                        write(f,","*string(landscape[(row,col)]["Fitness"]))
                    end
                end
                close(f)
            end
        end
    end
end

main()
