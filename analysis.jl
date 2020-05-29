"""
This file runs the analysis code for one treatment provided by the user.
This script does the following:
    1. Loads in  user-input parameters
    2. For each replicate in the input treatment, it loops through the
    different experiments for which output files exist and performs certain
    analysis calculations. All these calculations are are saved in a dictionary
    called results. Note: lanscape_analysis.jl must be run before this script.
    3. It saves all the data in the results dictionary to a csv file.
"""

include("analyze_functions.jl")

function main()
    """
    This function only exists to eliminate global variables
    """
    conf = ConfParse("analyze.cfg")
    parse_conf!(conf)

    user_input!(conf)

    analysis_dict = Dict{String,Any}()

    analysis_dict["treatment"] = retrieve(conf, "TREATMENT")
    analysis_dict["data_file_name"] = retrieve(conf,"DATA_FILE_NAME")

    results = Dict{String,Any}()

    treatment = analysis_dict["treatment"]

    replicate_names = [dir for dir in readdir(treatment) if occursin("replicate",dir)]
    for rep in replicate_names

        results[rep] = Dict{String,Any}()
        current_dir = treatment*"/"*rep
        println(current_dir)

        experiment_list = ["hybrid_ancestor","ancestor_1","ancestor_2",
                           "reduced_ancestor_1","reduced_ancestor_2",
                           "neutral_seperate_evolved",
                           "neutral_combined_evolved"]
        change_exp_basis = "hybrid_ancestor"
        if occursin("Control",treatment)
            experiment_list = ["ancestor","evolved"]
            change_exp_basis = "ancestor"
        end

        for experiment in experiment_list

            results[rep][experiment] = Dict{String,Any}()
            S_0 = load(current_dir*"/"*experiment*"_S_0.dat")
            grn = load(current_dir*"/"*experiment*"_GRN.dat")
            genes = length(S_0[1])

            p,n,z = measure_connectivity(grn)
            results[rep][experiment]["pos_interactions"] = p
            results[rep][experiment]["neg_interactions"] = n
            results[rep][experiment]["zero_interactions"] = z
            results[rep][experiment]["total_interactions"] = p+n
            results[rep][experiment]["connectivity_density"] = (p+n)/(genes^2)
            results[rep][experiment]["mean_connectivity"] = (p+n)/(genes)
            results[rep][experiment]["percent_pos_interactions"] = p/(p+n+z)
            results[rep][experiment]["percent_neg_interactions"] = n/(p+n+z)
            results[rep][experiment]["percent_zero_interactions"] = z/(p+n+z)

            results[rep][experiment]["change_pos_interactions"] = p-results[rep][change_exp_basis]["pos_interactions"]
            results[rep][experiment]["change_neg_interactions"] = n-results[rep][change_exp_basis]["neg_interactions"]
            results[rep][experiment]["change_zero_interactions"] = z-results[rep][change_exp_basis]["zero_interactions"]
            results[rep][experiment]["change_total_interactions"] = p+n-results[rep][change_exp_basis]["total_interactions"]
            results[rep][experiment]["change_connectivity_density"] = (p+n)/(genes^2)-results[rep][change_exp_basis]["connectivity_density"]
            results[rep][experiment]["change_mean_connectivity"] = (p+n)/(genes)-results[rep][change_exp_basis]["mean_connectivity"]
            results[rep][experiment]["change_percent_pos_interactions"] = (p/(p+n+z))-results[rep][change_exp_basis]["percent_pos_interactions"]
            results[rep][experiment]["change_percent_neg_interactions"] = (n/(p+n+z))-results[rep][change_exp_basis]["percent_neg_interactions"]
            results[rep][experiment]["change_percent_zero_interactions"] = (z/(p+n+z))-results[rep][change_exp_basis]["percent_zero_interactions"]

            neut,percent_neut = measure_redundancy(treatment*"/"*rep*"/"*experiment*"_landscape_interactions.csv")
            results[rep][experiment]["neut_interactions"] = neut
            results[rep][experiment]["essential_interactions"] = p+n-neut
            results[rep][experiment]["percent_essential_interactions"] = (p+n-neut)/(p+n)
            results[rep][experiment]["percent_neut_interactions"] = percent_neut
            results[rep][experiment]["change_neut_interactions"] = neut-results[rep][change_exp_basis]["neut_interactions"]
            results[rep][experiment]["change_percent_neut_interactions"] = percent_neut-results[rep][change_exp_basis]["percent_neut_interactions"]
            results[rep][experiment]["change_essential_interactions"] = p+n-neut-results[rep][change_exp_basis]["essential_interactions"]
            results[rep][experiment]["change_percent_essential_interactions"] = ((p+n-neut)/(p+n))-results[rep][change_exp_basis]["percent_essential_interactions"]

            gneut,percent_gneut = measure_redundancy(treatment*"/"*rep*"/"*experiment*"_landscape_genes.csv")
            results[rep][experiment]["neut_genes"] = gneut
            results[rep][experiment]["percent_neut_genes"] = percent_gneut
            results[rep][experiment]["change_neut_genes"] = gneut-results[rep][change_exp_basis]["neut_genes"]
            results[rep][experiment]["change_percent_neut_genes"] = percent_gneut-results[rep][change_exp_basis]["percent_neut_genes"]

            if experiment in ["ancestor_1","ancestor_2","reduced_ancestor_1","reduced_ancestor_2"]
                results[rep][experiment]["funcpleio_genes"] = -1
                results[rep][experiment]["percent_funcpleio_genes"] = -1
                results[rep][experiment]["change_funcpleio_genes"] = -1
                results[rep][experiment]["change_percent_funcpleio_genes"] = -1

                results[rep][experiment]["funcpleio_interactions"] = -1
                results[rep][experiment]["percent_funcpleio_interactions"] = -1
                results[rep][experiment]["change_funcpleio_interactions"] = -1
                results[rep][experiment]["change_percent_funcpleio_interactions"] = -1

                results[rep][experiment]["alt_pos_interactions"] = -1
                results[rep][experiment]["alt_neg_interactions"] = -1
                results[rep][experiment]["anc_neg_interactions"] = -1
                results[rep][experiment]["anc_pos_interactions"] = -1
                results[rep][experiment]["anc_total_interactions"] = -1
                results[rep][experiment]["alt_total_interactions"] = -1
                results[rep][experiment]["change_alt_pos_interactions"] = -1
                results[rep][experiment]["change_alt_neg_interactions"] = -1
                results[rep][experiment]["change_anc_neg_interactions"] = -1
                results[rep][experiment]["change_anc_pos_interactions"] = -1
                results[rep][experiment]["change_anc_total_interactions"] = -1
                results[rep][experiment]["change_alt_total_interactions"] = -1

                results[rep][experiment]["alt_neut_interactions"] = -1
                results[rep][experiment]["anc_neut_interactions"] = -1
                results[rep][experiment]["change_alt_neut_interactions"] = -1
                results[rep][experiment]["change_anc_neut_interactions"] = -1
                results[rep][experiment]["alt_essential_interactions"] = -1
                results[rep][experiment]["anc_essential_interactions"] = -1
                results[rep][experiment]["change_alt_essential_interactions"] = -1
                results[rep][experiment]["change_anc_essential_interactions"] = -1

                results[rep][experiment]["module_interaction_fitness"] = -1
            else
                num_genes = convert(Int64,length(grn[1,:])/2)
                trait_1_genes = collect(1:1:num_genes)
                trait_2_genes = collect(num_genes+1:1:num_genes*2)

                gpleio,percent_gpleio = measure_functional_pleiotropy(treatment*"/"*rep*"/"*experiment*"_landscape_genes.csv",trait_1_genes,trait_2_genes)
                results[rep][experiment]["funcpleio_genes"] = gpleio
                results[rep][experiment]["percent_funcpleio_genes"] = percent_gpleio
                results[rep][experiment]["change_funcpleio_genes"] = gpleio-results[rep][change_exp_basis]["funcpleio_genes"]
                results[rep][experiment]["change_percent_funcpleio_genes"] = percent_gpleio-results[rep][change_exp_basis]["percent_funcpleio_genes"]

                pleio,percent_pleio = measure_functional_pleiotropy(treatment*"/"*rep*"/"*experiment*"_landscape_interactions.csv",trait_1_genes,trait_2_genes)
                results[rep][experiment]["funcpleio_interactions"] = pleio
                results[rep][experiment]["percent_funcpleio_interactions"] = percent_pleio
                results[rep][experiment]["change_funcpleio_interactions"] = pleio-results[rep][change_exp_basis]["funcpleio_interactions"]
                results[rep][experiment]["change_percent_funcpleio_interactions"] = percent_pleio-results[rep][change_exp_basis]["percent_funcpleio_interactions"]

                panc,nanc,palt,nalt = measure_trait_connectivity(grn,trait_1_genes,trait_2_genes)
                @assert panc+nanc+palt+nalt == results[rep][experiment]["total_interactions"]
                results[rep][experiment]["alt_pos_interactions"] = palt
                results[rep][experiment]["alt_neg_interactions"] = nalt
                results[rep][experiment]["anc_neg_interactions"] = nanc
                results[rep][experiment]["anc_pos_interactions"] = panc
                results[rep][experiment]["anc_total_interactions"] = panc + nanc
                results[rep][experiment]["alt_total_interactions"] = palt + nalt
                results[rep][experiment]["change_alt_pos_interactions"] = palt-results[rep][change_exp_basis]["alt_pos_interactions"]
                results[rep][experiment]["change_alt_neg_interactions"] = nalt-results[rep][change_exp_basis]["alt_neg_interactions"]
                results[rep][experiment]["change_anc_neg_interactions"] = nanc-results[rep][change_exp_basis]["anc_neg_interactions"]
                results[rep][experiment]["change_anc_pos_interactions"] = panc-results[rep][change_exp_basis]["anc_pos_interactions"]
                results[rep][experiment]["change_alt_total_interactions"] = palt+nalt-results[rep][change_exp_basis]["alt_total_interactions"]
                results[rep][experiment]["change_anc_total_interactions"] = panc+nanc-results[rep][change_exp_basis]["anc_total_interactions"]

                anc_neut,p_anc_neut,alt_neut,p_alt_neut = measure_redundancy(treatment*"/"*rep*"/"*experiment*"_landscape_interactions.csv",trait_1_genes,trait_2_genes)
                @assert anc_neut+alt_neut == results[rep][experiment]["neut_interactions"]
                results[rep][experiment]["alt_neut_interactions"] = alt_neut
                results[rep][experiment]["anc_neut_interactions"] = anc_neut
                results[rep][experiment]["change_alt_neut_interactions"] = alt_neut - results[rep][change_exp_basis]["alt_neut_interactions"]
                results[rep][experiment]["change_anc_neut_interactions"] = anc_neut - results[rep][change_exp_basis]["anc_neut_interactions"]
                results[rep][experiment]["alt_essential_interactions"] = results[rep][experiment]["alt_total_interactions"] - results[rep][experiment]["alt_neut_interactions"]
                results[rep][experiment]["anc_essential_interactions"] = results[rep][experiment]["anc_total_interactions"] - results[rep][experiment]["anc_neut_interactions"]
                results[rep][experiment]["change_alt_essential_interactions"] = results[rep][experiment]["alt_essential_interactions"] - results[rep][change_exp_basis]["alt_essential_interactions"]
                results[rep][experiment]["change_anc_essential_interactions"] = results[rep][experiment]["anc_essential_interactions"] - results[rep][change_exp_basis]["anc_essential_interactions"]

                results[rep][experiment]["module_interaction_fitness"] = measure_module_interaction(S_0,grn,trait_1_genes,trait_2_genes)
            end
        end
    end


    f = open(treatment*"/"*analysis_dict["data_file_name"],"w")
    first_line = true
    for rep in keys(results)
        for experiment in keys(results[rep])
            line_string = treatment*","*rep*","*experiment*","
            first_line_string = "Treatment,Replicate,Experiment,"
            for stat in keys(results[rep][experiment])
                if first_line
                    if first_line_string == "Treatment,Replicate,Experiment,"
                        first_line_string *= stat
                    else
                        first_line_string *= ","*stat
                    end
                end
                if line_string == treatment*","*rep*","*experiment*","
                    line_string *= string(results[rep][experiment][stat])
                else
                    line_string *= ","*string(results[rep][experiment][stat])
                end
            end
            if first_line
                write(f,first_line_string*"\n")
                first_line = false
            end
            write(f,line_string*"\n")
        end
    end

    close(f)
end

main()
