"""
This file contains all of the functions to analyze the results of simulations
not including already in the sim_fucntions.jl file.
"""

include("sim_functions.jl")
using CSV

function measure_module_interaction(S_0::Array{Array{Float64,1},1},
                                    grn::Array{Float64,2},
                                    trait_1_genes::Array{Int64,1},
                                    trait_2_genes::Array{Int64,1})
    """
    Purpose: This function measures whether there are any functional
             interactions betwen two distinct modules by deleting all
             interactions between mdoules.

    Inputs:
            S_0: array of starting gene states
            grn: gene regulatory network of interest
            trait_1_genes: list of genes that are ancestrally involved in trait 1
            trait_2_genes: list of genes that are ancestrally involved in trait 2
    Outputs:
            1) Fitness when all between-module interactions are removed
    """

    mutant_grn = deepcopy(grn)

    for row in trait_1_genes
        for col in trait_2_genes
            mutant_grn[row,col] = 0.0
        end
    end
    for row in trait_2_genes
        for col in trait_1_genes
            mutant_grn[row,col] = 0.0
        end
    end

    fitness = 1.0
    for i=1:length(S_0)
        ancestor_stability, ancestor_S_final = assess_stability(grn,S_0[i])
        mutant_stability, mutant_S_final = assess_stability(mutant_grn,S_0[i])
        if mutant_stability == 0 || mutant_S_final != ancestor_S_final
            fitness = 0.0
            break
        end
    end

    return fitness
end

function measure_redundancy(filename::String)

    """
    Purpose: This function measures how many interactions/genes, and what
    percentage of interactions/genes, have no effect on fitness when knocked-out
    across all starting gene states.

    Inputs:
            filename: fitness landscape file
    Outputs:
            1) Number of interactions/genes with no effect on fitness when KO'd
            2) Percent of interactions/genes with no effect on fitness when K0'd
    """

    df = CSV.read(filename)
    df_mutants = df[df.Column .!= 0,:] #eliminate ancestor
    df_nonzero = df_mutants[df_mutants.Value .!= 0.0,:] #eliminate non-existent interactions/genes

    fitness_values = df_nonzero[:,:Fitness]

    total = length(fitness_values)
    neutral = count(x->x==1.0,fitness_values)

    return neutral,neutral/total
end

function measure_functional_pleiotropy(filename::String,
                                       trait_1_genes::Array{Int64,1},
                                       trait_2_genes::Array{Int64,1})

    """
    Purpose: This function measures how many interactions/genes, and what
    percentage of interactions/genes, are functionally pleiotropic, defined here
    as altering the final gene state of at least one gene in both traits.
    Percentage is taken as a percentage of all non-neutral interactions/genes.

    Inputs:
            filename: fitness landscape file
            trait_1_genes: list of genes that are ancestrally involved in trait 1
            trait_2_genes: list of genes that are ancestrally involved in trait 2
    Outputs:
            1) Number of interactions/genes with a pleiotropic effect
            2) Percent of interactions/genes with a pleiotropic effect
    """

    function get_finalgenestates(fgs::String)

        """
        Purpose: Gets FinalGeneStates into an array format

        Inputs:
                fgs: String represent final gene states (one per input)

        Outputs:
                phenotype: An array of arrays of strings corresponding to the
                final gene states
        """

        phenotype = Array{Array{String,1},1}()

        if occursin(";",fgs)
            temp_phenotype = split(fgs,";")
            temp_phenotype = [String(i) for i in temp_phenotype]
            for i=1:length(temp_phenotype)
                push!(phenotype,split(temp_phenotype[i],":"))
            end
        else
            push!(phenotype,split(fgs,":"))
        end
        return phenotype
    end

    df = CSV.read(filename)

    df_ancestor = df[df.Column .== 0,:]
    anc_finalstates = df_ancestor[:,:FinalGeneStates]
    ancestor_phenotype = get_finalgenestates(anc_finalstates[1])
    #println("1")
    df_mutants = df[df.Column .!= 0,:] #eliminate ancestor
    df_deleterious = df_mutants[df_mutants.Fitness .== 0.0,:] #eliminate all nonfunctional genes/interactions
    del_finalstates = df_deleterious[:,:FinalGeneStates]
    #println("2")
    total = length(del_finalstates)
    pleio = 0
    #println("3")
    for i=1:total
        trait_1_altered = false
        trait_2_altered = false
        mutant_phenotype = get_finalgenestates(del_finalstates[i])
        #println("4")
        for j=1:length(mutant_phenotype)
            for k=1:length(trait_1_genes)
                if mutant_phenotype[j][trait_1_genes[k]] != ancestor_phenotype[j][trait_1_genes[k]]
                    trait_1_altered = true
                end
            end
            #println("5")
            for k=1:length(trait_2_genes)
                if mutant_phenotype[j][trait_2_genes[k]] != ancestor_phenotype[j][trait_2_genes[k]]
                    trait_2_altered = true
                end
            end
            #println("6")
        end
        if trait_1_altered == true && trait_2_altered == true
            pleio += 1
        end
    end
    #println("7")

    return pleio, pleio/total
end

function measure_network_complexity(S_0::Array{Array{Float64,1},1},
                                    grn::Array{Float64,2},
                                    num_trials::Int64)
    """
    Purpose: This function estimates the minimal number of interactions
    required to maintain the same final gene state as a proxy of network
    complexity. It does this by removing neutral interactions one-by-one until
    no neutral interactions remain. Because these interactions are removed
    randomly, it does it for a specified number of trials.

    Inputs:
            S_0: initial gene states
            grn: interaction matrix
            num_trials: the number of trials to perform this estimation

    Outputs:
            min_complexity: the reduced network with the minimal interactions
            across all trials
            percent_complexity: the min_complexity divided by the total number
            of interactions
    """

    genes = length(S_0[1])
    S_final = Array{Array{Float64,1},1}()
    for i=1:length(S_0)
        ancestor_stability, ancestor_S_final = assess_stability(grn,S_0[i])
        push!(S_final,ancestor_S_final)
    end

    nonzero_interactions = Array{Tuple{Int64,Int64},1}()
    for row=1:genes
        for col=1:genes
            if grn[row,col] != 0.0
                push!(nonzero_interactions,(row,col))
            end
        end
    end

    complexity_values = Array{Int64,1}()
    for i=1:num_trials
        potential_interactions = deepcopy(nonzero_interactions)
        if length(potential_interactions) == 0
            println("No potential interactions to start")
            println(S_final)
            push!(complexity_values,0)
            continue
        end
        reduced_grn = deepcopy(grn)
        no_neutrals = false

        while no_neutrals == false

            neutral_found = false
            shuffle!(potential_interactions)

            for j=1:length(potential_interactions)
                mutant_grn = deepcopy(reduced_grn)
                row = potential_interactions[j][1]
                col = potential_interactions[j][2]
                @assert mutant_grn[row,col] != 0.0
                mutant_grn[row,col] = 0.0

                is_neutral = true
                for k=1:length(S_0)
                    mutant_stability, mutant_S_final = assess_stability(mutant_grn,S_0[k])
                    if mutant_S_final != S_final[k]
                        is_neutral = false
                        break
                    end
                end

                if is_neutral == true
                    neutral_found = true
                    reduced_grn = deepcopy(mutant_grn)
                    filter!(x->x!=(row,col),potential_interactions)
                    break
                end
            end

            if neutral_found == false || length(potential_interactions) == 0
                no_neutrals = true
            end

        end

        push!(complexity_values, length(potential_interactions))
        if length(potential_interactions) == 0
            println("No potential interactions left")
            println(S_final)
        end
    end
    #println(minimum(complexity_values)," ",minimum(complexity_values)/length(nonzero_interactions))
    println(mean(complexity_values)," ",std(complexity_values))
    return minimum(complexity_values),minimum(complexity_values)/length(nonzero_interactions)
end

function measure_trait_connectivity(grn::Array{Float64,2},
                                    trait_1_genes::Array{Int64,1},
                                    trait_2_genes::Array{Int64,1})
    """
    Purpose: This function measures how many interactions occur between genes
    ancestrally in the same trait/module and how many interactions occur between
    genes in ancestrally-different modules.

    Inputs:
            grn: gene regulatory network of interest
            trait_1_genes: list of genes that are ancestrally involved in trait 1
            trait_2_genes: list of genes that are ancestrally involved in trait 2

    Outputs:
            1) Number of positive interactions between genes in ancestral module
            2) Number of negative interactions between genes in ancestral module
            3) Number of positive interactions between genes in alternative module
            4) Number of negative interactions between genes in alternative module
    """

    function get_int_sign(value::Float64)
        if value < 0.0
            return "neg"
        elseif value > 0.0
            return "pos"
        else
            return "zero"
        end
    end

    genes = length(grn[1,:])
    anc_ints = Dict("pos"=>0,"neg"=>0)
    alt_ints = Dict("pos"=>0,"neg"=>0)

    for row=1:genes
        for col=1:genes
            if grn[row,col] != 0.0
                if row in trait_1_genes && col in trait_1_genes
                    anc_ints[get_int_sign(grn[row,col])] += 1
                elseif row in trait_2_genes && col in trait_2_genes
                    anc_ints[get_int_sign(grn[row,col])] += 1
                else
                    alt_ints[get_int_sign(grn[row,col])] += 1
                end
            end
        end
    end
    return anc_ints["pos"],anc_ints["neg"],alt_ints["pos"],alt_ints["neg"]
end

function load(filename::String)

    """
    Purpose: This function loads the files saved during a simulation run. Based,
    on the file name, it will load either an input state (S_0), a gene
    regulatory network (GRN), or a final state (S_final).

    Inputs:
            1) Name of file to load (must include directory names)

    Outputs:
            1) Either an S_0, a GRN, a S_final, or a list of fixed mutatins,
            dependent on the loaded file (inferred by file name)
    """

    if isfile(filename)
        f = open(filename,"r")
        input = read(f,String)
        close(f)
        if occursin("GRN",filename)
            input = split(input,"\n")
            GRN = zeros(Float64,length(input),length(split(input[1],",")))
            for i=1:length(input)
                GRN[i,:] = [parse(Float64,j) for j in split(input[i],",")]
            end
            return GRN
        elseif occursin("S_0",filename) || occursin("S_final",filename)
            input = split(input,"\n")
            S_values = Array{Array{Float64,1},1}()
            for i=1:length(input)
                temp = split(input[i],",")
                push!(S_values,[parse(Float64,j) for j in temp])
            end
            return S_values
        elseif occursin("Mutations",filename)
            input = split(input,"\n")
            mutation_list = Array{Mutation,1}()
            for i=1:length(input)
                temp = split(input[i],",")
                #println(i," ",temp)
                if temp == [""]
                    if i != 10001
                        println(filename," ",i," ",temp)
                    end
                    continue
                end
                push!(mutation_list,Mutation(parse(Int64,temp[1]),
                                             parse(Float64,temp[2]),
                                             parse(Float64,temp[3])))
            end
            return mutation_list
        else
            throw(error("Unknown file type"))
        end
    else
        println(filename)
        throw(error("File does not exist"))
    end
end

function get_ancestor_trials(filename::String)
    """
    Purpose: This function gets the number of trials needed to generate a stable
             ancestral GRN (specifically, the second GRN).

    Inputs:
            1) Name of file containing number(must include directories)

    Outputs:
            1) Number of ancestor trials
    """
    if isfile(filename)
        f = open(filename,"r")
        input = read(f,String)
        close(f)
        input = split(input,"\n")
        @assert length(input) == 2
        return parse(Int64,input[2])
    end
    return -1
end

function analyze_fitness_landscape(landscape_type::String,
                                   S_0::Array{Array{Float64,1},1},
                                   grn::Array{Float64,2},
                                   S_final)

   """
   Purpose: This function generates the fitness landscape of either single
            interaction or single gene mutations for a given genotype. If
            instruction is 0.0 or gene has no interactions, set fitness to 1
            and final gene state to ancestral gene state

   Inputs:
           1) landscape_type: Either genes or instructions
           2) S_0: initial gene state vector
           3) grn: matrix of gene interactions
           4) S_final: final gene state vector

   Outputs:
           1) List where each entry is a dictionary for a different
              gene/instruction and the dictionary entries are the
              instruction's/genes's value, it's vector of
              final gene states, and that mutation's fitness
   """

   @assert landscape_type == "genes" || landscape_type == "interactions"


   genes = length(grn[1,:])

   if landscape_type == "interactions"

       landscape = Dict{Tuple{Int64,Int64},Dict{String,Any}}()

       for row=1:genes
           for col=1:genes

               fitness = 1.0
               final_states = Array{Array{Float64,1},1}()

               if grn[row,col] == 0.0
                   final_states = deepcopy(S_final)
               else
                   mutant_grn = deepcopy(grn)
                   mutant_grn[row,col] = 0.0

                   for i=1:length(S_0)
                       mutant_stability, mutant_S_final = assess_stability(mutant_grn,S_0[i])

                       push!(final_states,mutant_S_final)
                       if mutant_stability != true || mutant_S_final != S_final[i]
                           fitness = 0.0
                       end
                   end
               end
               landscape[(row,col)] = Dict("Value" => grn[row,col],
                                           "Final_States" => final_states,
                                           "Fitness" => fitness)
           end
       end
   elseif landscape_type == "genes"

       landscape = Dict{Tuple{Int64,Int64},Dict{String,Any}}()

       for g=1:genes

           fitness = 1.0
           final_states = Array{Array{Float64,1},1}()
           value = 0.0

           #KO gene
           mutant_grn = deepcopy(grn)
           functional = false
           for h=1:genes
               if grn[g,h] != 0.0 || grn[h,g] != 0.0
                   functional = true
               end
               #remove interactions that regulate target gene
               mutant_grn[g,h] = 0.0
               #remove interactions that target gene regulates
               mutant_grn[h,g] = 0.0
           end

           if functional == false
               final_states = deepcopy(S_final)
           else
               value = 1.0
               for i=1:length(S_0)
                   mutant_stability, mutant_S_final = assess_stability(mutant_grn,S_0[i])

                   push!(final_states,mutant_S_final)
                   if mutant_stability != true || mutant_S_final != S_final[i]
                       fitness = 0.0
                   end
               end
           end
           landscape[(g,g)] = Dict("Value" => value,
                                   "Final_States" => final_states,
                                   "Fitness" => fitness)
       end
   end

   landscape[(0,0)] = Dict("Value" => 1.0,
                           "Final_States" => S_final,
                           "Fitness" => 1.0)

   return landscape
end

function measure_alt_dfe(S_0::Array{Array{Float64,1},1},grn::Array{Float64,2},
                         trait_1_genes::Array{Int64,1},
                         trait_2_genes::Array{Int64,1})
end
