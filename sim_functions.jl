"""
This file contains the functions to run simulations using Wagner's GRN model.
"""

using Random
using Distributions
using StatsBase
using ConfParser
using DataFrames
using GLM

struct Mutation
    location::Int64
    ancestral_value::Float64
    new_value::Float64
    fitness_diff::Float64
    incoming_interactions::Int64
    outgoing_interactions::Int64
end

struct Genotype
    grn::Array{Float64,2}
    fixed_mutations::Array{Mutation,1}
end

function user_input!(param::ConfParse)
    """
    Goal: This function reads in the command line arguments, checks to make sure
    they are formatted correctly, and changes the values inputted from sim.cfg
    to user provided values.
    """
    if length(ARGS) > 0
        @assert length(ARGS)%3 == 0 "Error in ARGS"
        for i=1:convert(Int64,length(ARGS)/3)
            @assert ARGS[(i-1)*3+1] == "-set" "Error with ARG "*string(i)*" -set"
            try
                commit!(param,ARGS[(i-1)*3+2],ARGS[(i-1)*3+3])
            catch
                println("Unable to set ARG "*string(ARGS[(i-1)*3+2])*" to "*ARGS[(i-1)*3+3])
            end
        end
    end
end

function create_ancestor(params::Dict{String,Any},dir)
    """
    Goal: Create a genotype by randomly generating an array of initial gene
    states, then repeatedly generate random grns, until discovers a grn and
    corresponding final gene state that is stable and non-zero.

    Inputs:
        1. params: dictionary of parameters for the simulation

    Outputs:
        1. S_0: randomly-generated inital gene state vector
        2. ancestor_grn: randomly-generated gene regulatory network that results
           in a non-zero stable final gene state
        3. ancestor_S_final: the final gene state vector
    """

    S_0 = rand([0.0,1.0],params["genes"])

    ancestor_grn = generate_random_grn(params)
    ancestor_grn_stability, ancestor_S_final = assess_stability(ancestor_grn,S_0)
    trial_ct = 1
    #sum(ancestor_S_final) ensures ancestor_S_final is not the 0 vector
    #This only works when gene values can be either 0.0 or 1.0
    while ancestor_grn_stability == false || sum(ancestor_S_final) == 0.0
        ancestor_grn = generate_random_grn(params)
        ancestor_grn_stability, ancestor_S_final = assess_stability(ancestor_grn,S_0)
        trial_ct += 1
    end
    save(dir,"num_ancestor_trials","Ancestor_Trials",trial_ct)
    return S_0,ancestor_grn,ancestor_S_final
end

function create_sep_ancestors(params::Dict{String,Any},dir)
    """
    Goal: Create a genotype by randomly generating an array of initial gene
    states, then repeatedly generate random grns, until discovers a grn and
    corresponding final gene state that is stable and non-zero.

    Inputs:
        1. params: dictionary of parameters for the simulation

    Outputs:
        1. S_0: randomly-generated inital gene state vector
        2. ancestor_grn: randomly-generated gene regulatory network that results
           in a non-zero stable final gene state
        3. ancestor_S_final: the final gene state vector
    """

    S_0_1 = rand([0.0,1.0],convert(Int64,params["genes"]/2))
    S_0_1 = vcat(S_0_1,[0.0 for i=1:params["genes"]/2])

    S_0_2 = [0.0 for i=1:params["genes"]/2]
    S_0_2 = vcat(S_0_2,rand([0.0,1.0],convert(Int64,params["genes"]/2)))

    ancestor_grn = generate_random_grn(params)
    ancestor_grn_stability_1, ancestor_S_final_1 = assess_stability(ancestor_grn,S_0_1)
    ancestor_grn_stability_2, ancestor_S_final_2 = assess_stability(ancestor_grn,S_0_2)
    trial_ct = 1
    #sum(ancestor_S_final) ensures ancestor_S_final is not the 0 vector
    #This only works when gene values can be either 0.0 or 1.0
    while ancestor_grn_stability_1 == false || sum(ancestor_S_final_1) == 0.0 || ancestor_grn_stability_2 == false || sum(ancestor_S_final_2) == 0.0
        ancestor_grn = generate_random_grn(params)
        ancestor_grn_stability_1, ancestor_S_final_1 = assess_stability(ancestor_grn,S_0_1)
        ancestor_grn_stability_2, ancestor_S_final_2 = assess_stability(ancestor_grn,S_0_2)
        trial_ct += 1
    end
    save(dir,"num_ancestor_trials","Ancestor_Trials",trial_ct)
    return [S_0_1,S_0_2],ancestor_grn,[ancestor_S_final_1,ancestor_S_final_2]
end

function save(dir::String,filename::String,S_0::Array{Float64,1},
              grn::Array{Float64,2},S_final::Array{Float64,1})
      """
      Goal: This function saves a provided initial gene state, gene regulatory
      network, and final gene state in a provided directory and with a provided
      filename.

      Inputs:
        1. dir: directory in which to save files
        2. filename: string to start file names (usually the experiment name)
        3. S_0: array of initial gene states
        4. grn: gene regulatory network matrix
        5. S_final: array of final gene states
      """

      @assert isdir(dir)
      N = length(S_0)

      f = open(dir*"/"*filename*"_S_0.dat","w")
      for j=1:N
          if j!=1
              write(f,",")
          end
          write(f,string(S_0[j]))
      end
      close(f)

      f = open(dir*"/"*filename*"_GRN.dat","w")
      for j=1:N
          if j!= 1
              write(f,"\n")
          end
          for k=1:N
              if k!=1
                  write(f,",")
              end
              write(f,string(grn[j,k]))
          end
      end
      close(f)

      f = open(dir*"/"*filename*"_S_final.dat","w")
      for j=1:N
          if j!=1
              write(f,",")
          end
          write(f,string(S_final[j]))
      end
      close(f)
end

function save(dir::String,filename::String,S_0::Array{Float64,1},
              grn::Array{Float64,2},S_final::Array{Float64,1},
              fixed_mutations::Array{Mutation,1})
      """
      Goal: This function extends the previous save function by including a list
      of fixed mutations. It calls the previous save function and then directly
      saves only the list of fixed mutations.

      Inputs:
        1. dir: directory in which to save files
        2. filename: string to start file names (usually the experiment name)
        3. S_0: array of initial gene states
        4. grn: gene regulatory network matrix
        5. S_final: array of final gene states
        6. fixed_mutations: an array of fixed mutations
      """

      save(dir,filename,S_0,grn,S_final)

      f = open(dir*"/"*filename*"_fixed_mutations.dat","w")
      for j=1:length(fixed_mutations)
          write(f,string(fixed_mutations[j].location)*",")
          write(f,string(fixed_mutations[j].ancestral_value)*",")
          write(f,string(fixed_mutations[j].new_value)*",")
          write(f,string(fixed_mutations[j].fitness_diff)*",")
          write(f,string(fixed_mutations[j].incoming_interactions)*",")
          write(f,string(fixed_mutations[j].outgoing_interactions)*"\n")
      end
      close(f)
end

function save(dir::String,filename::String,S_0::Array{Array{Float64,1}},
              grn::Array{Float64,2},S_final::Array{Array{Float64,1}},
              fixed_mutations::Array{Mutation,1})

      """
      Goal: This function combines the functionality of the previous two save
      functions but allows the array of initial gene states and the array of
      final gene states to contain multiple arrays (to account of multiple
      environments as in the separate experimental condition.)

      Inputs:
        1. dir: directory in which to save files
        2. filename: string to start file names (usually the experiment name)
        3. S_0: array of array of initial gene states
        4. grn: gene regulatory network matrix
        5. S_final: array of array of final gene states
        6. fixed_mutations: an array of fixed mutations
      """

      @assert isdir(dir)
      N = length(S_0[1])

      f = open(dir*"/"*filename*"_S_0.dat","w")
      for j=1:length(S_0)
          if j!= 1
              write(f,"\n")
          end
          for k=1:N
              if k!=1
                  write(f,",")
              end
              write(f,string(S_0[j][k]))
          end
      end
      close(f)

      f = open(dir*"/"*filename*"_GRN.dat","w")
      for j=1:N
          if j!= 1
              write(f,"\n")
          end
          for k=1:N
              if k!=1
                  write(f,",")
              end
              write(f,string(grn[j,k]))
          end
      end
      close(f)

      f = open(dir*"/"*filename*"_S_final.dat","w")
      for j=1:length(S_final)
          if j!= 1
              write(f,"\n")
          end
          for k=1:N
              if k!=1
                  write(f,",")
              end
              write(f,string(S_final[j][k]))
          end
      end
      close(f)

      f = open(dir*"/"*filename*"_Fixed_Mutations.dat","w")
      for j=1:length(fixed_mutations)
          if j!= 1
              write(f,"\n")
          end
          write(f,string(fixed_mutations[j].location)*",")
          write(f,string(fixed_mutations[j].ancestral_value)*",")
          write(f,string(fixed_mutations[j].new_value)*",")
          write(f,string(fixed_mutations[j].fitness_diff)*",")
          write(f,string(fixed_mutations[j].incoming_interactions)*",")
          write(f,string(fixed_mutations[j].outgoing_interactions))
      end
      close(f)
end

function save(dir::String,filename::String,S_0::Array{Array{Float64,1}},
              grn::Array{Float64,2},S_final::Array{Array{Float64,1}})

      """
      Goal: This function performs the save function as above but without
      requiring a list of saved mutations.

      Inputs:
        1. dir: directory in which to save files
        2. filename: string to start file names (usually the experiment name)
        3. S_0: array of array of initial gene states
        4. grn: gene regulatory network matrix
        5. S_final: array of array of final gene states
      """

      @assert isdir(dir)
      N = length(S_0[1])

      f = open(dir*"/"*filename*"_S_0.dat","w")
      for j=1:length(S_0)
          if j!= 1
              write(f,"\n")
          end
          for k=1:N
              if k!=1
                  write(f,",")
              end
              write(f,string(S_0[j][k]))
          end
      end
      close(f)

      f = open(dir*"/"*filename*"_GRN.dat","w")
      for j=1:N
          if j!= 1
              write(f,"\n")
          end
          for k=1:N
              if k!=1
                  write(f,",")
              end
              write(f,string(grn[j,k]))
          end
      end
      close(f)

      f = open(dir*"/"*filename*"_S_final.dat","w")
      for j=1:length(S_final)
          if j!= 1
              write(f,"\n")
          end
          for k=1:N
              if k!=1
                  write(f,",")
              end
              write(f,string(S_final[j][k]))
          end
      end
      close(f)
end

function save(dir::String,filename::String,stat_name::String,stat_value)
    """
    Goal: Another save function that saves a given stat and its value to a file.
    This was developed to save how many random GRNs were generated for ancestor
    creation.

    Inputs:
        1. dir: Directory in which file will be saved
        2. filename: name of file to save
        3. stat_name: Name of statistic to save
        4. stat_value: Value of stat to save
    """

    f = open(dir*"/"*filename*".dat","w")
    if stat_name != "none"
        write(f,stat_name*"\n")
    end
    write(f,string(stat_value))
    close(f)
end

function calc_fitness(sigma::Float64,stability::Bool,S_opt::Array{Float64,1},
                      S_final::Array{Float64,1})
    """
    Goal: Calculate the fitness of a genotype (representing by its final gene
    state) in comparison to some optimal final gene state. Uses fitness function
    previously used (e.g., Siegal & Bergman 2002, PNAS)

    Inputs:
        1. sigma: variable that controls strength of stabilizing selection, or
           how deleterious a mutation is away from the optimal state
        2. stability: whether the genotype results in a stable final gene state
        3. S_opt: the optimal final gene state
        4. S_final: the genotype's final gene state
    """

    #If not stable, fitness = 0
    if stability == false
        return 0.0
    end

    #If sigma has an infinite value, then fitness landscapes is such that all
    #stable GRNs have the same fitness.
    if sigma == Inf || sigma == -Inf
        return 1.0
    end
    #Otherwise, use symmetric expoential fitness function based on distance to
    #optimum
    return exp(-1*dist(S_opt,S_final)/sigma)
end

function generate_random_grn(params::Dict{String,Any})
    """
    Goal: Randomly generates a gene regulatory network matrix based on
    simulation parameters.

    Inputs:
        1. params: dictionary of simulation parameters
    Outputs:
        1. gene regulatory network matrix
    """

    if params["interaction_type"] == "binary"
        ancestor_grn = rand([-1.0,1.0],params["genes"],params["genes"])
    elseif params["interaction_type"] == "normal"
        ancestor_grn = rand(Normal(),params["genes"],params["genes"])
    elseif params["interaction_type"] == "uniform"
        ancestor_grn = rand(Uniform(-1,1),params["genes"],params["genes"])
    else
        throw(error("Unknown Interaction type"))
    end

    num_ints = params["genes"]^2
    num_unconnected_nodes = convert(Int64,num_ints*(1-params["ancestor_connectivity"]))
    unconnected_nodes = sample(1:num_ints,num_unconnected_nodes,replace=false)

    #println(unconnected_nodes)
    for i=1:length(unconnected_nodes)
        #1-> (1,1), 2-> (2,1), 11-> (1,2), 12 -> (2,2), etc...
        row = rem((unconnected_nodes[i]-1),params["genes"])+1
        col = div((unconnected_nodes[i]-1),params["genes"])+1
        ancestor_grn[row,col] = 0.0
    end
    return ancestor_grn
end

function dist(S1::Array{Float64,1},S2::Array{Float64,1})
    """
    Goal: Calculates the distance between two arrays. Used to estimate fitness
    of genotypes.

    Inputs:
        1. & 2. arrays to compares
    Outputs:
        1. the distance between the two arrays
    """

    @assert length(S1) == length(S2)
    return sum([(S1[i] - S2[i])^2 for i=1:length(S1)])/(4*length(S1))
end

function assess_stability(test_grn::Array{Float64,2},test_S_0::Array{Float64,1})
    """
    Goal: This function measures whether a given GRN results in a stable final
    gene state vector for a given initial gene state vector.

    Inputs:
        1. test_grn: gene regulatory network matrix to be tested
        2. test_S_0: vector of initial gene states

    Outputs:
        1. boolean variable on whether the grn is stable
        2. the vector of final gene states
        3. the number of steps it took to reach a stable state
    """

    function filter_function(x::Float64)
        """
        Goal: This function normalizes the value of a gene's state to either
        0.0 or 1.0.
        """

        if x > 0
            return 1.0
        else
            return 0.0
        end
    end

    S_prev = deepcopy(test_S_0)
    for t=1:100
        S_next_temp = test_grn*S_prev
        S_next = map(filter_function,S_next_temp)

        #If stable, return stability, final gene state, and number of iterations
        #it took to reach state
        if S_next == S_prev
            return true, S_next, t
        else
            S_prev = deepcopy(S_next)
        end
    end

    #If not stable, return 0 vector for stable state
    return false,zeros(length(test_S_0)),-1
end

function rand_mutate_grn(grn::Array{Float64,2},mutation_type::String,
                         int_type::String)
    """
    Goal: Given a gene regulatory network interaction matrix and a given
    mutation type, apply a random mutation of that type to the network.

    Inputs:
        1. grn: the gene regulatory network interaction matrix
        2. mutation_type: whether mutation should be a gain-, change-, or loss-
        of-function mutation
        3. int_type: the interaction type (binary or uniform)
    Outputs:
        1. mutated gene regulatory network interaction matrix
        2. which interaction was mutated
        3. ancestral value of that interaction
        4. mutant value of the interaction
    """

    new_grn = deepcopy(grn)
    network_dim = length(new_grn[1,:])
    mutant_found = false
    random_entry = -1
    row = -1
    col = -1
    while mutant_found == false
        random_entry = rand(1:network_dim^2)
        #println(random_entry)
        row = div((random_entry-1),network_dim)+1
        col = rem((random_entry-1),network_dim)+1
        #println(row,col)
        if mutation_type == "loss" && grn[row,col] != 0.0
            new_grn[row,col] = 0.0
            mutant_found = true
        elseif mutation_type == "change" && grn[row,col] != 0.0
            #this only works for binary interactions currently
            new_grn[row,col] *= -1.0
            mutant_found = true
        elseif mutation_type == "gain" && grn[row,col] == 0.0
            if int_type == "binary"
                new_grn[row,col] = rand([-1.0,1.0])
            elseif int_type == "normal"
                new_grn[row,col] = rand(Normal())
            elseif int_type == "uniform"
                new_grn[row,col] = rand(Uniform(-1,1))
            else
                throw(error("Unknown Interaction type"))
            end
            mutant_found = true
        end

    end
    @assert grn[row,col] != new_grn[row,col]
    return new_grn,random_entry,grn[row,col],new_grn[row,col]
end

function neutral_evolve(genotype::Genotype,
                        mu_values::Dict{String,Float64},
                        S_0::Array{Array{Float64,1},1},
                        S_opt::Array{Array{Float64,1},1},
                        int_type::String)

    """
    Goal: This function applies a neutral mutation to a given gene regulatory
    network.

    Inputs:
        1. genotype: the genotype (gene regualtory network) to mutate
        2. mu_values: the raw mutational values that determine the rates at
        which gain-, change-, and loss-of-function mutations occur (these are
        then adjusting for the number of interactions in the network and their
        relative rates to each other)
        3. S_0: array of vectors of initial gene states
        4. S_opt: array of vectors of gene states with greatest fitness
        5. int_type: type of interactions (binary, uniform, etc...)

    Outputs:
        1. Genotype with mutated gene regulatory network matrix and appended
        mutation list
        2. array of final gene states vector
    """

    size = length(S_0[1])
    p,n,z = measure_connectivity(genotype.grn)

    actual_mu_values = [mu_values["loss"]*(p+n),
                        mu_values["change"]*(p+n),
                        mu_values["gain"]*z]

    relative_mu_values = [actual_mu_values[i]/sum(actual_mu_values) for i=1:3]
    mutant_found = false
    while mutant_found == false
        mutant_type = sample(["loss","change","gain"],Weights(relative_mu_values))

        mutant_grn,loc,anc_val,mut_val = rand_mutate_grn(genotype.grn,mutant_type,int_type)
        mutant_found = true
        mutant_S_final_list = Array{Array{Float64,1},1}()
        for i=1:length(S_0)
            stability,mutant_S_final = assess_stability(mutant_grn,S_0[i])
            push!(mutant_S_final_list,mutant_S_final)
            if stability == false || mutant_S_final != S_opt[i]
                mutant_found = false
                break
            end
        end
        if mutant_found == true
            #println(mutant_S_final_list)

            row = div((loc-1),size)+1
            col = rem((loc-1),size)+1

            incoming_interacts = sum([abs(x) for x in mutant_grn[row,:]])
            outgoing_interacts = sum([abs(x) for x in mutant_grn[:,col]])

            new_mut = Mutation(loc,anc_val,mut_val,0,incoming_interacts,outgoing_interacts)
            new_mut_list = deepcopy((genotype.fixed_mutations))
            push!(new_mut_list,new_mut)
            return Genotype(mutant_grn,new_mut_list),mutant_S_final_list
        end
    end
    return genotype,[]
end

function hybridize_networks(S_0_1::Array{Float64,1},grn_1::Array{Float64,2},
                            S_final_1::Array{Float64,1},
                            S_0_2::Array{Float64,1},grn_2::Array{Float64,2},
                            S_final_2::Array{Float64,1})
    """
    Goal: Combine two genetic networks into one "hybrid" network, where the
    diagonal "matrices" are the entries of both gene regulatory network
    matrices and the off-diagonal "matrices" are the 0-matrix.

    Inputs:
        1&4. S_0_j: vector of initial gene states for gene network j
        2&5. grn_j: gene regulatory network interaction matrix for gene network j
        3&6. S_final_j: vector of final gene states for gene network j

    Outputs:
        1. new_S_0: vector of initial gene states for the combined gene network
        2. new_grn: gene regulatory network interaction matrix for the combined
        gene network
        3. new_S_final: vector of final gene states for the combined gene
        network
    """

    new_S_0 = vcat(S_0_1,S_0_2)
    new_S_final = vcat(S_final_1,S_final_2)

    dim = length(grn_1[1,:])
    zero_matrix12 = zeros(dim,dim)
    zero_matrix21 = zeros(dim,dim)

    grn_1_extended = vcat(grn_1,zero_matrix12)
    grn_2_extended = vcat(zero_matrix21,grn_2)
    new_grn = hcat(grn_1_extended,grn_2_extended)

    return new_S_0,new_grn,new_S_final
end

function neutral_walk(ancestor_genotype::Genotype,mu_values::Dict{String,Float64},
                      S_0::Array{Array{Float64,1},1},
                      S_opt::Array{Array{Float64,1},1},int_type::String,
                      dir::String,filename::String,walk_length::Int64)
    """
    Purpose: Performs a neutral walk for a given genotype. A neutral walk is
    the sequential fixation of a series of neutral mutations. (In the same way
    that an adaptive walk is the sequential fixation of a series of adaptive
    mutations.)

    Inputs:
        1) ancestor_genotype: Ancestral genotype
        2) mu_values: dictionary of mutation values
        3) S_0: array of gene starting states
        4) S_opt: array of optimal gene final states
        5) int_type: Type of interactions (normal or boolean)
        6) dir: location to store results
        7) filename: name of experiment
        8) walk_length: Amount of mutations to perform in neutral walk. If Inf,
                   perform walk until number of total interactions reaches an
                   equilibrium.

    Outputs:
        1) Evolved genotype
        2) Evolved S_final(s)
    """

    genotype = deepcopy(ancestor_genotype)
    evolved_S_final = []

    ct = 0
    interaction_values = Array{Int64,1}()
    while ct < walk_length
        genotype,evolved_S_final = neutral_evolve(genotype,mu_values,S_0,S_opt,
                                                  int_type)
        p,n,z = measure_connectivity(genotype.grn)
        push!(interaction_values,p+n)
        ct += 1
    end
    save(dir,filename,S_0,genotype.grn,evolved_S_final,genotype.fixed_mutations)
    return genotype,evolved_S_final
end

function adaptive_walk(ancestor_genotype::Genotype,mu_values::Dict{String,Float64},
                      S_0::Array{Array{Float64,1},1},
                      S_opt::Array{Array{Float64,1},1},
                      S_final::Array{Array{Float64,1},1},
                      int_type::String,
                      dir::String,filename::String,walk_length::Int64,
                      selection_method::String, sigma::Float64)
    """
    Purpose: Performs a neutral walk for a given genotype. A neutral walk is
    the sequential fixation of a series of neutral mutations. (In the same way
    that an adaptive walk is the sequential fixation of a series of adaptive
    mutations.)

    Inputs:
        1) ancestor_genotype: Ancestral genotype
        2) mu_values: dictionary of mutation values
        3) S_0: array of gene starting states
        4) S_opt: array of optimal gene final states
        5) int_type: Type of interactions (normal or boolean)
        6) dir: location to store results
        7) filename: name of experiment
        8) walk_length: Amount of mutations to attempt in adaptive walk
        9) selection_method: method used to determine if stable mutation is accepted
        10) sigma: constant used to calculate fitness

    Outputs:
        1) Evolved genotype
        2) Evolved S_final(s)
    """
    size = length(S_0[1])
    _, __, t = assess_stability(ancestor_genotype.grn,S_0[1])
    save(dir,"hybrid_ancestor_convergence_time","none",t)

    genotype = deepcopy(ancestor_genotype)
    evolved_S_final = []

    ct = 0
    # interaction_values = Array{Int64,1}()
    mutant_S_final_list = Array{Array{Float64,1},1}()
    push!(mutant_S_final_list,S_final[1])
    while ct < walk_length
        genotype,mutant_S_final_list, mutant_found = adaptive_evolve(genotype,mu_values,S_0,S_opt,
                                                  int_type,selection_method, sigma,
                                                  mutant_S_final_list)
        # p,n,z = measure_connectivity(genotype.grn)
        # push!(interaction_values,p+n)
        ct += 1
        if mutant_found
            if mutant_S_final_list[end] == S_opt[end]
                println("optimal state reached in ", ct)
                save(dir,filename*"_steps","none",ct)
                ct = walk_length
            end
        end
    end
    println("mutations accepted: ", length(mutant_S_final_list)-1)
    save(dir,filename,S_0,genotype.grn,mutant_S_final_list,genotype.fixed_mutations)
    save(dir,filename*"_num_mutations","none",length(mutant_S_final_list)-1)

    _, __, t = assess_stability(genotype.grn,S_0[1])
    save(dir,filename* "_convergence_time","none",t)

    avg_connectivity = sum([abs(x) for x in genotype.grn])/size
    save(dir,filename* "_avg_connectivity","none",avg_connectivity)

    if mutant_S_final_list[end] != S_opt[end]
        println("optimal state NOT REACHED")
        save(dir,filename*"_steps","none",walk_length)
        ct = walk_length
    end
    return genotype,evolved_S_final
end

function adaptive_evolve(genotype::Genotype,
                        mu_values::Dict{String,Float64},
                        S_0::Array{Array{Float64,1},1},
                        S_opt::Array{Array{Float64,1},1},
                        int_type::String,
                        selection_method::String,
                        sigma::Float64,
                        S_final_list::Array{Array{Float64,1},1})

    """
    Goal: This function creates a mutation, and applies the chosen selection
    method to determine if the mutation fixes

    Inputs:
        1. genotype: the genotype (gene regualtory network) to mutate
        2. mu_values: the raw mutational values that determine the rates at
        which gain-, change-, and loss-of-function mutations occur (these are
        then adjusted for the number of interactions in the network and their
        relative rates to each other)
        3. S_0: array of vectors of initial gene states
        4. S_opt: array of vectors of gene states with greatest fitness
        5. int_type: type of interactions (binary, uniform, etc...)

    Outputs:
        1. Genotype with mutated gene regulatory network matrix and appended
        mutation list
        2. array of final gene states vector
    """

    p,n,z = measure_connectivity(genotype.grn)

    actual_mu_values = [mu_values["loss"]*(p+n),
                        mu_values["change"]*(p+n),
                        mu_values["gain"]*z]

    relative_mu_values = [actual_mu_values[i]/sum(actual_mu_values) for i=1:3]
    mutant_type = sample(["loss","change","gain"],Weights(relative_mu_values))

    mutant_grn,loc,anc_val,mut_val = rand_mutate_grn(genotype.grn,mutant_type,int_type)

    mutant_found = true
    new_S_final_list = deepcopy(S_final_list)
    stability = true
    mutant_S_final = []
    size = length(S_0[1])
    # evaluate mutated GRN for each S_0
    for i=1:length(S_0)
        stability,mutant_S_final = assess_stability(mutant_grn,S_0[i])
        if stability == false
            mutant_found = false
            break
        end
    end
    # using the new final state, select based on selection method
    fitness_diff = 1
    if mutant_found
        old_fitness = 0
        fitness = 0
        for i in 1:length(S_0)
            old_fitness += calc_fitness(sigma, true, S_opt[i], S_final_list[end])
            fitness += calc_fitness(sigma, true, S_opt[i], mutant_S_final)
            fitness_diff = fitness/old_fitness
        end
        # separate by selection method
        if selection_method == "absolute"
            mutant_found = absolute_selection(fitness, old_fitness)
        end
    end
    # if mutant accepted, update s_final list, and update genotype, then return
    if mutant_found
        #println(mutant_S_final_list)
        push!(new_S_final_list, mutant_S_final)

        row = div((loc-1),size)+1
        col = rem((loc-1),size)+1

        incoming_interacts = sum([abs(x) for x in mutant_grn[row,:]])
        outgoing_interacts = sum([abs(x) for x in mutant_grn[:,col]])

        new_mut = Mutation(loc,anc_val,mut_val,fitness_diff,incoming_interacts,outgoing_interacts)
        new_mut_list = deepcopy((genotype.fixed_mutations))
        push!(new_mut_list,new_mut)
        return Genotype(mutant_grn,new_mut_list),new_S_final_list, true
    end
    return genotype,new_S_final_list, false
end

function absolute_selection(fitness, old_fitness)
    return fitness > old_fitness
end

function run_sim(params::Dict{String,Any})
    """
    Purpose: This function runs 1 replicate of this simulation for the following
    two scenarios: separate neutral evolution and combined neutral evolution.

    Inputs: A dictionary containing all experiments parameters from the sim.cfg
    file

    Outputs: Saves all in the folder 'parent_folder/treatment/replicate_{i},' where i is the
    replicate number set by the random seed
    """

    #Create treatment and replicate folder if necessary
    dir = params["parent_folder"]* "/" * params["treatment"]*"/replicate_"*
          string(params["random_seed"])

    if isdir(dir) == false
        mkpath(dir)
    end


    #Some setup
    Random.seed!(params["random_seed"])
    mu_values = Dict("gain" => params["mu_gain"],
                     "change" => params["mu_change"],
                     "loss" => params["mu_loss"])

    #create two ancestor gene regulatory networks (GRNs)
    S_0_1, ancestor_grn_1, ancestor_S_final_1 = create_ancestor(params,dir)
    save(dir,"ancestor_1",S_0_1,ancestor_grn_1,ancestor_S_final_1)
    S_0_2, ancestor_grn_2, ancestor_S_final_2 = create_ancestor(params,dir)
    save(dir,"ancestor_2",S_0_2,ancestor_grn_2,ancestor_S_final_2)



    #Reduce ancestor GRNs to number of interactions established by
    #mutation-selection balance
    genotype_1 = Genotype(ancestor_grn_1,Mutation[])
    genotype_1_evolved_S_final = []
    genotype_1,genotype_1_evolved_S_final = neutral_walk(genotype_1, mu_values,
                                                         [S_0_1],
                                                         [ancestor_S_final_1],
                                                         params["interaction_type"],
                                                         dir,"reduced_ancestor_1",
                                                         params["walk_length"])
    println("neutral walk complete")
    genotype_2 = Genotype(ancestor_grn_2,Mutation[])
    genotype_2_evolved_S_final = []
    genotype_2,genotype_2_evolved_S_final = neutral_walk(genotype_2, mu_values,
                                                         [S_0_2],
                                                         [ancestor_S_final_2],
                                                         params["interaction_type"],
                                                         dir,"reduced_ancestor_2",
                                                         params["walk_length"])
    println("neutral walk complete")
    #Hybridize both reduced networks together
    hybrid_S_0,hybrid_grn,hybrid_S_final = hybridize_networks(S_0_1,
                                                              genotype_1.grn,
                                                              genotype_1_evolved_S_final[1],
                                                              S_0_2,
                                                              genotype_2.grn,
                                                              genotype_2_evolved_S_final[1])
    hybrid_stability,hybrid_final = assess_stability(hybrid_grn,hybrid_S_0)
    @assert hybrid_stability == true
    @assert hybrid_final == hybrid_S_final
    save(dir,"hybrid_ancestor",hybrid_S_0,hybrid_grn,hybrid_S_final)

    #Create alternative S_0_1 & S_0_2 that contain the ancestral S_0_i and N
    #additional 0's, where N is the number of genes. This allows for the analysis
    #of the hybrid network when only one of the signals, represented by setting
    #one of the S_0's to zero, is present at a time.
    alt_S_0_1 = vcat(S_0_1,[0.0 for i=1:length(S_0_2)])
    alt_stability_1,alt_S_final_1 = assess_stability(hybrid_grn,alt_S_0_1)
    @assert alt_stability_1 == true
    alt_S_0_2 = vcat([0.0 for i=1:length(S_0_1)],S_0_2)
    alt_stability_2,alt_S_final_2 = assess_stability(hybrid_grn,alt_S_0_2)
    @assert alt_stability_2 == true

    if params["experiment_type"] == "adaptive"
        # Adaptively evolve hybrid genotype, where the optimal phenotype is
        # the combined
        if params["reverse_signal"]
            S_0 = [alt_S_0_1, alt_S_0_2]
            S_final = [alt_S_final_1, alt_S_final_2]
        else
            S_0 = [alt_S_0_1]
            S_final = [alt_S_final_1]
        end

        genotype_adaptive = Genotype(hybrid_grn,Mutation[])
        geontype_adaptive, genotype_adaptive_S_final = adaptive_walk(genotype_adaptive,
                                                                    mu_values,
                                                                    S_0,
                                                                    [hybrid_S_final],
                                                                    S_final,
                                                                    params["interaction_type"],
                                                                    dir, "adaptive_combined_evolved",
                                                                    params["walk_length"],
                                                                    params["selection_method"],
                                                                    params["sigma"])
    end

    if params["experiment_type"] == "neutral"
        #Neutrally evolve hybrid genotype when a mutation has to be neutral both
        #when signal 1 (i.e., S_0_1) is on and signal 2 (i.e., S_0_2) is off and
        #when the opposite is true. This scenario could be imagined as each signal
        #represents a separate envionmental condition.
        genotype_sep_evo = Genotype(hybrid_grn,Mutation[])
        genotype_sep_evo_S_final = []
        genotype_sep_evo, genotype_sep_evo_S_final = neutral_walk(genotype_sep_evo,
                                                                  mu_values,
                                                                  [alt_S_0_1,alt_S_0_2],
                                                                  [alt_S_final_1,alt_S_final_2],
                                                                  params["interaction_type"],
                                                                  dir,"neutral_separate_evolved",
                                                                  params["walk_length"])
        #Neutrally evolve hybrid genotype when a mutation has to be neutral when
        #both signal 1 (i.e., S_0_1) and signal 2 (i.e., S_0_2) is on at the same
        #time. This scenario can be thought of as two GRNs being activated in the
        #same environment when they used to be activated in different environments.
        genotype_neut_evo = Genotype(hybrid_grn,Mutation[])
        genotype_neut_evo_S_final = []
        genotype_neut_evo, genotype_neut_evo_S_final = neutral_walk(genotype_neut_evo,
                                                                     mu_values,
                                                                     [hybrid_S_0],
                                                                     [hybrid_S_final],
                                                                     params["interaction_type"],
                                                                     dir,"neutral_combined_evolved",
                                                                     params["walk_length"])
    end
end

function measure_connectivity(grn::Array{Float64,2})
    """
    Goal: Return the number of positive, negative, and non-existent interactions
    in a gene regulatory interaction matrix

    Inputs:
        1. grn: gene regulatory interaction matrix
    Outputs:
        1. number of positive interactions (grn[i][j] > 0.0)
        2. number of negative interactions (grn[i][j] < 0.0)
        3. number of zero interactions (grn[i][j] == 0.0)
    """

    network_dim = length(grn[1,:])
    num_pos = 0
    num_neg = 0
    num_zero = 0
    for entry=1:network_dim^2
        row = rem((entry-1),network_dim)+1
        col = div((entry-1),network_dim)+1
        #println(entry," ",row," ",col)
        if grn[row,col] == 0
            num_zero += 1
        elseif grn[row,col] > 0
            num_pos += 1
        else
            num_neg += 1
        end
    end
    return num_pos, num_neg, num_zero
end
