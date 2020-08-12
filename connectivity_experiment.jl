using Random
using StatsBase


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
        S_0 = rand([0.0,1.0],params["genes"])
        ancestor_grn = generate_random_grn(params)
        ancestor_grn_stability, ancestor_S_final = assess_stability(ancestor_grn,S_0)
        trial_ct += 1
    end
    save(dir,"num_ancestor_trials","Ancestor_Trials",trial_ct)
    return trial_ct
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

function trial(genes, connectivity, num, seed)
    params = Dict{String, Any}()
    params["genes"] = genes
    params["ancestor_connectivity"] = connectivity
    params["interaction_type"] = "binary"
    Random.seed!(seed)
    total = 0
    dir = "connectivity_tests/" * string(connectivity) * "/" * string(genes)
    for trial in 1:num
        print(trial)
        trial_dir = dir*"/"*string(trial)

        if isdir(trial_dir) == false
            mkpath(trial_dir)
        end

        attempts = create_ancestor(params, trial_dir)
        total += attempts
    end
    print("\n")
    println(total)
    save(dir,"00total","none",total)
    return total
end

function main()
    totals_array = Array{Any, 1}()
    for genes in 10:20:130
        total = trial(genes, 0.75, 100, 0)
        push!(totals_array, total)
    end
    return totals_array
end

m = main()
