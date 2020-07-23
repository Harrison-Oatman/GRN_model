using Plots
using StatsBase
using DataFrames
using GLM
using Statistics

include("sim_functions.jl")

function to_csv(series, parent_folder::String, filename::String, x_values)
    f = open(parent_folder*"/"*filename,"w")

    for (i,treatment) in enumerate(series)
        write(f, string(x_values[i]))
        for stat in treatment
            write(f, ","*string(stat))
        end
        write(f,"\n")
    end

    close(f)
end

function to_csv(array::Array{Float64,2}, path::String)
    N = length(array[1,:])
    f = open(path,"w")
    for i in 1:N
        for j in 1:N
            if j != 1
                write(f,",")
            end
            write(f,string(array[i,j]))
        end
        write(f,"\n")
    end

    close(f)
end

function assemble_series(parent_folder::String,treatment_list,filename::String)
    """
    Produces a two-dimensional array of values A[treatment,replicate]
    """
    results = Array{Array{Float64,1},1}()
    for treatment in treatment_list
        stat_list = Array{Float64,1}()
        treatment_dir = parent_folder * "/" * treatment

        replicate_list = [dir for dir in readdir(treatment_dir) if occursin("replicate",dir)]
        for replicate in replicate_list
            file = treatment_dir * "/" * replicate * "/" * filename * ".dat"
            open(file, "r") do f
                value = parse(Float64,readline(f))
                push!(stat_list,value)
            end
        end
        push!(results,stat_list)
    end
    return results
end

function get_mutation_lists(parent_folder::String,treatment_list,filename::String)
    """
    Produces a two-dimensional array of values with the first entry as the
    treatment, and the second entry as the replicate.
    """
    results = Array{Array{Array{Mutation,1},1},1}()
    for treatment in treatment_list
        stat_list = Array{Array{Mutation,1},1}()
        treatment_dir = parent_folder * "/" * treatment

        replicate_list = [dir for dir in readdir(treatment_dir) if occursin("replicate",dir)]
        for replicate in replicate_list
            mutation_list = Array{Mutation,1}()
            file = treatment_dir * "/" * replicate * "/" * filename * ".dat"
            f = open(file, "r")
            for ln in eachline(f)
                lnarr = split(ln,",")
                if lnarr != [""]
                    mutation = Mutation(parse(Int64,lnarr[1]),parse(Float64,lnarr[2]),parse(Float64,lnarr[3]),
                                        parse(Float64,lnarr[4]),parse(Int64,lnarr[5]),parse(Int64,lnarr[6]))
                    push!(mutation_list,mutation)
                end
            end
            close(f)
            push!(stat_list,mutation_list)
        end
        push!(results,stat_list)
    end
    return results
end

function produce_plots(series,x_values::Array{Int64,1},log::Bool,title::String)
    gr()
    println("here")
    mean_list = Array{Float64,1}()
    error_list = Array{Float64,1}()
    if log
        p = plot(yaxis=:log)
    else
        p = plot()
    end

    plot!(p,title=title)

    for (i, arr) in enumerate(series)
        sample_mean = mean(arr)
        sample_error = std(arr)/sqrt(length(arr))

        push!(mean_list,sample_mean)
        push!(error_list,sample_error)

        x = [x_values[i]+0.05*(j-50) for (j,_) in enumerate(arr)]

        plot!(p,x,arr,seriestype = :scatter,label="")
    end
    plot!(p,x_values,mean_list,label="Mean")
    display(p)
    savefig(p,"analysis/"*title*".pdf")
end

function produce_plots(series_list::Array{Array{Array{Float64,1},1}},x_values::Array{Int64,1},
                        log::Bool,title::String, label_list::Array{String,1})
    gr()
    println("here")

    if log
        p = plot(yaxis=:log)
    else
        p = plot()
    end

    plot!(p,title=title)
    for (series_inx,series) in enumerate(series_list)
        mean_list = Array{Float64,1}()
        error_list = Array{Float64,1}()
        for (i, arr) in enumerate(series)
            sample_mean = mean(arr)
            sample_error = std(arr)/sqrt(length(arr))

            push!(mean_list,sample_mean)
            push!(error_list,sample_error)

            x = [x_values[i]+0.05*(j-50) for (j,_) in enumerate(arr)]

            plot!(p,x,arr,seriestype = :scatter,label="")
        end
        plot!(p,x_values,mean_list,label=label_list[series_inx])
    end
    # display(p)
    savefig(p,"analysis/"*title*".pdf")
end

function series_from_mutators(mutation_list,size_list)
    quadrants = Array{Array{Int64,1},1}()
    incoming = Array{Array{Int64,1},1}()
    outgoing = Array{Array{Int64,1},1}()
    fitness_difference = Array{Array{Float64,1},1}()

    for (treatment_inx,treatment) in enumerate(mutation_list)
        size = size_list[treatment_inx]
        treatment_quadrant = Array{Int64,1}()
        treatment_incoming = Array{Int64,1}()
        treatment_outgoing = Array{Int64,1}()
        treatment_fitness_difference = Array{Float64,1}()

        for replicate in treatment
            for mutation in replicate
                push!(treatment_incoming,mutation.incoming_interactions)
                push!(treatment_outgoing,mutation.outgoing_interactions)
                push!(treatment_fitness_difference,mutation.fitness_diff)

                loc = mutation.location

                row = div((loc-1),size*2)+1
                col = rem((loc-1),size*2)+1

                if row < size
                    if col < size
                        quadrant = 1.0
                    else
                        quadrant = 2.0
                    end
                else
                    if col < size
                        quadrant = 3.0
                    else
                        quadrant = 4.0
                    end
                end

                push!(treatment_quadrant,float(quadrant))
            end
        end

        push!(quadrants,treatment_quadrant)
        push!(incoming,treatment_incoming)
        push!(outgoing,treatment_outgoing)
        push!(fitness_difference,treatment_fitness_difference)
    end
    return quadrants,incoming,outgoing,fitness_difference
end

function generate_landscape_csv(parent_folder::String,treatment_list)
    for treatment in treatment_list
        replicate_list = [dir for dir in readdir(parent_folder*"/"*treatment) if occursin("replicate",dir)]
        for replicate in replicate_list
            path = parent_folder*"/"*treatment*"/"*replicate*"/"
            print(path*"hybrid_ancestor_GRN.dat")
            hybrid_GRN = load(path*"hybrid_ancestor_GRN.dat")
            S_0 = load(path*"adaptive_combined_evolved_S_0.dat")[1]
            S_opt = load(path*"hybrid_ancestor_S_final.dat")[1]

            activation = hypothetical_fitness(hybrid_GRN, S_0, S_opt, 1.0)
            repression = hypothetical_fitness(hybrid_GRN, S_0, S_opt, -1.0)
            loss_of_interaction = hypothetical_fitness(hybrid_GRN, S_0, S_opt, 0.0)

            to_csv(activation, path*"activation.csv")
            to_csv(repression, path*"repression.csv")
            to_csv(loss_of_interaction, path*"loss_of_interaction.csv")
        end
    end
end

function hypothetical_fitness(GRN::Array{Float64,2}, S_0::Array{Float64,1},
                              S_opt::Array{Float64,1}, value::Float64)
    N = length(S_0)
    results = zeros(Float64, N, N)
    for i in 1:N
        for j in 1:N
            if GRN[i,j] == value
                results[i,j] = 100
            else
                mutated_GRN = deepcopy(GRN)
                mutated_GRN[i,j] = value
                stability, S_next, __ = assess_stability(mutated_GRN, S_0)
                fitness = calc_fitness(1.0, stability,S_opt, S_next)
                results[i,j] = fitness
            end
        end
    end
    return results
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
        else
            throw(error("Unknown file type"))
        end
    else
        println(filename)
        throw(error("File does not exist"))
    end
end
