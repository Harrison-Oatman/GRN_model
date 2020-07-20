using Plots
using StatsBase
using DataFrames
using GLM
using Statistics

struct Mutation
    location::Int64
    ancestral_value::Float64
    new_value::Float64
    fitness_diff::Float64
    incoming_interactions::Int64
    outgoing_interactions::Int64
end



function assemble_series(parent_folder::String,treatment_list,filename::String)
    """
    Produces a two-dimensional array of values with the first entry as the
    treatment, and the second entry as the replicate.
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

function produce_plots(series::Array{Array{Any,1},1},x_values::Array{Int64,1},log::Bool,title::String)
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

function produce_plots(series_list::Array{Array{Array{Any,1},1}},x_values::Array{Int64,1},
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
    display(p)
    savefig(p,"analysis/"*title*".pdf")
end

function main()
    parent_folder = "results"
    treatment_list = ["Size10","Size30","Size50"]

    x_values = [10,30,50]

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

    # produce_plots(quadrants,x_values,false,"Quadrants of mutations")

    produce_plots(fitness_difference,x_values,false,"Fitness difference of mutations")

    return [quadrants,incoming,outgoing,fitness_difference]
end

println("starting analysis")
m = main()
