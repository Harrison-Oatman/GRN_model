"""
This file combines all the output data from each treatment's analsis.jl run
into one csv file
"""

using CSV
using DataFrames

filenames = ["Genes10_Mu110/results.csv","Genes20_Mu110/results.csv",
             "Genes30_Mu110/results.csv","Genes40_Mu110/results.csv",
             "Genes50_Mu110/results.csv","Genes60_Mu110/results.csv",
             "Genes70_Mu110/results.csv","Genes80_Mu110/results.csv",
             "Genes90_Mu110/results.csv","Genes100_Mu110/results.csv"]
             """
             "Genes10_Mu1100/results.csv","Genes20_Mu1100/results.csv",
             "Genes30_Mu1100/results.csv","Genes40_Mu1100/results.csv",
             "Genes50_Mu1100/results.csv","Genes60_Mu1100/results.csv",
             "Genes70_Mu1100/results.csv","Genes80_Mu1100/results.csv",
             "Genes90_Mu1100/results.csv","Genes100_Mu1100/results.csv",
             "Genes10_Mu11/results.csv","Genes20_Mu11/results.csv",
             "Genes30_Mu11/results.csv","Genes40_Mu11/results.csv",
             "Genes50_Mu11/results.csv","Genes60_Mu11/results.csv",
             "Genes70_Mu11/results.csv","Genes80_Mu11/results.csv",
             "Genes90_Mu11/results.csv","Genes100_Mu11/results.csv]
             """

end_filename = "all_results.csv"

all_data = DataFrame()

for i=1:length(filenames)
    f = CSV.read(filenames[i])
    if i == 1
        global all_data = f
    else
        global all_data = vcat(all_data, f)
    end
end

CSV.write(end_filename,all_data)
