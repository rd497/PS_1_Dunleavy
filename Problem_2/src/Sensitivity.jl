include("Include.jl")

data_dictionary=DataDictionary(0.0,0.0,0.0)

path_to_senstivity_files= "./sensitivity"
file_pattern = "AdjSimulation"
time_skip = 1

SA=calculate_average_scaled_sensitivity_array(path_to_senstivity_files,file_pattern,data_dictionary)
writedlm("./sensitivity_results2/avg_scaled_sensitivity_repressed2.dat",SA)
SB = calculate_sensitivity_array(path_to_senstivity_files,file_pattern,time_skip,data_dictionary)
writedlm("./sensitivity_results2/SB.dat",SB)

include("PlotSensitivity.jl")

(U,S,V) = svd(SA,thin=false)

writedlm("./sensitivity_results2/U.dat",U)
writedlm("./sensitivity_results2/S.dat",S)
writedlm("./sensitivity_results2/V.dat",V)

IP = estimate_identifiable_parameters(SA,0.01)

A = SA[4,:]'
B = SA[9,:]'

D = zeros(2,24)
D[1,:] = A
D[2,:] = B

IP2 = estimate_identifiable_parameters(D,0.01)
writedlm("./sensitivity_results2/D.dat",D)

return IP, IP2
