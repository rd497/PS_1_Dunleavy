include("Include.jl")
tic()
time_start=0.0
time_stop=240.0
time_step_size=1.0

data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

local_data_dictionary = deepcopy(data_dictionary)

for parameter_index in 1:24

data_dictionary=deepcopy(local_data_dictionary)
(T,X) = adj_washout_simulation(time_start,time_stop,time_step_size,parameter_index,data_dictionary)
data_array = [T X]
  file_path = "./sensitivity2/AdjSimulation"*string(parameter_index)*".dat"
  writedlm(file_path,data_array)

end

toc()
