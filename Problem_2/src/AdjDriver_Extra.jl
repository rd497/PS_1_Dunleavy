function AdjDriver(time_start,time_stop,time_step_size,data_dictionary)

# include -

# Setup the timescale of the simulation -
number_of_timesteps = length(time_start:time_step_size:time_stop)

# Load the data dictionary (default parameter values) -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# What is the size of the system?
number_of_states = data_dictionary["number_of_states"]

# main loop -
parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
average_scaled_sensitivity_array = zeros(number_of_states,1)

data_array=zeros(2*number_of_states+1)

for parameter_index in 1

  # grab the dictionary -
  local_data_dictionary = deepcopy(data_dictionary)

  # Solve the adj simulation -
  # You need to point this to your specific adj simulation code -

  (T,X) = adj_washout_simulation(time_start,time_stop,time_step_size,parameter_index,data_dictionary)

  # dump the raw sensitivity arrays to disk -
  # you can modify this to point to some place on disk ...
  data_array = [T X]
  file_path = "./sensitivity/AdjSimulation"*string(parameter_index)*".dat"
  writedlm(file_path,data_array)

end

end
