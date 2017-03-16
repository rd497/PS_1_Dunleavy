function adj_washout_simulation(time_start,time_stop,time_step_size,parameter_index,data_dictionary)

  # First - run the model to steady-state w/no ATRA forcing -
  XSS = estimate_steady_state(0.001,data_dictionary)

  # Next, set the IC to the steady-state value -
  initial_condition_array = XSS;
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # What is the size of the system?
  number_of_states = data_dictionary["number_of_states"]

  adj_initial_condition_vector=zeros(2*number_of_states)
  adj_initial_condition_vector[1:number_of_states] = initial_condition_array
  adj_initial_condition_vector[number_of_states+1:2*number_of_states] = zeros(number_of_states)

  # Phase 1: Run the model 1/4 of the final time w/o ATRA
  # Run the model for a section of time w/no ATRA forcing -
  time_start_phase_1 = 0.0;
  time_stop_phase_1 = 10.0;

  # Solve the model equations -
  (T1,X1) = SolveAdjBalances(time_start_phase_1,time_stop_phase_1,time_step_size,adj_initial_condition_vector,parameter_index,data_dictionary);
  adj_initial_condition_vector=X1[end,:]
# (T1,X1) = SolveBalances(time_start_phase_1,time_stop_phase_1,time_step_size,data_dictionary);

  # Phase 2: Express gene_1 and run to end-time (ic = last point phase 1) -
  # Update the IC again -
  initial_condition_array = X1[end,1:number_of_states];
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # Grab the control parameters - turn on gene_1 =
  control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
  control_parameter_dictionary["W_gene_1_RNAP"] = 1.0;
  data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary;

  # Run the model for a section of time w/no ATRA forcing -
  time_start_phase_2 = time_stop_phase_1+time_step_size
  time_stop_phase_2 = time_start_phase_2 + 60

  # Solve the model equations -
  (T2,X2) = SolveAdjBalances(time_start_phase_2,time_stop_phase_2,time_step_size,adj_initial_condition_vector,parameter_index,data_dictionary);
  adj_initial_condition_vector=X2[end,:]
#  (T2,X2) = SolveBalances(time_start_phase_2,time_stop_phase_2,time_step_size,data_dictionary);

  # Washout -
  initial_condition_array = X2[end,:];
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # Grab the control parameters - turn off gene_1 =
  control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
  control_parameter_dictionary["W_gene_1_RNAP"] = 0.0;
  data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary;

  time_start_phase_3 = time_stop_phase_2+time_step_size
  time_stop_phase_3 = time_stop
  (T3,X3) = SolveAdjBalances(time_start_phase_3,time_stop_phase_3,time_step_size,adj_initial_condition_vector,parameter_index,data_dictionary);
  adj_initial_condition_vector=X3[end,1:2*number_of_states]
#  (T3,X3) = SolveBalances(time_start_phase_3,time_stop_phase_3,time_step_size,data_dictionary);

  # Package the two phases together -
  T = [T1; T2; T3]
  X = [X1; X2; X3]

  return (T,X)

end
