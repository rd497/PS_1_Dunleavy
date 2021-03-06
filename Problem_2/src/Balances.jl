# ----------------------------------------------------------------------------------- #
# AdjBalances: Evaluates adjoint model equations given time, state, parameter_index,
# and the data_dictionary.
#
# Type: GRN-JULIA
# Version: 1.0
#
# Input arguments:
# t  - current time
# x  - state array
# parameter_index - index of the parameter that we are calculating sc's
# data_dictionary  - Data dictionary instance (holds model parameters)
#
# Return arguments:
# dxdt - derivative array at current time step
# ----------------------------------------------------------------------------------- #
function AdjBalances(t,x,parameter_index,data_dictionary)

  # look up the number of states -
  number_of_states = data_dictionary["number_of_states"]

  # partition that state -
  state_array = x[1:number_of_states]
  sensitivity_array = x[(number_of_states+1):2*number_of_states]

  # call -
  dxdt_array = Balances(t,state_array,data_dictionary)

  # Calculate the sensitivity states -
  local_data_dictionary = deepcopy(data_dictionary)
  JM = calculate_jacobian(t,state_array,local_data_dictionary)
  BM = calculate_bmatrix(t,state_array,local_data_dictionary)

  # calulate the sensitivity state -
  dsdt_array = JM*sensitivity_array+BM[:,parameter_index]

  r_array = [dxdt_array ; dsdt_array]

  # return -
  return r_array

end


# ----------------------------------------------------------------------------------- #
# Balances: Evaluates model equations given time, state and the data_dictionary.
# Type: GRN-JULIA
# Version: 1.0
#
# Input arguments:
# t  - current time
# x  - state array
# data_dictionary  - Data dictionary instance (holds model parameters)
#
# Return arguments:
# dxdt - derivative array at current time step
# ----------------------------------------------------------------------------------- #
function Balances(t,x,data_dictionary)

  # correct for negatives -
  idx_small = find(x.<0)
  x[idx_small] = 0.0

  # Get model matricies and other required data from the data_dictionary -
  stoichiometric_matrix = data_dictionary["stoichiometric_matrix"]
	dilution_matrix = data_dictionary["dilution_matrix"]
	degradation_matrix = data_dictionary["degradation_matrix"]
  mugmax = data_dictionary["maximum_specific_growth_rate"]

  # Calculate the kinetics -
  transcription_rate_array = calculate_transcription_rates(t,x,data_dictionary)
  translation_rate_array = calculate_translation_rates(t,x,data_dictionary)
  mRNA_degradation_rate_array = calculate_mRNA_degradation_rates(t,x,data_dictionary)
  protein_degradation_rate_array = calculate_protein_degradation_rates(t,x,data_dictionary)
  background_transcription_rate_array = calculate_background_transcription_rates(t,x,transcription_rate_array,data_dictionary)
  input_array = calculate_input_array(t,x,data_dictionary)

  # Call my control function -
  control_array = Control(t,x,data_dictionary)

  # Compute the modified rate -
  transcription_rate_array = transcription_rate_array.*control_array;

  # rate array -
  txtl_rate_array = [transcription_rate_array ; translation_rate_array]
  degradation_rate_array = [mRNA_degradation_rate_array ; protein_degradation_rate_array]

  # Evaluate the balance equations -
  dxdt_array = stoichiometric_matrix*txtl_rate_array+degradation_matrix*degradation_rate_array+mugmax*dilution_matrix*x+background_transcription_rate_array+input_array

  # return -
  return dxdt_array
end
