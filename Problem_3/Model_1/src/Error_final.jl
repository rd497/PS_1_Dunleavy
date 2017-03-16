function Error_final(p)

data_dictionary = DataDictionary(0.0,0.0,0.0)

binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
data_dictionary["n_gene_2_gene_1"] = p[1]
data_dictionary["K_gene_2_gene_1"] = p[2]
data_dictionary["n_gene_3_gene_2"] = p[3]
data_dictionary["K_gene_3_gene_2"]= p[4]
data_dictionary["n_gene_3_gene_1"]= p[5]
data_dictionary["K_gene_3_gene_1"]= p[6]

control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
data_dictionary["W_gene_1_RNAP"]= p[7]
data_dictionary["W_gene_2_RNAP"]= p[8]
data_dictionary["W_gene_2_gene_1"]= p[9]
data_dictionary["W_gene_3_RNAP"]= p[10]
data_dictionary["W_gene_3_gene_2"]= p[11]
data_dictionary["W_gene_3_gene_1"]= p[12]

misc_parameter_dictionary = data_dictionary["misc_parameter_dictionary"]
data_dictionary["rnapII_concentration"] = p[13]
data_dictionary["ribosome_concentration"] = p[14]
#data_dictionary["degradation_constant_mRNA"] = p[15]
#data_dictionary["degradation_constant_protein"] = p[16]
#data_dictionary["kcat_transcription"] = p[17]
#data_dictionary["kcat_translation"] = p[18]
#data_dictionary["maximum_specific_growth_rate"] = p[19]
#data_dictionary["avg_gene_concentration"] = p[20]
#data_dictionary["saturation_constant_transcription"] = p[21]
#data_dictionary["saturation_constant_translation"] = p[22]

# Call washout
include("Washout_final.jl")
include("ErrorCalc.jl")

X = Washout_final(data_dictionary)

return X

end
