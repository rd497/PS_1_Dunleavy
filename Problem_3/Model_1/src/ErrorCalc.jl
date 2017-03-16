function ErrorCalc(P1::Array{Float64,1},P2::Array{Float64,1},P3::Array{Float64,1})

# Read data from file
path_to_data_file = "./data/Washout-Data-PS1-Q3.dat"
X = readdlm(path_to_data_file)

time_array=X[:,1]

P1_data=X[:,2]; P1_std=X[:,3]
P3_data=X[:,4]; P3_std=X[:,5]
P2_data=X[:,6]; P2_std=X[:,7]

totalE=(1000*P1_data-(P1)).^2+(1000*P2_data-(P2)).^2+(1000*P3_data-(P3)).^2

sumE=sum(totalE)
RMS=(sumE/length(P1))^(0.5)

return RMS
end
