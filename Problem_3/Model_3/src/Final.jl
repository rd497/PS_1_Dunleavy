using NLopt
tic()

include("Include.jl")
include("Error.jl")

p0 = readdlm("./initial_p.dat")

p0 = p0[1:14]
N_param = length(p0)

opt = Opt(:LN_COBYLA, N_param)
min_objective!(opt, Error)
lower_bounds!(opt, zeros(N_param))
upper_bounds!(opt, [200, 200, 200, 200, 200, 200, 1, 1, 1, 1, 1, 1, 1e7, 9e9])
#upper_bounds!(opt, [200, 200, 200, 200, 200, 200, 1, 1, 1, 1, 1, 1, 1e7, 9e9, 10^3, 10^2, 10^3, 10^3, 10^3, 10^3, 10^6, 10^7])
#xtol_rel!(opt,1e-10)
maxtime!(opt, 300)

(minf,minx,ret) = optimize(opt, p0)

p_final = minx

include("Error_final.jl")

X_nM = Error_final(p_final) #Data is in nM, convert to uM
X = X_nM/1000

################################### Plot Final States with Optimized Parameters
const shaded_color_value = "lightgray"
const mean_color_value = "dimgray"
const experimental_color_value = "black"

const P1_color = "blue"
const P1_shaded_color="skyblue"

const P2_color = "orange"
const P2_shaded_color="navajowhite"

const P3_color = "green"
const P3_shaded_color="lightgreen"

time_array=X_nM[:,1]
P1_mean=X[:,2]; P1_std=X[:,3]
P2_mean=X[:,4]; P2_std=X[:,5]
P3_mean=X[:,6]; P3_std=X[:,7]

SF = (1.96/sqrt(10))

clf(); subplot(1,2,1)

P1_lower_bound = P1_mean - SF*P1_std
P1_upper_bound = P1_mean + SF*P1_std
plot(time_array,P1_mean,lw=2,color=P1_color)
fill_between(time_array,vec(P1_lower_bound),vec(P1_upper_bound),color=P1_shaded_color,lw=3)

P2_lower_bound = P2_mean - SF*P2_std
P2_upper_bound = P2_mean + SF*P2_std
plot(time_array,P2_mean,lw=2,color=P2_color)
fill_between(time_array,vec(P2_lower_bound),vec(P2_upper_bound),color=P2_shaded_color,lw=3)

P3_lower_bound = P3_mean - SF*P3_std
P3_upper_bound = P3_mean + SF*P3_std
plot(time_array,P3_mean,lw=2,color=P3_color)
fill_between(time_array,vec(P3_lower_bound),vec(P3_upper_bound),color=P3_shaded_color,lw=3)

xlabel("Time [hr]"); xlim(0.0, 24.0)
ylabel("Protein [uM]"); ylim(0.0, 20.0)
title("Model 3")

###################################### Plot Experimental Values
path_to_data_file = "./data/Washout-Data-PS1-Q3.dat"
X = readdlm(path_to_data_file)

time_array=X[:,1]
P1_mean=X[:,2]; P1_std=X[:,3]
P3_mean=X[:,4]; P3_std=X[:,5]
P2_mean=X[:,6]; P2_std=X[:,7]

subplot(1,2,2)

P1_lower_bound = P1_mean - SF*P1_std
P1_upper_bound = P1_mean + SF*P1_std
plot(time_array,P1_mean,lw=2,color=P1_color)
fill_between(time_array,vec(P1_lower_bound),vec(P1_upper_bound),color=P1_shaded_color,lw=3)

P2_lower_bound = P2_mean - SF*P2_std
P2_upper_bound = P2_mean + SF*P2_std
plot(time_array,P2_mean,lw=2,color=P2_color)
fill_between(time_array,vec(P2_lower_bound),vec(P2_upper_bound),color=P2_shaded_color,lw=3)

P3_lower_bound = P3_mean - SF*P3_std
P3_upper_bound = P3_mean + SF*P3_std
plot(time_array,P3_mean,lw=2,color=P3_color)
fill_between(time_array,vec(P3_lower_bound),vec(P3_upper_bound),color=P3_shaded_color,lw=3)

xlabel("Time [hr]"); xlim(0.0, 24.0)
ylim(0.0, 20.0)
title("Experimental")
####################################### Write Final Parameter
writedlm("./final_p.dat",p_final)

return p_final, Error(p_final,[1, 1]), Error(p0,[1, 1])
