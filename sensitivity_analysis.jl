# load the include -
include("Include.jl")

# build performance function -
function performance(κ, model::BSTModel, visit_df::DataFrame, i::Int64)

    # main simulation -
    SF = 1e9

    # setup static -
    sfa = model.static_factors_array
    sfa[1] = visit_df[i,:TFPI]       # 1 TFPI
    sfa[2] = visit_df[i,:AT]         # 2 AT
    sfa[3] = (5e-12) * SF            # 3 TF
    sfa[6] = 0.005                   # 6 TRAUMA
     
    # setup dynamic -
    # grab the multiplier from the data -
    ℳ = model.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = visit_df[i, :II]         # 1 FII 
    xₒ[2] = visit_df[i, :VII]        # 2 FVII 
    xₒ[3] = visit_df[i, :V]          # 3 FV
    xₒ[4] = visit_df[i, :X]          # 4 FX
    xₒ[5] = visit_df[i, :VIII]       # 5 FVIII
    xₒ[6] = visit_df[i, :IX]         # 6 FIX
    xₒ[7] = visit_df[i, :XI]         # 7 FXI
    xₒ[8] = visit_df[i, :XII]        # 8 FXII 
    xₒ[9] = (1e-14)*SF               # 9 FIIa
    xₒ[19] = 0.000274                # 19 PL
    model.initial_condition_array = xₒ
    
    #get the parameters -
    tmp_alpha = κ[1:10]
    g = κ[11:end]

    # set new parameters -
    α = model.α
    α = tmp_alpha
    model.α = α

    # set G values -
    G = model.G;

   # update G -
   G = model.G
   AT_idx = findfirst(x->x=="AT",model.total_species_list)
   TFPI_idx = findfirst(x->x=="TFPI",model.total_species_list)
   FIIa_idx = findfirst(x->x=="FIIa",model.total_species_list)
   AP_idx = findfirst(x->x=="AP",model.total_species_list)
   FVIIa_idx = findfirst(x->x=="FVIIa",model.total_species_list)
   FXa_idx = findfirst(x->x=="FXa",model.total_species_list)
 
    # adjusting parameters for r1
    G[TFPI_idx, 1] = -1*g[1]
 
    # adjusting parameters for r2
    G[AP_idx,2] = g[2]  
 
    # adjusting parameters for r4
    G[AP_idx,4] = g[3]   
    
    # adjusting parameters for r5
    G[AP_idx,5] = g[4]
    G[FVIIa_idx,5] = g[5]

    # adjusting parameters for r6
    G[FXa_idx,6] = g[6]

    # adjusting parameters for r9
    G[AT_idx,9] = g[7]

    # adjusting parameters for r10
    G[FIIa_idx,10] = g[8]

    # put it back -
    model.G = G;

    # solve -
    # run the model -
    global (T,U) = evaluate(model,tspan=(0.0,120.0))
    data = [T U]
    # test -
    return integrate(T,U[:,9])    # AUC
end

# build the model structure -
path_to_model_file = joinpath(pwd(), "model", "Feedback.bst")

# build the default model structure -
model = build(path_to_model_file)

# load the training data -
_PATH_TO_DATA = joinpath(pwd(),"data")
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Thrombin-TF-Transformed-w-Labels.csv")
training_df = CSV.read(path_to_training_data, DataFrame)

# which visit?
visit = 2;

# let's filter visit 4s since we look to train using that visit
visit_df = filter(:visitid => x->(x==visit), training_df) 


# size of training set -
(R,C) = size(visit_df)

a = ones(10)
a[1] = 0.7
a[9] = 0.2

#update G -
g = [0.65, 0.01, 0.25, 0.05, 0.9, 0.75, 0.045, 0.01] # look at sample_ensemble.jl for specific G values

# fusion -
parameters = vcat(a,g)

np = length(parameters)

L = zeros(np)
U = zeros(np)
for pᵢ ∈ 1:(np)
    L[pᵢ] = 0.01*parameters[pᵢ]
    U[pᵢ] = 10.0*parameters[pᵢ]
end
#L[end] = -3.0;
#U[end] = 0.0;

patient_index = 1;
samples = 1000;


# setup call to Morris method -
F(parameters) =  performance(parameters, model, visit_df, patient_index)
m = gsa(F, Morris(num_trajectory=samples), [[L[i],U[i]] for i in 1:np], relative_scale = false);
means = transpose(m.means)
means_star =  transpose(m.means_star)
variances = transpose(m.variances)
results_array = hcat(means, means_star, variances)

# dump sensitivity data to disk -
 CSV.write(joinpath(pwd(),"sensitivity","Sensitivity-Morris-test-$(samples).csv"), Tables.table(results_array), header = vcat("mean", "mean_star","variance"))