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
    xₒ[19] = visit_df[i, :PLT]       # 19 PL
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
    PL_idx = findfirst(x->x=="PL",model.total_species_list)
    FVIIa_idx = findfirst(x->x=="FVIIa",model.total_species_list)
    FXa_idx = findfirst(x->x=="FXa",model.total_species_list)
    FVa_idx = findfirst(x->x=="FVa",model.total_species_list)
 
    G[TFPI_idx,1] = -1*g[1];

    # r2 -
    G[AP_idx,2] = g[2];
 
    # r4 -
    G[AP_idx,4] = g[3];
 
    # r5 -
    G[AP_idx,5] = g[4];
    G[FVIIa_idx,5] = g[5];
 
    # r6 -
    G[FXa_idx,6] = g[6];
    G[FVa_idx,6] = g[7];
 
    # r9 -
    G[AT_idx,9] = g[8];
 
    # r10 -
    G[FIIa_idx,10] = g[9];
 
    # put it back -
    model.G = G;

    # solve -
    # run the model -
    global (T,U) = evaluate(model,tspan=(0.0,180.0))
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
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Composition-Transformed-w-Labels.csv")
training_df = CSV.read(path_to_training_data, DataFrame)

# which visit?
visit = 4;

# let's filter visit 4s since we look to train using that visit
visit_df = filter(:Visit => x->(x==visit), training_df) 


# size of training set -
(R,C) = size(visit_df)

a = ones(10)
a[1] = 0.7
a[9] = 0.2

#update G -
g = [0.65, 0.01, 0.25, 0.2, 0.9, 0.95, 0.9, 0.045, 0.01] # look at sample_ensemble.jl for specific G values

# fusion -
parameters = vcat(a,g)

np = length(parameters)

L = zeros(np)
U = zeros(np)
for pᵢ ∈ 1:(np)
    L[pᵢ] = 0.5*parameters[pᵢ]
    U[pᵢ] = 2.0*parameters[pᵢ]
end
#L[end] = -3.0;
#U[end] = 0.0;

patient_index = 3;
samples = 1000;
bootreps = 100;
# initialize -
sampler = SobolSample()
    
# generate a sampler -
(A,B) = QuasiMonteCarlo.generate_design_matrices(samples,L,U,sampler)

# setup call to Sobol method -
F(parameters) =  performance(parameters, model, visit_df, patient_index)
m = gsa(F,Sobol(order = [0,1], nboot = bootreps, conf_level = 0.95), A,B)

# dump -
results_array = hcat(m.ST, m.S1, m.ST_Conf_Int, m.S1_Conf_Int)
CSV.write(joinpath(pwd(),"sobol","Sensitivity-Sobol-$(samples)-boot-$(bootreps).csv"), Tables.table(results_array), header = vcat("Total_order", "First_order", "Total_order_CI", "First_order_CI"))