# load the include -
include("Include.jl")

# how many samples?
number_of_samples = 2
number_of_parameters = 19
ensemble_archive = zeros(number_of_parameters+1,1); # first row is the fitness 

# load the training data -
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Composition-Transformed-w-Labels.csv")
training_df = CSV.read(path_to_training_data,DataFrame)

# which visit?
visit = 2;

#let's filter visit 4s since we look to train using that visit
visit_df = filter(:Visit => x->(x==visit), training_df) 

# size of training set -
(R,C) = size(visit_df)

# build the model structure -
model_file = joinpath(_PATH_TO_MODEL, "Feedback.bst")
model = build(model_file)

# let us compute experimental Y here -
path_to_exp_data = joinpath(_PATH_TO_DATA,"thrombin-platelet-exp-data.csv")
exp_df = CSV.read(path_to_exp_data,DataFrame)
T = exp_df[:,1]
X = exp_df[:,2]

Y = Array{Float64,1}(undef,5)

# thrombin curve properties -

Y[1] = τlag(x->x>10.0, T, X) # lag time
Y[2] = maximum(X)            # FIIa peak
Y[3] = τpeak(T, X)           # time to peak   
Y[4] = max_FIIa_rate_exp(T, X)   # max rate
Y[5] = auc(T, X)             # area under curve

# main loop -
p_previous = nothing
for i ∈ 2:number_of_samples
    # run the learn routine -
    (p, T, U, Yₘ) = learn_optim(i, model, visit_df, Y; pₒ = nothing)

    # compute the fitness -
    fitness = norm((Yₘ .- Y).^2)
    global ensemble_archive[1] = fitness
    
    # cache the parameters -
    for k ∈ 1:number_of_parameters
        global ensemble_archive[k+1] = p[k]
    end

    # dump parameters to disk (just in case) -
    CSV.write(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "PSET-Actual-P$(i).csv"), 
        Tables.table(ensemble_archive); header=["parameters"])
    
    # dump output to disk -
    data_output = [Y Yₘ]
    data_output_header = ["actual", "simulated"]
    CSV.write(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "OUT-Actual-P$(i).csv"), 
        Tables.table(data_output); header = data_output_header)

    # dump state to disk -
    data_state = [T U]
    data_state_header = ["T", "FII", "FVII", "FV", "FX", "FVIII", "FIX", "FXI", "FXII", "FIIa", "FVIIa", "FVa", "FXa", "FVIIIa", "FIXa", "FXIa", "FXIIa", "FIIa_inactive", "AP", "PL"]
    CSV.write(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "SIM-Actual-P$(i).csv"),  
        Tables.table(hcat(data_state)); header=data_state_header)

    # update/clean up for the next patient -
    global p_previous = p
end