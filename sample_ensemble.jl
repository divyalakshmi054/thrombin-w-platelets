# load the include -
include("Include.jl")

# load the model -
model_file = joinpath(_PATH_TO_MODEL,"Feedback.bst")

# build the model -
model = build(model_file)

# load the training data -
data_file = joinpath(_PATH_TO_DATA,"Training-Composition-Transformed-w-Labels.csv")
training_df = CSV.read(data_file,DataFrame)

# which visit?
visit = 1;

visit_df = filter(:Visit => x->(x==visit), training_df) 

# size of training set -
(R,C) = size(visit_df)

# main simulation -
SF = 1e9
for i ∈ 1:R

    # build new model -
    dd = deepcopy(model)
    # setup static -
    sfa = dd.static_factors_array
    sfa[1] = visit_df[i,:TFPI]       # 1 TFPI
    sfa[2] = visit_df[i,:AT]         # 2 AT
    sfa[3] = (5e-12) * SF            # 3 TF
    sfa[6] = 0.005                   # 6 TRAUMA

    # grab the multiplier from the data -
    ℳ = dd.number_of_dynamic_states
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
    # xₒ[10] = 0.2
    xₒ[19] = visit_df[i, :PLT]       # 19 PL
    dd.initial_condition_array = xₒ

    # update α -
    α = dd.α
    α[1] = 0.7
    α[9] = 0.2

    # update G -
    G = dd.G
    AT_idx = findfirst(x->x=="AT",dd.total_species_list)
    TFPI_idx = findfirst(x->x=="TFPI",dd.total_species_list)
    FIIa_idx = findfirst(x->x=="FIIa",dd.total_species_list)
    AP_idx = findfirst(x->x=="AP",dd.total_species_list)
    PL_idx = findfirst(x->x=="PL",dd.total_species_list)
    FVIIa_idx = findfirst(x->x=="FVIIa",dd.total_species_list)
    FXa_idx = findfirst(x->x=="FXa",dd.total_species_list)
    FVa_idx = findfirst(x->x=="FVa",dd.total_species_list)

    
   # r1 -
   G[TFPI_idx,1] = -0.65;

   # r2 -
   G[AP_idx,2] = 0.01;

   # r4 -
   G[AP_idx,4] = 0.25;

   # r5 -
   G[AP_idx,5] = 0.2;
   G[FVIIa_idx,5] = 0.9;

   # r6 -
   G[FXa_idx,6] = 0.95;
   G[FVa_idx,6] = 0.9;

   # r9 -
   G[AT_idx,9] = 0.045;

   # r10 -
   G[FIIa_idx,10] = 0.01;

    # run the model -
    global (T,U) = evaluate(dd,tspan=(0.0,120.0))
    data = [T U]

    path_to_sim_data = joinpath(_PATH_TO_TMP, "SIM-visit-$(visit)-TF-$(i).csv")
    CSV.write(path_to_sim_data, Tables.table(hcat(data),header=vcat("Time",dd.list_of_dynamic_species)))

    # make plots -
    _PATH_TO_FIGS = joinpath(pwd(),"figs")
    path_to_thrombin_figs = joinpath(_PATH_TO_FIGS, "Thrombin_plot_visit_$(visit)_run$(i).png")
    savefig(plot(T, U[:,9], xticks=0.0:10:180,xlabel="Time (min)", ylabel="FIIa (nM)", title="[Thrombin] vs. time, patient $(i), visit $(visit)"), path_to_thrombin_figs)

end


