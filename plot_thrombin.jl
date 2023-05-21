# include -
include("Include.jl")
size_t = 12001;

visit_1_FIIa = Array{Float64}(undef,(size_t,1))
sd_visit_1_FIIa = Array{Float64}(undef,(size_t,1))
visit_4_FIIa = Array{Float64}(undef,(size_t,1))
sd_visit_4_FIIa = Array{Float64}(undef,(size_t,1))

SIM_FIIa = Array{Float64}(undef,(size_t,12))
Z = zeros(size_t)

T_df = CSV.read(joinpath(_PATH_TO_ACTUAL_ENSEMBLE,"SIM-Pregnancy-Visit-4-1.csv"),DataFrame)
Time = T_df[1:size_t,:Time]
for i ∈ 1:12
    FIIa_df = CSV.read(joinpath(_PATH_TO_ACTUAL_ENSEMBLE,"SIM-Pregnancy-Visit-4-$(i).csv"),DataFrame)
    Temp_Thrombin = FIIa_df[1:size_t,:FIIa]
    SIM_FIIa[:,i] = Temp_Thrombin
end

visit_4_FIIa = mean(SIM_FIIa,dims = 2)
sd_visit_4_FIIa = std(SIM_FIIa,dims = 2)
U = sd_visit_4_FIIa
plot(Time,visit_4_FIIa,color=colorant"#0077BB",label="",lw=2,fillrange=(visit_4_FIIa-U,U+visit_4_FIIa),fillalpha=0.4)
xlabel!("Time (min)")
ylabel!("Activated Thrombin FIIa (nM)")
plot!([-5],[-5],xlim=(0.0,120.0),line=:scatter,markerstrokecolor=colorant"#0077BB",color=colorant"#0077BB",label = "Visit 1",foreground_color_legend = nothing)

#exp_df = CSV.read(joinpath(_PATH_TO_DATA,"thrombin-platelet-exp-data.csv"),DataFrame)
#y = scatter!(Time,exp_df[:,:FIIa],color=colorant"#0077BB",markerstrokecolor=colorant"#0077BB",label="")

#for i ∈ 1:10
#    FIIa_df = CSV.read(joinpath(_PATH_TO_ACTUAL_ENSEMBLE,"SIM-Actual-visit-3-P$(i).csv"),DataFrame)
#    Temp_Thrombin = FIIa_df[1:size_t,:FIIa]
#    SIM_FIIa[:,i] = Temp_Thrombin
#end

#visit_6_FIIa = mean(SIM_FIIa,dims = 2)
#sd_visit_6_FIIa = std(SIM_FIIa,dims = 2)
#U = sd_visit_6_FIIa
#plot!(Time,visit_6_FIIa,color=colorant"#EE7733",label="",lw=2,fillrange=(visit_6_FIIa-U,U+visit_6_FIIa),fillalpha=0.4)
#y = plot!([-5],[-5],xlim=(0.0,20.0),line=:scatter,color=colorant"#EE7733",markerstrokecolor = colorant"#EE7733",label = "Visit 6",foreground_color_legend = nothing)
savefig(y,joinpath(pwd(),"figs","Actual-N10-FIIa.pdf"))