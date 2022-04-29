using CSV
using DataFrames
using Distributions
using Survival: EventTime, CoxModel, fit

"""
Right-censored time-to-event data with predictors (size N_events by N_pred)
"""
struct SurvivalData
    events::AbstractVector{EventTime{<:Real}}
    predictors::AbstractMatrix{<:Real}
    SurvivalData(events, predictors) = begin
        if length(events) != size(predictors)[1]
            error("Number of rows of predictors must match the number of events")
        else
            new(events, predictors)
        end
    end
end
SurvivalData(times::AbstractVector{<:Real},
             observed::AbstractVector{Bool},
             predictors::AbstractMatrix{<:Real}) = begin
                 if length(times) != length(observed)
                     error("Number of times must match length of observed flags")
                 else
                     events = [EventTime(times[i], observed[i]) for i in 1:length(times)]
                     return SurvivalData(events, predictors)
                 end
             end
Base.length(data::SurvivalData) = length(data.events)

"""
Fit a Cox proportional hazards model to the data.
"""
function coxph(data::SurvivalData; predictor_idx = nothing)
    if predictor_idx == nothing
        x = data.predictors
    else
        x = reshape(data.predictors[:, predictor_idx],
                    (length(data), length(predictor_idx)))
    end
    return fit(CoxModel, x, data.events)
end

struct TrialDesign
    N_pt::Int
    N_bm::Int
    N_tx::Int
    duration::Real
    bm_generator
    tx_generator
end

function simulate(design::TrialDesign,
                  β_0, β_bm, β_tx, β_int, log_σ_u, log_α)::SurvivalData
    N_pt = design.N_pt
    N_bm = design.N_bm
    N_tx = design.N_tx
    x_bm = zeros(N_pt, N_bm)
    x_tx = zeros(N_pt, N_tx)
    x_int = zeros(N_pt, N_bm * N_tx) 
    for pt in 1:N_pt
        x_bm[pt, :] = design.bm_generator()
        x_tx[pt, :] = design.tx_generator()
        x_int[pt, :] = [bm * tx for bm in x_bm[pt, :] for tx in x_tx[pt, :]]
    end
    x = hcat(ones(N_pt), x_bm, x_tx, x_int)
    β = vcat(β_0, β_bm, β_tx, β_int)
    u = exp(log_σ_u) * randn(N_pt)
    α = exp(log_α)    
    θ = exp.(x * β .+ u) .^ (-1 / α)
    times = [rand(Weibull(α, θ[pt])) for pt in 1:N_pt]
    observed = times .<= design.duration
    censored = .~ observed
    if any(censored)
        times[censored] = design.duration * ones(count(censored))
    end
    return SurvivalData(times, observed, x)
end

function simulate(design::TrialDesign,
                  θ::AbstractVector{<:Real})::SurvivalData
    N_pt = design.N_pt
    N_bm = design.N_bm
    N_tx = design.N_tx
    @assert size(θ)[1] == 1 + N_bm + N_tx + N_bm * N_tx + 2
    parameter_sizes = (1, N_bm, N_tx, N_bm * N_tx, 1, 1)
    parameters = []
    total = 1
    for (i, N) in enumerate(parameter_sizes)
        if N == 1
            p = θ[total]
        else
            p = θ[total:(total + N - 1)]
        end
        push!(parameters, p)
        total += N
    end
    return simulate(design, parameters...)
end

"""
Design for Temozolomide vs Radiotherapy clinical trial presented in Malmstrom+2012

http://www.sciencedirect.com/science/article/pii/S1470204512702656
DOI: 10.1016/S1470-2045(12)70265-6
zotero: https://www.zotero.org/groups/2336486/gcta/items/itemKey/MJ99D3L4

There are four biomarker dummy variables, corresponding to the variables age and
MGMT status.  Age is coded as 0 for age ∈ (60, 70) years and as 1 for 
age > 70 years.  MGMT status is coded as 0 for unmethylated (uMGMT) and 1 for 
methyldated (mMGMT)

There are three treatment dummy variables with one-hot-encoding such that one
and only one treatment is indicated.  The treatments in order are
1. Temozolomide (TMZ)
2. Hypofractionated radiotherapy (HFR)
3. Standard radiotherapy (SR)
"""
function construct_tmz_design()::TrialDesign
    N_pt = 93 + 98 + 100
    N_bm = 4
    N_tx = 3
    duration = 36.0 # months
    bm_generator() = begin
        x_bm = zeros(4)
        if rand() < 0.5
            x_bm[1:2] = [1, 0]
        else
            x_bm[1:2] = [0, 1]
        end
        if rand() < 0.5
            x_bm[3:4] = [1, 0]
        else
            x_bm[3:4] = [0, 1]
        end
        return x_bm
    end
    tx_generator() = begin
        x_tx = zeros(3)
        x_tx[rand(1:N_tx)] = 1
        return x_tx
    end
    return TrialDesign(N_pt, N_bm, N_tx, duration, bm_generator, tx_generator)
end

"""
Construct a fairly wide, albeit somewhat informative prior distribution.

Parameter space is (β_0, β_bm, β_tx, β_int, log_σ_u, log_α)
"""
function construct_prior(design::TrialDesign;
                         scale_β = 10.0,
                         scale_u = 10.0)
    N_bm = design.N_bm
    N_tx = design.N_tx
    N_int = N_bm * N_tx
    N_pred = 1 + N_bm + N_tx + N_int
    prior = Product(append!(Normal.(zeros(N_pred), scale_β),  # fixed-effects
                            [Normal(log(scale_u), log(10)),   # random-effects
                             Normal(0.0, 1.0)]))       # baseline Weibull shape
    return prior
end

"""
For regression testing againt R implementation.

R code to reproduce:

> library("survival")
> data("lung")
> coxph(Surv(time, status) ~ sex, data = lung)

Should produce

Call:
coxph(formula = Surv(time, status) ~ sex, data = lung)

       coef exp(coef) se(coef)      z       p
sex -0.5310    0.5880   0.1672 -3.176 0.00149

Likelihood ratio test=10.63  on 1 df, p=0.001111
n= 228, number of events= 165 

While the coxph function here:

julia> data = get_lung_data();
julia> coxph(data)

Should produce the same coefficient:

CoxModel{Float64}

Coefficients:
────────────────────────────────────────────
     Estimate  Std.Error   z value  Pr(>|z|)
────────────────────────────────────────────
x1  -0.531024   0.167179  -3.17638    0.0015
────────────────────────────────────────────
"""
function get_lung_data()::SurvivalData
    df = dropmissing(CSV.read("lung.csv", normalizenames=true))
    times = df.time
    # status == 1 -> right censored, status == 2 -> death observed
    observed = df.status .== 2
    predictors = reshape(df.sex, (nrow(df), 1))
    return SurvivalData(times, observed, predictors)
end
