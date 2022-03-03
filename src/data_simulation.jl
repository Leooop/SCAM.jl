using DataFormatter
using DataFrames
using OrdinaryDiffEq
using NamedArrays

# functions defined to allow simulation of DataSets from DataFormatter.jl

function simulate_ponctual(p, data_params, target)
    if any(target .∈ Ref((:Δσ, :Δσ_peak)))
        df_sol = get_peak_values(p, data_params)
    elseif any(target .∈ Ref((:ϵ̇_dev,)))
    end
    return df_sol
end

function simulate_timeseries(p, data_params, target, time_vec=nothing)
    if any(target .∈ Ref((:Δσ,)))
        df_sol = integrate_csr(p, data_params ; time_vec)
    elseif any(target .∈ Ref((:ϵ̇_dev,)))
        #model = build_model(mp, data.pc ; control = ConstantStress(data.Δσ), geom = p.geom)
    end
    return df_sol
end

function get_peak_values(p, data)
    df = integrate_csr(p, data ; stop_at_peak=true)
    row_peak = df[argmin(df.Δσ_sim),:]
    return DataFrame(["Δσ_peak_sim", "D_peak_sim"] .=> [row_peak."Δσ_sim",row_peak."D_sim"])
end

function integrate_csr(p, data ; time_vec=Float64[], stop_at_peak=false)
    # define model
    model = build_model(p, data)
    
    #define derivative update functions
    update_func!(du,u,p,t) = update_derivatives!(du,u,p,t,model)

    # get initial conditions and tspan
    Dᵢ = (:Dᵢ ∈ names(p.mp)[1]) ? p.mp[:Dᵢ] : getp(:D₀,p.mp,p.mp_default)
    u0 = init_variables(model, Dᵢ)
    tspan = (0.0, p.model.ϵmax/data.ϵ̇_dev)
    
    # create callbacks
    conditionD(u,t,i) = u[2] - p.model.Dmax # if damage overshoots 1
    affectD!(integrator) = (integrator.u[2] = p.model.Dmax) # set it to 1
    cbD = ContinuousCallback(conditionD, affectD!, save_positions=(false,false)) # check that continuously

    condition_last_time(u,t,i) = t > tspan[2]
    cb_last_time = DiscreteCallback(condition_last_time, terminate! ; save_positions=(false,false))

    if stop_at_peak
        conditionpeak(u,t,i) = (abs(i.u[1]) < abs(i.uprev[1])) # if axial stress decreases
        cbpeak = DiscreteCallback(conditionpeak, terminate! ; save_positions=(true,true))
        cb = CallbackSet(cbD, cbpeak, cb_last_time)
    else
        cb = CallbackSet(cbD, cb_last_time)
    end


    # create ODEProblem and solve
    prob = ODEProblem(update_func!, convert.(eltype(p.mp), u0), tspan, p.model)
    sol = solve(prob, Tsit5();
            abstol=p.solver.abstol,
            reltol=p.solver.reltol,
            callback = cb,
            maxiters=p.solver.maxiters,
            saveat=time_vec)

    # find s_peak and convert to Δσ
    s_vec = sol[1,:]
    Δσ_vec = sₐₓ2Δσ.(Ref(model),s_vec)
    return DataFrame(["t","Δσ_sim", "D_sim"].=> [sol.t, Δσ_vec, sol[2,:]])
end

function build_model(p, data)
    mp = p.mp
    mpd = p.mp_default
    geom = p.model.geom
    pc = data.pc
    ϵ̇ = data.ϵ̇_dev

    r = Rheology(; μ=getp(:μ,mp,mpd), ψ=getp(:ψ,mp,mpd), a=getp(:a,mp,mpd), D₀=getp(:D₀,mp,mpd), n=getp(:n,mp,mpd), K₁c=getp(:K₁c,mp,mpd), l̇₀=getp(:l̇₀,mp,mpd))
    plasticity = Plasticity(MinViscosityThreshold(), CoulombYieldStress(μ=getp(:μ,mp,mpd),C=0))
    cm = ConstitutiveModel(;weakening = LinearWeakening(γ=getp(:γ,mp,mpd)),
                       damage = PrincipalKICharlesLaw(r),
                       elasticity = IncompressibleElasticity(G=getp(:G,mp,mpd)), #3.146e9
                       plasticity)
    setup = TriaxialSetup(;geom,
                           control = ConstantStrainRate(ϵ̇),
                           pc,
                           )
    return Model(;setup, cm)
end

getp(sym::Symbol,A::NamedArray,B::NamedArray) = haskey(A.dicts[1],sym) ? A[sym] : B[sym]