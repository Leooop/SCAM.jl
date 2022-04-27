using DataFormatter
using DataFrames
using OrdinaryDiffEq
using NamedArrays

# functions defined to allow simulation of DataSets from DataFormatter.jl

function simulate_ponctual(p, data_params, target)
    if target == [:Δσ_peak]
        df_sol = get_peak_values(p, data_params)
    elseif target == [:ϵ̇c]
        df_sol = get_creep_ϵ̇_min(p,data_params)
    elseif target == [:ϵ_peak, :Δσ_peak]
        df_sol = get_peak_ϵ_and_Δσ(p, data_params)
    else
        error("targets not recognized for ponctual dataset")
    end
    return df_sol
end

function simulate_timeseries(p, data_params, target, time_vec=Float64[] ; stop_at_peak=false)
    if target == [:Δσ]
        df_sol = integrate_csr(p, data_params ; time_vec, stop_at_peak)
    elseif target == [:ϵ̇c]
        df_sol = get_creep_ϵ̇(p, data_params ; time_vec)
    end
    return df_sol
end

function get_peak_values(p, data)
    df = integrate_csr(p, data ; stop_at_peak=true)
    row_peak = df[argmin(df.Δσ_sim),:]
    return DataFrame(["Δσ_peak_sim", "D_peak_sim"] .=> [row_peak.Δσ_sim, row_peak.D_sim])
end

function get_peak_ϵ_and_Δσ(p, data)
    df = integrate_csr(p, data ; stop_at_peak=true)
    row_peak = df[end,:] #argmax(abs.(df.Δσ_sim))
    ϵ_peak = row_peak.t * data.ϵ̇
    Δσ_peak = row_peak.Δσ_sim
    D_peak = row_peak.D_sim
    return DataFrame(["ϵ_peak_sim", "Δσ_peak_sim", "D_peak_sim"] .=> [ϵ_peak, Δσ_peak, D_peak])
end

function get_creep_ϵ̇_min(p, data)
    df_creep = get_creep_ϵ̇(p,data)
    id_min = argmin(abs.(df_creep.ϵ̇c_sim))
    ϵ̇cmin = clamp(df_creep.ϵ̇c_sim[id_min],-1,-1e-16)
    Dc = df_creep.Dc_sim[id_min]
    return DataFrame(["ϵ̇c_sim", "Dc_sim"] .=> [ϵ̇cmin, Dc])
end
function get_creep_ϵ̇(p, data ; time_vec=Float64[])

    Δσc = data.Δσc
    sc = Δσ2sₐₓ(p.model.geom,Δσc)

    # define model
    model = build_model(p, data; control_type=ConstantStrainRate)
    
    #define derivative update functions
    update_func!(du,u,p,t) = update_derivatives!(du,u,p,t,model)

    # get initial conditions and tspan
    Dᵢ = (:Dᵢ ∈ names(p.mp)[1]) ? p.mp[:Dᵢ] : getp(:D₀,p.mp,p.mp_default)
    u0 = init_variables(model, Dᵢ)
    tspan = (0.0, p.model.ϵmax/data.ϵ̇_dev) ; @assert tspan[2] >= 0
    
    # create callbacks
    # conditionD(u,t,i) = u[2] - p.model.Dmax # if damage overshoots 1
    # affectD!(integrator) = (integrator.u[2] = p.model.Dmax) # set it to 1
    # cbD = ContinuousCallback(conditionD, affectD!, save_positions=(false,false)) # check that continuously

    condition_sc(u,t,i) = u[1] - sc # if axial dev stress equals target creep stress
    cbsc = ContinuousCallback(condition_sc, terminate! ; save_positions=(true,false)) # check that continuously

    conditionpeak(u,t,i) = (abs(i.u[1]) < abs(i.uprev[1])) # if axial stress decreases
    cbpeak = DiscreteCallback(conditionpeak, terminate! ; save_positions=(false,false))
    cb = CallbackSet(cbpeak, cbsc)


    # create ODEProblem and solve
    prob = ODEProblem(update_func!, convert.(eltype(p.mp), u0), tspan, p.model)
    sol = solve(prob, Tsit5();#Rosenbrock23();#Tsit5();
            abstol=p.solver.abstol,
            reltol=p.solver.reltol,
            callback = cb,
            maxiters=p.solver.maxiters,
            save_everystep=false)

    @assert length(sol.t) == 2
    sc_eff = sol[1,2]
    Dᵢc = sol[2,2]
    if abs(sc - 0.001*sc) <= abs(sc_eff) <= abs(sc + 0.001*sc)
        #@show Dᵢc
    else
        @warn "axial stress could not be integrated until data creep stress $sc, instead stoped at peak stress $sc_eff, proceeding with peak damage value $Dᵢc"
    end

    # INTEGRATE AT CONSTANT STRESS
    model_creep = build_model(p, data; control_type=ConstantStress)
    update_func_creep!(du,u,p,t) = update_derivatives!(du,u,p,t,model_creep)
    thousand_years = 3600*24*365*1e3

    u0 = [0.0, Dᵢc]
    prob_creep = ODEProblem(update_func_creep!, u0, (0.0,thousand_years), p.model)

    # stop at Dmax
    # conditionD(u,t,i) = u[2] - p.model.Dmax # if damage overshoots 1
    # affectD!(integrator) = (integrator.u[2] = p.model.Dmax) # set it to 1
    # cbD = ContinuousCallback(conditionD, affectD!, save_positions=(false,false)) # check that continuously

    conditionD(u,t,i) = u[2] - p.model.Dmax # if damage overshoots 1
    cbD = ContinuousCallback(conditionD, terminate!, save_positions=(false,false)) # check that continuously
    solc = solve(prob_creep, Tsit5() ; #KenCarp4(autodiff=false)
                    reltol=p.solver.abstol,
                    abstol=p.solver.reltol,
                    callback=cbD,
                    maxiters=p.solver.maxiters,
                    saveat=time_vec)

    D_vec = solc[2,:] 
    Ḋ_vec =  [Ḋ_func(sc_eff, D, p.model, model_creep) for D in D_vec]
    η_dam_vec = [η_dam_func(model_creep,D_vec[i],Ḋ_vec[i]) for i in eachindex(D_vec)]
    ϵ̇c_vec = sc_eff ./ (2 .* η_dam_vec)

    # append vectors if integration terminated without sampling all times :
    if !isempty(time_vec)
        len_solc = length(solc.t)
        len_ts = length(time_vec)
        if len_solc < len_ts
            missing_length = len_ts - len_solc
            missing_ϵ̇c = fill(-1.0, missing_length)
            missing_D = fill(p.model.Dmax, missing_length)
            append!(ϵ̇c_vec,missing_ϵ̇c)
            append!(D_vec,missing_D)
            @assert all(solc.t .≈ time_vec[begin:len_solc])
        end
        tvec = time_vec
    else
        tvec = solc.t
    end
    
    df_sol = DataFrame(["t", "ϵ̇c_sim", "Dc_sim"].=> [tvec, ϵ̇c_vec, D_vec])
    return df_sol
end

function integrate_csr(p, data ; time_vec=Float64[], stop_at_peak=false)
    # define model
    model = build_model(p, data ; control_type=ConstantStrainRate)
    
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

    # condition_last_time(u,t,i) = t > tspan[2]
    # cb_last_time = DiscreteCallback(condition_last_time, terminate! ; save_positions=(false,false))

    if stop_at_peak
        conditionpeak(u,t,i) = (abs(i.u[1]) < abs(i.uprev[1])) # if axial stress decreases
        savepos = isempty(time_vec) ? (true,true) : (false,false)
        cbpeak = DiscreteCallback(conditionpeak, terminate! ; save_positions=savepos)
        cb = CallbackSet(cbD, cbpeak)
    else
        cb = CallbackSet(cbD)
    end


    # create ODEProblem and solve
    prob = ODEProblem(update_func!, convert.(eltype(p.mp), u0), tspan, p.model)
    sol = solve(prob, AutoTsit5(Rosenbrock23()); #Rodas5(); #Tsit5();#Rosenbrock23();#Tsit5();
            abstol=p.solver.abstol,
            reltol=p.solver.reltol,
            callback = cb,
            maxiters=p.solver.maxiters,
            saveat=time_vec)

    # find s_peak and convert to Δσ
    s_vec = sol[1,:]
    Δσ_vec = sₐₓ2Δσ.(Ref(model),s_vec)
    df_sol = DataFrame(["t","Δσ_sim", "D_sim"].=> [sol.t, Δσ_vec, sol[2,:]])
    #@show(df_sol)
    return df_sol
end


function build_model(p, data ; control_type=ConstantStrainRate)
    mp = p.mp
    mpd = p.mp_default
    geom = p.model.geom
    DamType = p.model.damage_type
    pc = data.pc
    ϵ̇ = data.ϵ̇_dev
    if control_type == ConstantStrainRate
        control = ConstantStrainRate(ϵ̇)
    elseif control_type == ConstantStress
        Δσc = data.Δσc
        sc = Δσ2sₐₓ(p.model.geom,Δσc)
        control = ConstantStress(sc)
    end

    r = Rheology(; μ=getp(:μ,mp,mpd), ψ=getp(:ψ,mp,mpd), a=getp(:a,mp,mpd), D₀=getp(:D₀,mp,mpd), n=getp(:n,mp,mpd), K₁c=getp(:K₁c,mp,mpd), l̇₀=getp(:l̇₀,mp,mpd))
    plasticity = Plasticity(MinViscosityThreshold(), CoulombYieldStress(μ=getp(:μ,mp,mpd),C=0))
    cm = ConstitutiveModel(;weakening = LinearWeakening(γ=getp(:γ,mp,mpd)),
                       damage = DamType(r),
                       elasticity = IncompressibleElasticity(G=getp(:G,mp,mpd)), #3.146e9
                       plasticity)
    setup = TriaxialSetup(;geom,
                           control,
                           pc,
                           )
    return Model(;setup, cm)
end

getp(sym::Symbol,A::NamedArray,B::NamedArray) = haskey(A.dicts[1],sym) ? A[sym] : B[sym]