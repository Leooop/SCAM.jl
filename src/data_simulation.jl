using DataFormatter
using DataFrames
using OrdinaryDiffEq

# functions defined to allow simulation of DataSets from DataFormatter.jl

function simulate_ponctual(p, data, target)
    if any(target .∈ Ref((:Δσ, :Δσ_peak)))
        df_sim = get_peak_values(p, data)
    elseif any(target .∈ Ref((:ϵ̇_dev,)))
    end
    return sim_target
end

function simulate_timeseries(p, data_params, target, time_vec=nothing)
    if any(target .∈ Ref((:Δσ,)))
        df_sol = integrate_csr(p, data_params ; time_vec)
    elseif any(target .∈ Ref((:ϵ̇_dev,)))
        #model = build_model(mp, data.pc ; control = ConstantStress(data.Δσ), geom = p.geom)
    end
    return df_sol
end


function integrate_csr(p, data ; time_vec=nothing)
    # define model
    model = build_model(p, data)
    #define derivative update function
    update_func!(du,u,p,t) = update_derivatives!(du,u,p,t,model)

    # get initial conditions and tspan
    Dᵢ = (:Dᵢ ∈ names(p.mp)[1]) ? p.mp[:Dᵢ] : p.mp[:D₀]
    u0 = init_variables(model, Dᵢ)
    tspan = (0, p.model.ϵmax/data.ϵ̇_dev)
    saveat = !isnothing(time_vec) ? time_vec : []
    
    # create callbacks
    conditionD(u,t,i) = u[2] - p.model.Dmax # if damage overshoots 1
    affectD!(integrator) = (integrator.u[2] = p.model.Dmax) # set it to 1
    cbD = ContinuousCallback(conditionD,affectD!) # check that continuously
    conditionpeak(u,t,i) = (abs(i.u[1]) < abs(i.uprev[1])) # if axial stress decreases
    cbpeak = DiscreteCallback(conditionpeak,terminate! ; save_positions=(false,false))
    cb = CallbackSet(cbD,cbpeak)

    # create ODEProblem and solve
    prob = ODEProblem(update_func!, u0, tspan, p.model)
    sol = solve(prob, Tsit5();
            abstol=p.solver.abstol,
            reltol=p.solver.reltol,
            callback = cb,
            maxiters=p.solver.maxiters,
            saveat)

    # find s_peak and convert to Δσ
    s_vec = sol[1,:]
    Δσ_vec = sₐₓ2Δσ.(Ref(model),s_vec)
    return DataFrame(["Δσ_sim", "D_sim"].=> [Δσ_vec, sol[2,:]])
end

function build_model(p, data)
    mp = deepcopy(p.mp_default)
    params_names = names(p.mp)[1]
    for name in params_names
        val = getindex(p.mp,name)
        (name != :Dᵢ) && setindex!(mp,val,name)
    end
    geom = p.model.geom
    pc = data.pc
    ϵ̇ = data.ϵ̇_dev

    r = Rheology(; μ=mp[:μ], ψ=mp[:ψ], a=mp[:a], D₀=mp[:D₀], n=mp[:n], K₁c=mp[:K₁c], l̇₀=mp[:l̇₀])
    plasticity = Plasticity(MinViscosityThreshold(), CoulombYieldStress(μ=mp[:μ],C=0))
    cm = ConstitutiveModel(;weakening = LinearWeakening(γ=mp[:γ]),
                       damage = PrincipalKICharlesLaw(r),
                       elasticity = IncompressibleElasticity(G=mp[:G]), #3.146e9
                       plasticity)
    setup = TriaxialSetup(;geom,
                           control = ConstantStrainRate(ϵ̇),
                           pc,
                           )
    return Model(;setup, cm)
end


