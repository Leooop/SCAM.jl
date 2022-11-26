using OrdinaryDiffEq

function simulate(model, tspan; 
        solver = Tsit5(),
        saveat = [],
        abstol = 1e-6,
        reltol = 1e-4,
        maxiters = 1e5,
        Dᵢ=nothing,
        Dmax=0.95,
        stop_at_peak = false,
        cb=nothing)

    isnothing(Dᵢ) &&  (Dᵢ = model.cm.damage.r.D₀)

    update!(du,u,p,t) = update_derivatives!(du,u,p,t,model)
    p = init_params(model ; Dmax)
    prob = ODEProblem(update!, init_variables(model,Dᵢ), tspan, p)

    cb = !isnothing(cb) ? cb : get_callbacks(model, p ; stop_at_peak)
    sol = solve(prob, solver;#Tsit5() ;#Tsit5(); SSPRK83(), KenCarp4(autodiff=false)
            saveat,
            abstol,
            reltol,
            maxiters,
            callback=cb)
    return sol
end

function damage_limiter_callback(m, p)
    condition(u,t,integrator) = u[2] - p.Dmax # if damage overshoots 1
    affect!(integrator) = (integrator.u[2] = p.Dmax) # set it to 1
    ContinuousCallback(condition,affect!) # check that continuously
end

function terminate_at_Dmax_callback(m, p)
    condition(u,t,integrator) = u[2] - p.Dmax # if damage overshoots 1
    ContinuousCallback(condition,terminate!) # check that continuously
end

function stop_at_peak_callback(m, p)
    D₀ = m.cm.damage.r.D₀
    conditionpeak(u,t,i) = (abs(u[1]) < abs(i.uprev[1])) && (u[2] > D₀+0.5*(p.Dmax-D₀)) # if axial stress decreases
    DiscreteCallback(conditionpeak, terminate! ; save_positions=(true,false))
end

function get_callbacks(model, p ; stop_at_peak=false)
    stop_at_peak && return stop_at_peak_callback(model, p)
    simtype = typeof(control(model))
    if simtype <: ConstantStrainRate
        cb = isnothing(model.cm.plasticity) ? terminate_at_Dmax_callback(model, p) : damage_limiter_callback(model, p)
    elseif simtype <: ConstantStress
        cb = terminate_at_Dmax_callback(model, p)
    end
    return cb
end