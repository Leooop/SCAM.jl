using SCAM
using OrdinaryDiffEq
using Plots

##########################
## CONSTANT STRAIN RATE ##
##########################

ϵ̇ = -1e-5

setup = TriaxialSetup(
    geom = Geom3D(),
    control = ConstantStrainRate(ϵ̇),
    pc = 50e6
) 

elast = IncompressibleElasticity(G = 30e9)
mmp = MicroMechanicalParameters(
        μ=0.7, 
        ψ=45, 
        a=0.5e-3, 
        D₀=0.2, 
        n=10, 
        K₁c=2e6, 
        l̇₀=1e-2
)
damage_growth = PrincipalKICharlesLaw(mmp)
weak = LinearWeakening(0.5)
cm = ConstitutiveModel(
    weakening = weak,
    damage = damage_growth,
    elasticity = elast,
    plasticity = nothing
)
model = Model(cm,setup)

tspan = (0.0, 530)
sol = simulate(model, tspan; 
        solver = Tsit5(), # ODE solver
        saveat = range(0, tspan[2]; length=500), 
        abstol = 1e-6,
        reltol = 1e-4,
        maxiters = 1e5,
        Dᵢ=nothing, 
        Dmax=0.95,
        stop_at_peak = false,
        cb=nothing 
) 

ϵs = -sol.t.*ϵ̇.*100 # in %
σs = sol[1,:]
Ds = sol[2,:]
plot(ϵs, -σs./1e6,
    c=:black,
    lw=2,
    label="",
    xlabel = "axial strain (%)",
    ylabel = "axial deviatoric stress (MPa)"
)
plot!(twinx(), ϵs, Ds,
    c = :firebrick,
    lw=2,
    label = "",
    ylabel= "damage"
)

#####################
## CONSTANT STRESS ##
#####################

creep_stress = -200e6
σs_to_peak = σs[1:findfirst(diff(σs).>= 0)]
id = argmin(abs.(σs_to_peak .- creep_stress))

σc = σs[id]
Dᵢc = Ds[id]


setup_creep = TriaxialSetup(
    geom = Geom3D(),
    control = ConstantStress(σc),
    pc = 50e6
) 
model_creep = Model(cm,setup_creep)

tspan = (0.0, 1500)
sol_creep = simulate(model_creep, tspan; 
        solver = Tsit5(), # ODE solver
        saveat = range(0, tspan[2]; length=500), 
        abstol = 1e-6,
        reltol = 1e-4,
        maxiters = 1e5,
        Dᵢ = Dᵢc, # if nothing D(t0) = D0
        Dmax=0.95,
        stop_at_peak = false,
        cb=nothing # whatever DiffEq Callback. If nothing, uses the appropriate callbacks
) 

ts_creep = sol_creep.t
ϵs_creep = -sol_creep[1,:].*100 # in %
Ds_creep = sol_creep[2,:]

plot(ts_creep, ϵs_creep,
    c=:black,
    lw=2,
    label="strain",
    legend =:topleft,
    xlabel = "time (s)",
    ylabel = "axial creep strain (%)"
)
plot!(twinx(), ts_creep, Ds_creep,
    c = :firebrick,
    lw=2,
    label = "",
    ylabel= "damage"
)