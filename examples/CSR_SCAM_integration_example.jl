using SCAM
using OrdinaryDiffEq
using Plots

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

tspan = (0.0, 550)
sol = simulate(model, tspan; 
        solver = Tsit5(), # ODE solver
        saveat = [], 
        abstol = 1e-6,
        reltol = 1e-4,
        maxiters = 1e5,
        Dᵢ=nothing, # if nothing D(t0) = D0
        Dmax=0.95,
        stop_at_peak = false,
        cb=nothing # whatever DiffEq Callback. If nothing, uses the appropriate callbacks
) 

ϵs = -sol.t.*ϵ̇.*100 # in %
σs = -sol[1,:]./1e6 # in MPa
Ds = sol[2,:]
plot(ϵs, σs,
    c=:black,
    lw=2,
    label="",
    xlabel = "axial strain (%)",
    ylabel = "axial stress (MPa)"
)
plot!(twinx(), ϵs, Ds,
    c = :firebrick,
    lw=2,
    label = "",
    ylabel= "damage"
)