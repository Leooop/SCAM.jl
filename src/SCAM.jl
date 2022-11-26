module SCAM

using Base: @kwdef
using StaticArrays
using IfElse
using Requires

function __init__()
    @require OrdinaryDiffEq="1dea7af3-3e70-54e6-95c3-0bf5283fa5ed" begin
        include("simulation_helpers.jl")
        include("dataset_simulation.jl")
    end
end

include("types.jl")
include("damage_functions.jl")
include("functions.jl")

export Model, Rheology, NumericalSetup, TriaxialSetup, DeformationControl, ConstantStrainRate, ConstantStress, ConstitutiveModel, LinearWeakening, AsymptoticWeakening, EnergyBasedWeakening
export ConstantYieldStress, CoulombYieldStress, StrainWeakenedCoulombYieldStress, SimpleCharlesLaw, InvariantsKICharlesLaw, PrincipalKICharlesLaw, MicroMechanicalCharlesLaw, Elasticity, IncompressibleElasticity, Geom2D, Geom3D
export Plasticity, MinViscosityThreshold, DamageThreshold
export update_derivatives, update_derivatives!
export init_variables, init_params, simulate

end#module
