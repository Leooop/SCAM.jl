module SCAM

using Base: @kwdef
using StaticArrays
using IfElse

using Reexport
#@reexport using DamagedShearBand: Rheology

include("types.jl")
include("damage_functions.jl")
include("functions.jl")
include("dataset_simulation.jl")
include("simulation_helpers.jl")

export Model, Rheology, NumericalSetup, TriaxialSetup, DeformationControl, ConstantStrainRate, ConstantStress, ConstitutiveModel, LinearWeakening, AsymptoticWeakening, EnergyBasedWeakening
export ConstantYieldStress, CoulombYieldStress, StrainWeakenedCoulombYieldStress, SimpleCharlesLaw, InvariantsKICharlesLaw, PrincipalKICharlesLaw, MicroMechanicalCharlesLaw, Elasticity, IncompressibleElasticity, Geom2D, Geom3D
export Plasticity, MinViscosityThreshold, DamageThreshold
export update_derivatives, update_derivatives!
export init_variables, init_params, simulate

end#module
