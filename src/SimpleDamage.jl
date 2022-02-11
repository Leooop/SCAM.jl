module SimpleDamage

using Base: @kwdef
using StaticArrays
import DamagedShearBand as DSB

using Reexport
@reexport using DamagedShearBand: Rheology

include("types.jl")
include("damage_functions.jl")
include("functions.jl")
include("data_simulation.jl")

export Model, NumericalSetup, TriaxialSetup, DeformationControl, ConstantStrainRate, ConstantStress, ConstitutiveModel, LinearWeakening, AsymptoticWeakening, EnergyBasedWeakening
export ConstantYieldStress, CoulombYieldStress, SimpleCharlesLaw, InvariantsKICharlesLaw, PrincipalKICharlesLaw, MicroMechanicalCharlesLaw, Elasticity, IncompressibleElasticity, Geom2D, Geom3D
export Plasticity, MinViscosityThreshold, DamageThreshold
export update_derivatives, update_derivatives!
export init_variables, init_params

end#module
