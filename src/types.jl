
Maybe(T)=Union{T,Nothing}
# geometric types
abstract type Geom end
struct Geom2D <: Geom end
struct Geom3D <: Geom end


# RHEOLOGY type :
const N = Number
@kwdef struct Rheology{t1<:N,t2<:N,t3<:N,t4<:N,t5<:N,t6<:N,t7<:N,t8<:N,t9<:N}
  G::t1 = 30e9
  μ::t2 = 0.6 # Friction coef
  β::t3 = 0.1 # Correction factor
  K₁c::t4 = 1.74e6 # Critical stress intensity factor (Pa.m^(1/2))
  a::t5 = 1e-3 # Initial flaw size (m)
  ψ::t6 = atand(0.6)# crack angle to the principal stress (degree)
  D₀::t7 = 0.1# Initial flaw density
  n::t8 = 34.0 # Stress corrosion index
  l̇₀::t9 = 0.24 # Ref. crack growth rate (m/s)
end
const MicroMechanicalParameters = Rheology

function Rheology(r::Rheology, kw::NamedTuple)
  values_dict = Dict{Symbol,Float64}()
  for sym in propertynames(r)
    values_dict[sym] = isdefined(kw,sym) ? getproperty(kw,sym) : getproperty(r,sym)
  end
  return Rheology(; G = values_dict[:G],
                    μ = values_dict[:μ],
                    β = values_dict[:β],
                    K₁c = values_dict[:K₁c],
                    a = values_dict[:a],
                    ψ = values_dict[:ψ],
                    D₀ = values_dict[:D₀],
                    n = values_dict[:n],
                    l̇₀ = values_dict[:l̇₀])
end

Base.Broadcast.broadcastable(r::Rheology) = Ref(r)

function Base.show(io::IO, ::MIME"text/plain", r::Rheology)
  print(io, "Rheology instance with fields :\n",
  "\t├── G (shear modulus)                             : $(r.G)\n",
  "\t├── μ (flaws friction coefficient)                : $(r.μ)\n",
  "\t├── β (correction factor)                         : $(r.β)\n",
  "\t├── K₁c (fracture toughness)                      : $(r.K₁c)\n",
  "\t├── a (flaws radius)                              : $(r.a)\n",
  "\t├── ψ (flaws angle wrt σ₁ in degree)              : $(r.ψ)\n",
  "\t├── D₀ (initial damage, constrains flaws density) : $(r.D₀)\n",
  "\t├── n (stress corrosion index)                    : $(r.n)\n",
  "\t├── l̇₀ (ref crack growth rate in m/s)             : $(r.l̇₀)\n")
end


# Smoother types
abstract type Smoother end

"""
    SwitchSmoother <: Smoother

switch smoother based on tanh function

# Fields
- 'steepness': steepness of the transition 

"""
@kwdef struct SwitchSmoother <: Smoother
    x_one_percent::Float64 = 1e-3
end

"tanh smoother : goes from 0 to 1 when x goes from negative to positive, steepness is a field of the instanciated type"
function (ss::SwitchSmoother)(x) 
    k = atanh(2*(0.99-0.5))/ss.x_one_percent
    0.5 + 0.5*tanh(k * x)
end

"tanh smoother, switches from v1 to v2 when x goes from negative to positive, steepness is a field of the instanciated type"
(ss::SwitchSmoother)(v1,v2,x) = (1 - ss(x))*v1 + ss(x)*v2 

"tanh smoother, switches from v1 to v2 when x = (v2-v1)/v2 goes from negative to positive, steepness is a field of the instanciated type"
(ss::SwitchSmoother)(v1,v2) = (1 - ss(x))*v1 + ss(x)*v2 

# Numerical setup types
abstract type DeformationControl end
@kwdef struct ConstantStrainRate{tE<:Real} <: DeformationControl
    ϵ̇::tE = 1e-5
end
@kwdef struct ConstantStress{tS<:Real} <: DeformationControl
    s::tS = 1e8
end

abstract type NumericalSetup end
@kwdef struct TriaxialSetup{tG<:Geom,tC<:DeformationControl,tP<:Real} <: NumericalSetup 
    geom::tG = Geom3D() # geometric approximation of the problem
    control::tC = ConstantStrainRate()
    pc::tP = 1e5
end
control(ns::NumericalSetup) = ns.control

# weakening types
abstract type Weakening end
@kwdef struct LinearWeakening{T<:Real} <: Weakening
    γ::T = 1.0 # max total weakening ∈ [0,1]
end
@kwdef struct AsymptoticWeakening{T<:Real} <: Weakening
    γ::T = 1.0
end
struct EnergyBasedWeakening <: Weakening end

# yield stress types
abstract type YieldStress end
@kwdef struct ConstantYieldStress{T<:Real} <: YieldStress
    val::T = 0.0
end
(σy::ConstantYieldStress)() = σy.val
(σy::ConstantYieldStress)(p) = σy()

@kwdef struct CoulombYieldStress{tμ<:Real,tC<:Real} <: YieldStress 
    μ::tμ = 0.6
    C::tC = 0.0
end
(σy::CoulombYieldStress)(p) = (ϕ = atan(σy.μ) ; p*sin(ϕ) + σy.C*cos(ϕ))

@kwdef struct StrainWeakenedCoulombYieldStress{tμ<:Real,tC<:Real,tμ2<:Real,tC2<:Real,tϵpc<:Real} <: YieldStress 
    μ::tμ = 0.6
    C::tC = 0.0
    μweak::tμ2 = 0.6
    Cweak::tC2 = 0.0
    ϵₚcrit::tϵpc = 0.0
end
μeff(σy::StrainWeakenedCoulombYieldStress,ϵₚ) = abs(ϵₚ) < σy.ϵₚcrit ? σy.μweak + ((σy.ϵₚcrit-abs(ϵₚ))/σy.ϵₚcrit)*(σy.μ-σy.μweak) : σy.μweak
Ceff(σy::StrainWeakenedCoulombYieldStress,ϵₚ) = abs(ϵₚ) < σy.ϵₚcrit ? σy.Cweak + ((σy.ϵₚcrit-abs(ϵₚ))/σy.ϵₚcrit)*(σy.C-σy.Cweak) : σy.Cweak
function (σy::StrainWeakenedCoulombYieldStress)(p,ϵₚ)
    μ, C = μeff(σy,ϵₚ), Ceff(σy,ϵₚ)
    ϕ = atan(μ)
    p*sin(ϕ) + C*cos(ϕ)
end

@kwdef struct DamagedViscoPlasticYieldStress{tμ<:Real,tC<:Real} <: YieldStress 
    μ::tμ = 0.6
    C::tC = 0.0
end
(σy::DamagedViscoPlasticYieldStress)(p) = (ϕ = atan(σy.μ) ; p*sin(ϕ) + σy.C*cos(ϕ))

# helper func:
(σy::YieldStress)(p,ϵₚ) = σy(p)

# damage growth types
abstract type DamageGrowth end
@kwdef struct SimpleCharlesLaw{tTd<:Real,tS<:Real,tM<:Real,tY<:YieldStress} <: DamageGrowth
    Td::tTd = 4.5e6
    S::tS = 4.5e6
    m::tM = 4.0
    σy::tY = CoulombYieldStress() # instance of YieldStress subtype
end

abstract type MicroMechanicalCharlesLaw <: DamageGrowth end
"Damage mechanics where KI is formulated using stress invariants"
@kwdef struct InvariantsKICharlesLaw <: MicroMechanicalCharlesLaw
    r::Rheology = Rheology() # Type containing micromechanical parameters defined in DamagedShearBand.jl along with damage functions used to evaluate KI 
end
"Damage mechanics where KI is formulated using principal stresses, only works in the 2D case"
@kwdef struct PrincipalKICharlesLaw <: MicroMechanicalCharlesLaw
    r::Rheology = Rheology() # Type containing micromechanical parameters defined in DamagedShearBand.jl along with damage functions used to evaluate KI 
end



abstract type Elasticity end
@kwdef struct IncompressibleElasticity{tG<:Real} <: Elasticity
    G::tG = 4e9
end

# Plasticity related types 

abstract type PlasticityThreshold end
@kwdef struct MinViscosityThreshold{tS<:Maybe(SwitchSmoother)} <: PlasticityThreshold 
    smoother::tS = nothing
end


@kwdef struct DamageThreshold{tS<:Maybe(SwitchSmoother),T<:Real} <: PlasticityThreshold
    smoother::tS = nothing
    value::T = 0.99
end

@kwdef struct Plasticity{tC<:Maybe(PlasticityThreshold),tY<:YieldStress}
    threshold::tC = DamageThreshold(0.95)
    σy::tY = CoulombYieldStress() # instance of YieldStress subtype
end


# Constitutive model type
@kwdef struct ConstitutiveModel{tW<:Maybe(Weakening),tD<:Maybe(DamageGrowth),tE<:Elasticity,tP<:Maybe(Plasticity)} 
    weakening::tW = LinearWeakening()# elastic modulus weakening type
    damage::tD = SimpleCharlesLaw()# Damage growth type
    elasticity::tE = IncompressibleElasticity() # Elasticity
    plasticity::tP = nothing
end

# Model type
@kwdef struct Model{tCM<:ConstitutiveModel,tNS<:NumericalSetup}
    cm::tCM = ConstitutiveModel()
    setup::tNS = TriaxialSetup()
end
Model(m::Model,cm::ConstitutiveModel) = Model(cm,m.setup)
Model(m::Model,setup::NumericalSetup) = Model(m.cm,setup)

control(m::Model) = m.setup.control
#Model{tCM,tNS}(cm,setup) where {tCM,tNS} = Model(cm,setup)

# ALIASES
const M = Model
const NS = NumericalSetup
const TS = TriaxialSetup
const DC = DeformationControl
const CSR = ConstantStrainRate
const CS = ConstantStress
const CM = ConstitutiveModel
const W = Weakening
const LW = LinearWeakening
const AW = AsymptoticWeakening
const EBW = EnergyBasedWeakening
const ConstY = ConstantYieldStress
const CoulombY = CoulombYieldStress
const SWCoulombY = StrainWeakenedCoulombYieldStress
const DVPY = DamagedViscoPlasticYieldStress
const DG = DamageGrowth
const SD = SimpleCharlesLaw
const MD = MicroMechanicalCharlesLaw
const IKI = InvariantsKICharlesLaw
const PKI = PrincipalKICharlesLaw
const E = Elasticity
const IE = IncompressibleElasticity
const P = Plasticity
const MVT = MinViscosityThreshold
const DT = DamageThreshold