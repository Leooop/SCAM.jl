
Maybe(T)=Union{T,Nothing}
# geometric types
abstract type Geom end
struct Geom2D <: Geom end
struct Geom3D <: Geom end

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
    r::Rheology = DSB.Rheology() # Type containing micromechanical parameters defined in DamagedShearBand.jl along with damage functions used to evaluate KI 
end
"Damage mechanics where KI is formulated using principal stresses, only works in the 2D case"
@kwdef struct PrincipalKICharlesLaw <: MicroMechanicalCharlesLaw
    r::Rheology = DSB.Rheology() # Type containing micromechanical parameters defined in DamagedShearBand.jl along with damage functions used to evaluate KI 
end



abstract type Elasticity end
@kwdef struct IncompressibleElasticity{tG<:Real} <: Elasticity
    G::tG = 4e9
end

# Plasticity related types 
#TODO : !!!!!!!!!
abstract type PlasticityThreshold end
@kwdef struct MinViscosityThreshold{tS<:Maybe(SwitchSmoother)} <: PlasticityThreshold 
    smoother::tS = nothing
end


@kwdef struct DamageThreshold{tS<:Maybe(SwitchSmoother),T<:Real} <: PlasticityThreshold
    smoother::tS = nothing
    value::T = 0.99
end

@kwdef struct Plasticity{tC<:PlasticityThreshold,tY<:YieldStress}
    threshold::tC = DamageThreshold(0.95)
    σy::tY = CoulombYieldStress() # instance of YieldStress subtype
end


# Constitutive model type
@kwdef struct ConstitutiveModel{tW<:Weakening,tD<:DamageGrowth,tE<:Elasticity,tP<:Maybe(Plasticity)} 
    weakening::tW = LinearWeakening()# elastic modulus weakening type
    damage::tD = SimpleCharlesLaw()# Damage growth type
    elasticity::tE = IncompressibleElasticity() # Elasticity
    plasticity::tP = nothing
end

# Model type
@kwdef struct Model{tCM<:ConstitutiveModel,tNS<:NumericalSetup}
    cm::tCM = ConstitutiveModel()
    setup::tNS = TriaxialSetup()
    function Model(cm::tCM,setup::tNS) where {tCM,tNS}
        if (cm.damage isa PrincipalKICharlesLaw) && is3D(setup.geom)
            @error "geometry can only be 2D when using `PrincipalKICharlesLaw` damage growth law"
        end
        new{tCM,tNS}(cm,setup)
    end
end
Model{tCM,tNS}(cm,setup) where {tCM,tNS} = Model(cm,setup)

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
const SD = SimpleCharlesLaw
const MD = MicroMechanicalCharlesLaw
const IKI = InvariantsKICharlesLaw
const PKI = PrincipalKICharlesLaw
const IE = IncompressibleElasticity
const P = Plasticity
const MVT = MinViscosityThreshold
const DT = DamageThreshold