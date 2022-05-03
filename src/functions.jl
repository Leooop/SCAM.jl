##########################
#### HELPER FUNCTIONS ####
##########################

confining_pressure(s::NumericalSetup) = s.pc
confining_pressure(m::Model) = confining_pressure(m.setup)

plastic_yield_stress(plas::Plasticity,p,ϵₚ) = plas.σy(p,ϵₚ)
plastic_yield_stress(m::Model,p,ϵₚ) = plastic_yield_stress(m.cm.plasticity,p,ϵₚ)

ϵ2ϵdev(ϵax,E,ν,pc) = 2/3 * (1+ν)/E * (E*ϵax + (1 - 2ν)*pc)
ϵ̇2ϵ̇dev(ϵ̇,ν) = 2/3 * (1+ν) * ϵ̇

is3D(g::Geom) = g isa Geom3D ? true : false
is3D(setup::TriaxialSetup) = is3D(setup.geom) 
is3D(m::Model) = is3D(m.setup) 
# variables initializers : [s_ax,D] or [ϵ_ax, D] depending on the control variable type (ConstantStress or ConstantStrainRate)
# no strain weakening :
init_variables(::SD) = [0.0,0.0]
init_variables(dam::MD) = [0.0,dam.r.D₀]
init_variables(::Nothing) = [0.0,0.0]
init_variables(dam::MD, Dᵢ) = [0.0,max(Dᵢ,dam.r.D₀)]
init_variables(::Nothing, Dᵢ) = [0.0,0.0]
init_variables(dam::Maybe(MD), ::YieldStress) = init_variables(dam::MD)
init_variables(dam::Maybe(MD), ::YieldStress, Dᵢ) = init_variables(dam::MD, Dᵢ)

# with strain weakening :
init_variables(::SD, ::SWCoulombY) = [0.0,0.0,0.0]
init_variables(dam::MD, ::SWCoulombY) = [0.0, dam.r.D₀, 0.0]
init_variables(::Nothing, ::SWCoulombY) = [0.0, 0.0, 0.0]
init_variables(dam::MD, ::SWCoulombY, Dᵢ) = [0.0, max(Dᵢ,dam.r.D₀), 0.0]

init_variables(m::Model) = init_variables(m.cm.damage, m.cm.plasticity.σy)
init_variables(m::Model, Dᵢ) = init_variables(m.cm.damage, m.cm.plasticity.σy, Dᵢ)


const DMAX_DEFAULT = 0.95
init_params(::Model{<:CM{tW,tD,tE,Nothing}}    ; Dmax=DMAX_DEFAULT) where {tW,tD,tE} = (;Dmax)
init_params(::Model ; Dmax=DMAX_DEFAULT) where {tW,tD,tE} = (;Dmax)
# init_params(::Model{<:CM{tW,tD,tE,<:P{<:DT}}}  ; Dmax=DMAX_DEFAULT) where {tW,tD,tE} = (;Dmax)
# init_params(::Model{<:CM{tW,tD,tE,<:P{<:MVT}}} ; Dmax=DMAX_DEFAULT) where {tW,tD,tE} = (;Dmax)

# get shear and normal stress on a defect oriented ψ degree wrt σ1, negative in compression by convention
"shear stress from σ₁, σ₃ and ψ"
flaw_shear_stress(σ₁,σ₃,ψ) = 0.5 * (σ₃ - σ₁) * sind(2*ψ)
"normal stress from σ₁, σ₃ and ψ"
flaw_normal_stress(σ₁,σ₃,ψ) =  0.5 *((σ₃ + σ₁) + (σ₃ - σ₁)*cosd(2*ψ))

####################
#### INVARIANTS ####
####################
# geometric factor linking second invariants of the deviatoric stress or strain tensor and s_ax or ϵ'_ax in triaxial setup
geom_τ_factor(setup::TriaxialSetup) = is3D(setup) ? sqrt(3)/2 : 1.0
geom_τ_factor(model::M) = geom_τ_factor(model.setup)

# stress invariants for a triaxial setup,
"stress invariants τ and σₘ, sₐₓ negative in compression"
function stress_invariants(sₐₓ, setup::TriaxialSetup) 
    pc = confining_pressure(setup)
	τ_fact = geom_τ_factor(setup)
	if is3D(setup.geom) # all 3D cases
		τ = τ_fact * abs(sₐₓ)
    	σₘ = 0.5*sₐₓ - pc
	else # 2D plane strain incompressible or inplane compressibility only.
		τ = τ_fact * abs(sₐₓ)
    	σₘ = sₐₓ - pc
	end
	return τ, σₘ
end
stress_invariants(s, model::M{<:CM,<:TS}) = stress_invariants(s, model.setup) 

# conversions between axial
sₐₓ2σₐₓ(setup::TriaxialSetup,s) = (Δσ=sₐₓ2Δσ(setup,s) ; Δσ - setup.pc)
sₐₓ2σₐₓ(m::Model,s) = sₐₓ2σₐₓ(m.setup,s)

σₐₓ2sₐₓ(setup::TriaxialSetup{Geom3D},σₐₓ) = (2/3)*(σₐₓ + setup.pc)
σₐₓ2sₐₓ(setup::TriaxialSetup{Geom2D},σₐₓ) = 0.5*(σₐₓ + setup.pc)
σₐₓ2sₐₓ(m::Model,σₐₓ) = σₐₓ2sₐₓ(m.setup,σₐₓ)

sₐₓ2Δσ(geom::Geom,s) = is3D(geom) ? (3/2)*s : 2.0*s
sₐₓ2Δσ(setup::TriaxialSetup,s) = sₐₓ2Δσ(setup.geom,s)
sₐₓ2Δσ(model::M,s) = sₐₓ2Δσ(model.setup,s)

Δσ2sₐₓ(geom::Geom,Δσ) = is3D(geom) ? (2/3)*Δσ : 0.5*Δσ
Δσ2sₐₓ(setup::TriaxialSetup,Δσ) = Δσ2sₐₓ(setup.geom,Δσ)
Δσ2sₐₓ(model::M,Δσ) = Δσ2sₐₓ(model.setup,Δσ)

# ϵII for incompressible triaxial experiment as defined in SiStER (ϵII = sqrt(0.5*eij*eij))
ϵII_func(ϵ, model::M{<:CM{tW,tD,<:IE},<:TS}) where {tW,tD} = geom_τ_factor(model)*abs(ϵ)
ϵ̇II_func(ϵ̇, model::M{<:CM{tW,tD,<:IE},<:TS}) where {tW,tD} = ϵII_func(ϵ̇, model) # same for the time derivative 

#############################
#### WEAKENING FUNCTIONS ####
#############################
#cm.weakening.max is the residual fraction of G0 remaining at D=1
# Linear simple damage
weakening_func(D,cm::CM{<:LW,<:SD}) = (1-(1-cm.weakening.γ)*D)
weakening_derivative(D,cm::CM{<:LW,<:SD}) = -(1-cm.weakening.γ*one(D))
# Linear simple micromechanical damage
weakening_func(D,cm::CM{<:LW,<:MD}) = ((cm.weakening.γ-1)*D + 1 - cm.damage.r.D₀*cm.weakening.γ) / (1-cm.damage.r.D₀)
weakening_derivative(D,cm::CM{<:LW,<:MD}) = (cm.weakening.γ-1)/(1-cm.damage.r.D₀)
# asymptotic simple damage 
weakening_func(D,cm::CM{<:AW,<:SD}) = 1/(1+D)^cm.weakening.γ #TODO: TESTS
weakening_derivative(D,cm::CM{<:AW,<:SD}) = -cm.weakening.γ/(1+D)^(cm.weakening.γ+1) #TODO: TESTS
# asymptotic micromechanical damage untested
weakening_func(D,cm::CM{<:AW,<:MD}) = ((D<cm.damage.r.D₀)&&(D=cm.damage.r.D₀) ; ((1+cm.damage.r.D₀)/(1+D))^cm.weakening.γ) #TODO: TESTS
weakening_derivative(D,cm::CM{<:AW,<:MD}) = -cm.weakening.γ*(1+cm.damage.r.D₀)^(cm.weakening.γ)/(1+D)^(cm.weakening.γ+1) #TODO: TESTS


weakening_func(D,m::Model) = weakening_func(D,m.cm)
weakening_derivative(D,m::Model) = weakening_derivative(D,m.cm)

# when material in non weakenable
weakening_func(D,m::M{<:CM{Nothing}}) = 1.0
weakening_derivative(D,m::M{<:CM{Nothing}}) = 0.0

function get_energy_based_G(D, τ, σ, s₁₂, model::M{<:CM{<:EBW,<:IKI}})
	r = model.cm.damage.r
	G = model.cm.elasticity.G
	A, B = DSB.compute_AB(r,D)
	A1, B1 = DSB.compute_A1B1(r,A,B)
	#return 1 / ( 1/(4G) * ( 1 + (B1/2) * (B1 + A1*σ/τ - A1*σ*s₁₂^2/τ^3) ) )
	return 1/(1 + B1^2/2 + A1*B1*σ /(2*τ))
end

#############
#### KI  ####
#############
"KI with yield stress from stress invariants (sign convention agnostic)"
function get_KI_from_s(s, D, model::M{<:CM{<:W,<:IKI}})
	r = model.cm.damage.r
    #D₀ = r.D₀ #τ_heal = 60
	τ, σₘ = stress_invariants(s, model)
	D = IfElse.ifelse(D < r.D₀, r.D₀, D)
	KI = DSB.compute_KI(r,σₘ,τ,D)
	return KI
end

"KI with yield stress from principal stresses (sign convention agnostic)"
function get_KI_from_s(s, D, model::M{<:CM{<:W,<:PKI}})
	r = model.cm.damage.r
	pc = confining_pressure(model)
	D = IfElse.ifelse(D < r.D₀, r.D₀, D)
	KI = compute_KI(r, sₐₓ2σₐₓ(model,s), -pc, D)
	return KI
end

#######################
#### DAMAGE GROWTH ####
#######################
Ḋ_func(s, D, p, model::M{<:CM{tW,Nothing}}) where {tW} = 0.0

"damage growth rate when damage type is SimpleCharlesLaw"
function Ḋ_func(s, D, p, model::M{<:CM{tW,<:SD}}) where {tW}
    d = model.cm.damage
    #τ_heal = 0
	τ, σₘ = stress_invariants(s, model)
	pin = -σₘ
    σy = model.cm.damage.σy(pin)
	Ḋ = IfElse.ifelse(((τ-σy)>0) & ((D < p.Dmax) | (tW<:AW)),
		1/(d.Td*(d.S)^d.m) * (τ-σy)^d.m,
		zero(D))
	return Ḋ #- (D-D₀)/τ_heal
end


"damage growth rate when damage type is MicroMechanicalCharlesLaw"
# function Ḋ_func(s, D, p, model::M{<:CM{tW,<:MD}}, τ_heal=Inf) where {tW}
# 	r = model.cm.damage.r
# 	(D >= p.Dmax) && (return -(p.Dmax-r.D₀)/τ_heal)
# 	(D <= 0) && (return 0.0)

# 	KI = get_KI_from_s(s, D, model)
# 	(KI <= 0) && (return -(D-r.D₀)/τ_heal)

# 	#(KI >= 50000) && @show D KI 

#   	dDdl = DSB.compute_dDdl(r,D) # damage derivative wrt crack length
# 	#(KI >= 50000) && @show dDdl 
	

#   	dldtmax = 2000.0
#   	dldt = min(r.l̇₀*(KI/r.K₁c)^(r.n),dldtmax)  #Vr cracks growth rate
#   	#@assert dDdl * dldt >= 0
# 	#(KI >= 50000) && @show dDdl*dldt ; println()
#   	return dDdl * dldt - (D-r.D₀)/τ_heal
# end

"damage growth rate when damage type is MicroMechanicalCharlesLaw"
function Ḋ_func(s, D, p, model::M{<:CM{tW,<:MD}}, τ_heal=Inf) where {tW}
	r = model.cm.damage.r
	D = IfElse.ifelse(D > p.Dmax, p.Dmax, D)
	D = IfElse.ifelse(D < r.D₀, r.D₀, D)

	KI = get_KI_from_s(s, D, model)
	KI = IfElse.ifelse(KI <= 0, zero(D), KI)
	#(KI >= 50000) && @show D KI 

  	dDdl = DSB.compute_dDdl(r,D) # damage derivative wrt crack length
	#(KI >= 50000) && @show dDdl 
	

  	dldtmax = 2000.0
  	dldt = min(r.l̇₀*(KI/r.K₁c)^(r.n),dldtmax)  #Vr cracks growth rate
	dldt = IfElse.ifelse(D == p.Dmax, zero(D), dldt)
	
  	#@assert dDdl * dldt >= 0
	#(KI >= 50000) && @show dDdl*dldt ; println()
  	return dDdl * dldt - (D-r.D₀)/τ_heal
end
#############################
#### VISCOSITY FUNCTIONS ####
#############################

function η_dam_func(model,D,Ḋ)
	G = model.cm.elasticity.G
	f = weakening_func(D,model)
	df = weakening_derivative(D,model)
	return - G * f^2 / (Ḋ * df)
end

# plastic viscosity for constant strain rate setup
function η_plas_func(model, s, ϵ̇, ϵₚ)
	τ, σₘ = stress_invariants(s, model)
	pin = -σₘ
    σy = plastic_yield_stress(model, pin, ϵₚ)
	ϵ̇II = ϵ̇II_func(ϵ̇, model)
	if τ >= σy
		η = abs(σy/(2*ϵ̇II))
	else
		η = 1e30
	end
	return η
end

# residual function from equilibrium between γ̇elast that corrects plastic viscosity and γ̇elast from current stress. (blue notebook p200-201)
# expressed in log(eta)

# fn(logη,ϵ̇,ϵ̇II,s,σy) = sqrt(0.75*(ϵ̇^2 - s*ϵ̇/exp10(logη)) + 3/16*(s^2/exp10(logη)^2)) - 0.5*(ϵ̇II - σy/exp10(logη))
# ∇fn(f::F, x; δx=1e-6) where F = (f(x+δx)-f(x-δx)) / 2δx
# function find_best_η(model, η, ϵ̇, s, σy ; maxiter=100, reltol=1e-8, α=5e2)
# 	ηlog = log10(η)
# 	ϵ̇II = ϵ̇II_func(ϵ̇, model)
# 	f(ηlog) = fn(ηlog,ϵ̇,ϵ̇II,s,σy)
# 	value_prev = f(ηlog)
# 	for i in 1:maxiter
# 		δη = - α * ∇fn(f,ηlog)
# 		ηlog += δη
# 		value = f(ηlog)
# 		abs(value - value_prev/value_prev) <= reltol && break
# 		value_prev = value
# 	end
# 	return exp10(η)
# end

# when no plasticity, viscosity is always linked to damage
get_viscosity(m::Model{<:CM{tW,tD,tE,Nothing}}, s, ϵ̇, D, Ḋ, ϵₚ, p) where {tW,tD,tE} = η_dam_func(m, D, Ḋ)
# when elasto-plastic only
get_viscosity(m::Model{<:CM{Nothing,Nothing,tE,tP}}, s, ϵ̇, D, Ḋ, ϵₚ, p) where {tE,tP} = η_plas_func(m, s, ϵ̇, ϵₚ)
# when plastic depends on the switch condition from damaged-elastic to elastic-plastic
function get_viscosity(m::Model{<:CM{tW,tD,tE,<:Plasticity{<:MinViscosityThreshold{Nothing}}}}, s, ϵ̇, D, Ḋ, ϵₚ, p) where {tW,tD,tE}
	#D₀ = m.cm.damage.r.D₀
	ηp =  η_plas_func(m, s, ϵ̇, ϵₚ)

	#τ, σₘ = stress_invariants(s, m)
	#pin = -σₘ
	#σy = plastic_yield_stress(m,pin)
	ηd = η_dam_func(m, D, Ḋ)

	#isaboveyield_τ = IfElse.ifelse(τ >= σy, true, false)
	#Dmincond = IfElse.ifelse(D >= D₀ + 0.5*(p.Dmax-D₀), true, false)
	η = ((ηd <= ηp) || (D >= p.Dmax)) ? ηp : ηd
	#@show D Ḋ ηd ηp η
	#println()

	return η
end

function get_viscosity(m::Model{<:CM{tW,tD,tE,<:Plasticity{<:DamageThreshold}}}, s, ϵ̇, D, Ḋ, ϵₚ, p) where {tW,tD,tE}
	return IfElse.ifelse(D >= m.cm.plasticity.threshold.value, η_plas_func(m, s, ϵ̇, ϵₚ), η_dam_func(m, D, Ḋ))
end


#is_plastic(m::Model{<:CM{tW,tD,tE,<:Plasticity{DamageThreshold}}},s,D) = (D >= m.cm.plasticity.threshold.value)
#is_plastic(m::Model{<:CM{tW,tD,tE,<:Plasticity{MinViscosityThreshold}}},s,D) = 

##########################
#### UPDATE FUNCTIONS ####
##########################

## INPLACE functions ##
#######################

### Constant strain rate :
"inplace derivative update function for constant strain rate test, variables are the deviatoric axial stress and the damage"
function update_derivatives!(du,u,p,t,model::M{tCM,<:TS{tG,<:CSR}}) where {tCM,tG}
    G = model.cm.elasticity.G
    ϵ̇ = model.setup.control.ϵ̇
	s, D = u
	Ḋ = Ḋ_func(s, D, p, model)
	η = get_viscosity(model, s, ϵ̇, D, Ḋ, 0.0, p)
	fw = weakening_func(D,model)
	du[1] = 2*G*fw * (ϵ̇ - s/(2η))
	du[2] = Ḋ 
    du
end

### Constant strain rate elasto-plastic only:
"inplace derivative update function for constant strain rate test, variables are the deviatoric axial stress and the damage"
function update_derivatives!(du,u,p,t,model::M{<:CM{tW,tD,tE,<:P{tPT,<:SWCoulombY}},<:TS{tG,<:CSR}}) where {tW,tD,tE,tPT,tG}
    G = model.cm.elasticity.G
    ϵ̇ = model.setup.control.ϵ̇
	s, D, ϵₚ = u
	Ḋ = Ḋ_func(s, D, p, model)
	η = get_viscosity(model, s, ϵ̇, D, Ḋ, ϵₚ, p)
	fw = weakening_func(D,model)
	du[1] = ṡ = 2*G*fw * (ϵ̇ - s/(2η))
	du[2] = Ḋ 
	τ̇, _ = stress_invariants(ṡ, model)
	du[3] = ϵ̇II_func(ϵ̇, model) - τ̇/(2*G*fw)
    du
end

### Constant stress :
"inplace derivative update function for constant strain rate test, variables are the axial strain and the damage"
function update_derivatives!(du,u,p,t,model::M{tCM,<:TS{tG,<:CS}}) where {tCM,tG}
    s = model.setup.control.s
	ϵ, D = u
	ϵ̇_prev = du[1]
	Ḋ = Ḋ_func(s, D, p, model)
	η = get_viscosity(model, s, ϵ̇_prev, D, Ḋ, 0.0, p)
	du[1] = s/(2η)
	du[2] = Ḋ
    du
end

## OUT OF PLACE functions ##
############################

"derivative update function for constant strain rate test, variables are the deviatoric axial stress and the damage"
function update_derivatives(u,p,t,model::M{tCM,<:TS{tG,<:CSR}}) where {tCM,tG}
    G = model.cm.elasticity.G
    ϵ̇ = model.setup.control.ϵ̇
	s, D = u
	Ḋ = Ḋ_func(s, D, p, model)
	η = get_viscosity(model, s, ϵ̇, D, Ḋ, 0.0, p)
	fw = weakening_func(D,model)
	ṡ = 2*G*fw * (ϵ̇ - s/(2η))
    return SA[ṡ, Ḋ]
end

"derivative update function for constant strain rate test, variables are the axial strain and the damage"
function update_derivatives(u,p,t,model::M{tCM,<:TS{tG,<:CS}}) where {tCM,tG}
    s = model.setup.control.s
	ϵ, D = u
	ϵ̇_prev = du[1]
	Ḋ = Ḋ_func(s, D, p, model)
	η = get_viscosity(model, s, ϵ̇_prev, D, Ḋ, 0.0, p)
	ϵ̇ = s/(2η)
    return SA[ϵ̇, Ḋ]
end