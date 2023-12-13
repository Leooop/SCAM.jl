##########################
#### HELPER FUNCTIONS ####
##########################

confining_pressure(s::NumericalSetup) = s.pc
confining_pressure(m::Model) = confining_pressure(m.setup)

plastic_yield_stress(plas::Plasticity, p, ϵₚ) = plas.σy(p, ϵₚ)
plastic_yield_stress(m::Model, p, ϵₚ) = plastic_yield_stress(m.cm.plasticity, p, ϵₚ)

ϵ2ϵdev(ϵax, E, ν, pc) = 2 / 3 * (1 + ν) / E * (E * ϵax + (1 - 2ν) * pc)
ϵ̇2ϵ̇dev(ϵ̇, ν) = 2 / 3 * (1 + ν) * ϵ̇

is3D(g::Geom) = g isa Geom3D ? true : false
is3D(setup::TriaxialSetup) = is3D(setup.geom)
is3D(m::Model) = is3D(m.setup)
# variables initializers : [s_ax,D] or [ϵ_ax, D] depending on the control variable type (ConstantStress or ConstantStrainRate)
# no strain weakening :
init_variables(::SD) = [0.0, 0.0]
init_variables(dam::MD) = [0.0, dam.r.D₀]
init_variables(::Nothing) = [0.0, 0.0]
init_variables(dam::MD, Dᵢ) = [0.0, max(Dᵢ, dam.r.D₀)]
init_variables(::Nothing, Dᵢ) = [0.0, 0.0]
init_variables(dam::Maybe(MD), ::Maybe(YieldStress)) = init_variables(dam::MD)
init_variables(dam::Maybe(MD), ::Maybe(YieldStress), Dᵢ) = init_variables(dam::MD, Dᵢ)

# with strain weakening :
init_variables(::SD, ::SWCoulombY) = [0.0, 0.0, 0.0]
init_variables(dam::MD, ::SWCoulombY) = [0.0, dam.r.D₀, 0.0]
init_variables(::Nothing, ::SWCoulombY) = [0.0, 0.0, 0.0]
init_variables(dam::MD, ::SWCoulombY, Dᵢ) = [0.0, max(Dᵢ, dam.r.D₀), 0.0]

init_variables(m::Model) = init_variables(m.cm.damage, m.cm.plasticity.σy)
init_variables(m::Model, Dᵢ) = init_variables(m.cm.damage, m.cm.plasticity.σy, Dᵢ)
init_variables(m::Model{<:CM{<:W,<:DG,<:E,Nothing}}) = init_variables(m.cm.damage, nothing)
init_variables(m::Model{<:CM{<:W,<:DG,<:E,Nothing}}, Dᵢ) = init_variables(m.cm.damage, nothing, Dᵢ)


const DMAX_DEFAULT = 0.95
init_params(::Model{<:CM{tW,tD,tE,Nothing}}; Dmax=DMAX_DEFAULT) where {tW,tD,tE} = (; Dmax)
init_params(::Model; Dmax=DMAX_DEFAULT) = (; Dmax)
# init_params(::Model{<:CM{tW,tD,tE,<:P{<:DT}}}  ; Dmax=DMAX_DEFAULT) where {tW,tD,tE} = (;Dmax)
# init_params(::Model{<:CM{tW,tD,tE,<:P{<:MVT}}} ; Dmax=DMAX_DEFAULT) where {tW,tD,tE} = (;Dmax)

# get shear and normal stress on a defect oriented ψ degree wrt σ1, negative in compression by convention
"shear stress from σ₁, σ₃ and ψ"
flaw_shear_stress(σ₁, σ₃, ψ) = 0.5 * (σ₃ - σ₁) * sind(2 * ψ)
"normal stress from σ₁, σ₃ and ψ"
flaw_normal_stress(σ₁, σ₃, ψ) = 0.5 * ((σ₃ + σ₁) + (σ₃ - σ₁) * cosd(2 * ψ))

####################
#### INVARIANTS ####
####################
# geometric factor linking second invariants of the deviatoric stress or strain tensor and s_ax or ϵ'_ax in triaxial setup
geom_τ_factor(setup::TriaxialSetup) = is3D(setup) ? sqrt(3) / 2 : 1.0
geom_τ_factor(model::M) = geom_τ_factor(model.setup)

# stress invariants for a triaxial setup,
"stress invariants τ and σₘ, sₐₓ negative in compression"
function stress_invariants(sₐₓ, setup::TriaxialSetup)
  pc = confining_pressure(setup)
  τ_fact = geom_τ_factor(setup)
  if is3D(setup.geom) # all 3D cases
    τ = τ_fact * abs(sₐₓ)
    σₘ = 0.5 * sₐₓ - pc
  else # 2D plane strain incompressible or inplane compressibility only.
    τ = τ_fact * abs(sₐₓ)
    σₘ = sₐₓ - pc
  end
  return τ, σₘ
end
stress_invariants(s, model::M{<:CM,<:TS}) = stress_invariants(s, model.setup)

# conversions between axial
sₐₓ2σₐₓ(setup::TriaxialSetup, s) = (Δσ = sₐₓ2Δσ(setup, s);
Δσ - setup.pc)
sₐₓ2σₐₓ(m::Model, s) = sₐₓ2σₐₓ(m.setup, s)

σₐₓ2sₐₓ(setup::TriaxialSetup{Geom3D}, σₐₓ) = (2 / 3) * (σₐₓ + setup.pc)
σₐₓ2sₐₓ(setup::TriaxialSetup{Geom2D}, σₐₓ) = 0.5 * (σₐₓ + setup.pc)
σₐₓ2sₐₓ(m::Model, σₐₓ) = σₐₓ2sₐₓ(m.setup, σₐₓ)

sₐₓ2Δσ(geom::Geom, s) = is3D(geom) ? (3 / 2) * s : 2.0 * s
sₐₓ2Δσ(setup::TriaxialSetup, s) = sₐₓ2Δσ(setup.geom, s)
sₐₓ2Δσ(model::M, s) = sₐₓ2Δσ(model.setup, s)

Δσ2sₐₓ(geom::Geom, Δσ) = is3D(geom) ? (2 / 3) * Δσ : 0.5 * Δσ
Δσ2sₐₓ(setup::TriaxialSetup, Δσ) = Δσ2sₐₓ(setup.geom, Δσ)
Δσ2sₐₓ(model::M, Δσ) = Δσ2sₐₓ(model.setup, Δσ)

# ϵII for incompressible triaxial experiment as defined in SiStER (ϵII = sqrt(0.5*eij*eij))
ϵII_func(ϵ, model::M{<:CM{tW,tD,<:IE},<:TS}) where {tW,tD} = geom_τ_factor(model) * abs(ϵ)
ϵ̇II_func(ϵ̇, model::M{<:CM{tW,tD,<:IE},<:TS}) where {tW,tD} = ϵII_func(ϵ̇, model) # same for the time derivative 

#############################
#### WEAKENING FUNCTIONS ####
#############################
#cm.weakening.max is the residual fraction of G0 remaining at D=1
# Linear simple damage
weakening_func(D, cm::CM{<:LW,<:SD}) = (1 - (1 - cm.weakening.γ) * D)
weakening_derivative(D, cm::CM{<:LW,<:SD}) = -(1 - cm.weakening.γ * one(D))
# Linear simple micromechanical damage
weakening_func(D, cm::CM{<:LW,<:MD}) = ((cm.weakening.γ - 1) * D + 1 - cm.damage.r.D₀ * cm.weakening.γ) / (1 - cm.damage.r.D₀)
weakening_derivative(D, cm::CM{<:LW,<:MD}) = (cm.weakening.γ - 1) / (1 - cm.damage.r.D₀)
# asymptotic simple damage 
weakening_func(D, cm::CM{<:AW,<:SD}) = 1 / (1 + D)^cm.weakening.γ #TODO: TESTS
weakening_derivative(D, cm::CM{<:AW,<:SD}) = -cm.weakening.γ / (1 + D)^(cm.weakening.γ + 1) #TODO: TESTS
# asymptotic micromechanical damage untested
weakening_func(D, cm::CM{<:AW,<:MD}) = ((D < cm.damage.r.D₀) && (D = cm.damage.r.D₀); ((1 + cm.damage.r.D₀) / (1 + D))^cm.weakening.γ) #TODO: TESTS
weakening_derivative(D, cm::CM{<:AW,<:MD}) = -cm.weakening.γ * (1 + cm.damage.r.D₀)^(cm.weakening.γ) / (1 + D)^(cm.weakening.γ + 1) #TODO: TESTS

#test weakening with G(Di) = G0
#weakening_func(D,cm::CM{<:LW,<:MD}) = ((cm.weakening.γ-1)*D + 1 - cm.damage.r.D₀*cm.weakening.γ) / (1-cm.damage.r.D₀)
#weakening_derivative(D,cm::CM{<:LW,<:MD}) = (cm.weakening.γ-1)/(1-cm.damage.r.D₀)


weakening_func(D, m::Model) = weakening_func(D, m.cm)
weakening_derivative(D, m::Model) = weakening_derivative(D, m.cm)

# when material in non weakenable
weakening_func(D, m::M{<:CM{Nothing}}) = 1.0
weakening_derivative(D, m::M{<:CM{Nothing}}) = 0.0

function get_energy_based_G(ϵ, D, model) #::M{<:CM{<:EBW,<:IKI}}
  r = model.cm.damage.r
  G = model.cm.elasticity.G
  A1, B1 = compute_A1B1_alt(r, D)
  Γ = compute_Γ(r, A1, B1)
  ϵkk = 0.0 # TODO : explore what happens with damage induced dilatancy
  eII = is3D(model) ? sqrt(3) * abs(ϵ) : 2 * abs(ϵ1) # strain invariant from harsha = √(2 eij eij)
  return G / Γ * (3 * (1 - 2r.ν) / (2 * (1 + r.ν)) + A1^2 / 2 - A1 * B1 * ϵkk / 2eII) # Cμ from Bhat, == G in 0D
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
  KI = compute_KI_invariants(r, σₘ, τ, D)
  return KI
end

"KI with yield stress from principal stresses (sign convention agnostic)"
function get_KI_from_s(s, D, model::M{<:CM{<:W,<:PKI}})
  r = model.cm.damage.r
  pc = confining_pressure(model)
  D = IfElse.ifelse(D < r.D₀, r.D₀, D)
  KI = compute_KI(r, sₐₓ2σₐₓ(model, s), -pc, D)
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
  Ḋ = IfElse.ifelse(((τ - σy) > 0) & ((D < p.Dmax) | (tW <: AW)),
    1 / (d.Td * (d.S)^d.m) * (τ - σy)^d.m,
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

  dDdl = compute_dDdl(r, D) # damage derivative wrt crack length
  #(KI >= 50000) && @show dDdl 


  dldtmax = 2000.0
  dldt = min(r.l̇₀ * (KI / r.K₁c)^(r.n), dldtmax)  #Vr cracks growth rate
  dldt = IfElse.ifelse(D == p.Dmax, zero(D), dldt)

  #@assert dDdl * dldt >= 0
  #(KI >= 50000) && @show dDdl*dldt ; println()
  return dDdl * dldt - (D - r.D₀) / τ_heal
end
#############################
#### VISCOSITY FUNCTIONS ####
#############################

function η_dam_func(model, D, Ḋ)
  G = model.cm.elasticity.G
  f = weakening_func(D, model)
  df = weakening_derivative(D, model)
  return -G * f^2 / (Ḋ * df)
end


# plastic viscosity for constant strain rate setup
function η_plas_func(model, s, ϵ̇, ϵₚ)
  Δσ = sₐₓ2Δσ(model, s)
  pc = confining_pressure(model)

  # radius and center of mohr coulomb circle : 0.5*abs(Δσ) >= sin(ϕ)(0.5*abs(Δσ) + pc) + cos(ϕ)*C 
  τmc = 0.5 * abs(Δσ)
  σmc = 0.5Δσ - pc

  σy = plastic_yield_stress(model, -σmc, ϵₚ)

  if τmc >= σy
    ϵ̇_eff = is3D(model) ? abs(0.75 * ϵ̇) : abs(ϵ̇) # sqrt((ϵ11 - ϵ22)^2 + 4ϵ12^2) (2D Mohr-Coulomb)
    η = abs(0.5 * σy / ϵ̇_eff)
  else
    η = Inf
  end

  # if model.cm.plasticity.σy isa StrainWeakenedCoulombYieldStress
  # 	@show τmc sy η
  # 	println()
  # end

  return η
end

# when no plasticity, viscosity is always linked to damage
get_viscosity(m::Model{<:CM{tW,tD,tE,Nothing}}, s, ϵ̇, D, Ḋ, ϵₚ, p) where {tW,tD,tE} = η_dam_func(m, D, Ḋ)

# when elasto-plastic only
get_viscosity(m::Model{<:CM{Nothing,Nothing,tE,tP}}, s, ϵ̇, D, Ḋ, ϵₚ, p) where {tE,tP} = η_plas_func(m, s, ϵ̇, ϵₚ)

# when plastic depends on the switch condition from damaged-elastic to elastic-plastic
function get_viscosity(m::Model{<:CM{tW,tD,tE,<:Plasticity{<:MinViscosityThreshold{Nothing}}}}, s, ϵ̇, D, Ḋ, ϵₚ, p) where {tW,tD,tE}

  ηp = η_plas_func(m, s, ϵ̇, ϵₚ)
  ηd = η_dam_func(m, D, Ḋ)

  # condition on interaction of cracks
  if D < p.Dmax - 1e-4
    dKI = get_KI_from_s(s, D + 1e-4, m) - get_KI_from_s(s, D, m)
    interaction_is_dominant = dKI > 0
  else
    interaction_is_dominant = true
  end

  if D >= p.Dmax
    η = ηp
  elseif (ηd <= ηp) & interaction_is_dominant
    η = ηp
  else
    η = ηd
  end
  #η = ((ηd <= ηp) || (D >= p.Dmax)) ? ηp : ηd
  #println("D = ", D, ", ηd = ", ηd, ", η = ", η)

  return η
end

function get_viscosity(m::Model{<:CM{tW,tD,tE,<:Plasticity{Nothing,<:DamagedViscoPlasticYieldStress}}}, s, ϵ̇, D, Ḋ, ϵₚ, p) where {tW,tD,tE}

  ηp = η_plas_func(m, s, ϵ̇, ϵₚ)
  ηd = η_dam_func(m, D, Ḋ)

  (D >= p.Dmax) * ηd

  return ηp + ηd
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
function update_derivatives!(du, u, p, t, model::M{tCM,<:TS{tG,<:CSR}}) where {tCM,tG}
  G = model.cm.elasticity.G
  ϵ̇ = model.setup.control.ϵ̇
  s, D = u
  Ḋ = Ḋ_func(s, D, p, model)
  η = get_viscosity(model, s, ϵ̇, D, Ḋ, 0.0, p)
  fw = weakening_func(D, model)
  du[1] = 2 * G * fw * (ϵ̇ - s / (2η))
  du[2] = Ḋ
  du
end

### Constant strain rate elasto-plastic only:
"inplace derivative update function for constant strain rate test, variables are the deviatoric axial stress and the damage"
function update_derivatives!(du, u, p, t, model::M{<:CM{tW,tD,tE,<:P{tPT,<:SWCoulombY}},<:TS{tG,<:CSR}}) where {tW,tD,tE,tPT,tG}
  G = model.cm.elasticity.G
  ϵ̇ = model.setup.control.ϵ̇
  s, D, ϵₚ = u
  Ḋ = Ḋ_func(s, D, p, model)
  η = get_viscosity(model, s, ϵ̇, D, Ḋ, ϵₚ, p)
  fw = weakening_func(D, model)
  du[1] = ṡ = 2 * G * fw * (ϵ̇ - s / (2η))
  du[2] = Ḋ
  du[3] = ϵ̇ - ṡ / (2 * G * fw)
  du
end

### Constant stress :
"inplace derivative update function for constant stress test, variables are the axial strain and the damage"
function update_derivatives!(du, u, p, t, model::M{tCM,<:TS{tG,<:CS}}) where {tCM,tG}
  s = model.setup.control.s
  ϵ, D = u
  Ḋ = Ḋ_func(s, D, p, model)
  #η = get_viscosity(model, s, 0.0, D, Ḋ, 0.0, p)
  ηd = η_dam_func(model, D, Ḋ)
  du[1] = s / (2ηd)
  du[2] = Ḋ
  du
end

geom_ϵII_factor(setup::TriaxialSetup) = is3D(setup) ? 0.75 : 1.0
geom_ϵII_factor(model::M) = geom_ϵII_factor(model.setup)

"inplace derivative update function for constant strain rate test, variables are the axial strain and the damage"
# function update_derivatives!(du,u,p,t,model::M{<:CM{tW,tD,tE,<:P{tPT,<:DVPY}},<:TS{tG,<:CS}}) where {tW,tD,tE,tPT,tG}
#     s = model.setup.control.s
# 	ϵ, D = u
# 	Ḋ = Ḋ_func(s, D, p, model)
# 	η = get_viscosity(model, s, 0.0, D, Ḋ, 0.0, p)

# 	Δσ = sₐₓ2Δσ(model,s)
# 	pc = confining_pressure(model)
# 	# radius and center of mohr coulomb circle : 0.5*abs(Δσ) >= sin(ϕ)(0.5*abs(Δσ) + pc) + cos(ϕ)*C 
# 	τmc = 0.5*abs(Δσ)
# 	σmc = 0.5Δσ - pc

# 	ηd = η_dam_func(model,D,Ḋ)
# 	# TODO :σy = plastic_yield_stress(model, -σmc, 0.0) + 2*

# 	if τmc >= σy
# 		du[1] = (s - (sign(s)*σy/geom_ϵII_factor(model)))/(2ηd)
# 	else
# 		du[1] = s/(2ηd)
# 	end
# 	du[2] = Ḋ
#     du
# end

## OUT OF PLACE functions ##
############################

"derivative update function for constant strain rate test, variables are the deviatoric axial stress and the damage"
function update_derivatives(u, p, t, model::M{tCM,<:TS{tG,<:CSR}}) where {tCM,tG}
  G = model.cm.elasticity.G
  ϵ̇ = model.setup.control.ϵ̇
  s, D = u
  Ḋ = Ḋ_func(s, D, p, model)
  η = get_viscosity(model, s, ϵ̇, D, Ḋ, 0.0, p)
  fw = weakening_func(D, model)
  ṡ = 2 * G * fw * (ϵ̇ - s / (2η))
  return SA[ṡ, Ḋ]
end

"derivative update function for constant strain rate test, variables are the axial strain and the damage"
function update_derivatives(u, p, t, model::M{tCM,<:TS{tG,<:CS}}) where {tCM,tG}
  s = model.setup.control.s
  ϵ, D = u
  ϵ̇_prev = du[1]
  Ḋ = Ḋ_func(s, D, p, model)
  η = get_viscosity(model, s, ϵ̇_prev, D, Ḋ, 0.0, p)
  ϵ̇ = s / (2η)
  return SA[ϵ̇, Ḋ]
end