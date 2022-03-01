using SimpleDamage ; import SimpleDamage as SD
using NamedArrays
import DamagedShearBand as DSB
using Test

pp =   [8.4e9, 0.8, 0.55, 2e-3, 0.4, 20.0, 1.014e6, 0.001, 0.4]
pp_names = [:G,   :γ,   :μ,   :a,  :D₀, :n,   :K₁c,   :l̇₀,   :Dᵢ]
pn = NamedArray(pp,pp_names)
r = DSB.Rheology(μ=pn[:μ], ψ=45, a=pn[:a], D₀=pn[:D₀], n=pn[:n], K₁c=pn[:K₁c], l̇₀=pn[:l̇₀])
pc = 1e6
ϵ̇ = 1e-5
s = 1e8
C = 0.0
τy_const = 2e6

plas_coulomb_minthres = Plasticity(MinViscosityThreshold(nothing),
                                    CoulombYieldStress(;μ=pn[:μ],C)
                                    )
plas_const_minthres = Plasticity(MinViscosityThreshold(nothing),
                                    ConstantYieldStress(τy_const)
                                    )
damp = PrincipalKICharlesLaw(r)
dami = InvariantsKICharlesLaw(r)
sd = SimpleCharlesLaw()

lw = LinearWeakening(pn[:γ])

elas = IncompressibleElasticity(pn[:G])

cmp = ConstitutiveModel(;weakening=lw, damage=damp, elasticity=elas, plasticity=plas_coulomb_minthres)
cmi = ConstitutiveModel(;weakening=lw, damage=dami, elasticity=elas, plasticity=plas_coulomb_minthres)



csr = ConstantStrainRate(ϵ̇)
cs =ConstantStress(s)

loading2D = TriaxialSetup(; geom=Geom2D(), control=csr, pc) 
loading3D = TriaxialSetup(; geom=Geom3D(), control=csr, pc) 
creep2D = TriaxialSetup(; geom=Geom2D(), control=cs, pc) 
creep3D = TriaxialSetup(; geom=Geom3D(), control=cs, pc) 


mpl2D = Model(;setup=loading2D, cm=cmp)
mil2D = Model(;setup=loading2D, cm=cmi)
mil3D = Model(;setup=loading3D, cm=cmi)



@testset "helper functions" begin
    @test SD.confining_pressure(loading2D) == SD.confining_pressure(mpl2D) == pc

    @test SD.plastic_yield_stress(plas_coulomb_minthres,0) == plas_coulomb_minthres.σy.C * cos(atan(plas_coulomb_minthres.σy.μ))
    @test SD.plastic_yield_stress(plas_const_minthres,0) == plas_const_minthres.σy.val 

    Di_less = r.D₀ - 0.1
    Di_more = r.D₀ + 0.1
    @test init_variables(damp,Di_less)[2] == damp.r.D₀
    @test init_variables(damp,Di_more)[2] == Di_more
end

@testset "stress invariants" begin
    @test SD.stress_invariants(0, loading2D) == (0, -loading2D.pc)
    @test SD.stress_invariants(0, loading3D) == (0, -loading3D.pc)
    @test SD.stress_invariants(-1, loading2D) == (1, -1 - loading2D.pc)
    @test SD.stress_invariants(-1, loading3D) == (sqrt(3)/2, -0.5 - loading3D.pc)
end

@testset "strain invariants" begin
    @test SD.ϵII_func(1, mil2D)  == 1
    @test SD.ϵII_func(1, mil3D) == sqrt(3)/2
end

@testset "weakening functions and derivatives" begin
    # simple damage
    let cm = ConstitutiveModel(;weakening=lw, damage=sd, elasticity=elas, plasticity=plas_coulomb_minthres)
        γ = cm.weakening.γ
        @test SD.weakening_func(0,cm) == 1
        @test SD.weakening_func(1,cm) == γ
        @test SD.weakening_derivative(0,cm) == SD.weakening_derivative(1,cm) == γ - 1
    end

    # micromechanical damage
    let cm = ConstitutiveModel(;weakening=lw, damage=damp, elasticity=elas, plasticity=plas_coulomb_minthres)
        γ = cm.weakening.γ
        D₀ = cm.damage.r.D₀
        @test SD.weakening_func(D₀,cm) == 1
        @test SD.weakening_func(1,cm) == γ
        @test SD.weakening_derivative(0,cm) == SD.weakening_derivative(1,cm) == (γ-1)/(1-D₀)
    end

end

@testset "axial stress proxies conversions" begin
    
    let geom = Geom2D()
        # pc = 0
        let sₐₓ = -1, pc = 0
            setup = TriaxialSetup(; geom, control=csr, pc) 
            σₐₓtrue = -2
            Δσtrue = σₐₓtrue

            @test SD.sₐₓ2σₐₓ(setup,sₐₓ) == σₐₓtrue
            @test SD.σₐₓ2sₐₓ(setup,σₐₓtrue) == -1
            @test SD.sₐₓ2Δσ(setup,sₐₓ) == Δσtrue
            @test SD.Δσ2sₐₓ(setup,Δσtrue) == sₐₓ
        end
        # pc = 2
        let sₐₓ = -1, pc = 2
            setup = TriaxialSetup(; geom, control=csr, pc) 
            σₐₓtrue = -4
            Δσtrue = -2

            @test SD.sₐₓ2σₐₓ(setup,sₐₓ) == σₐₓtrue
            @test SD.σₐₓ2sₐₓ(setup,σₐₓtrue) == sₐₓ
            @test SD.sₐₓ2Δσ(setup,sₐₓ) == Δσtrue
            @test SD.Δσ2sₐₓ(setup,Δσtrue) == sₐₓ
        end
    end

    let geom = Geom3D()
        # pc = 0
        let sₐₓ = -1, pc = 0
            setup = TriaxialSetup(; geom, control=csr, pc)
            σₐₓtrue = -3/2
            Δσtrue = σₐₓtrue
            @test SD.sₐₓ2σₐₓ(setup,sₐₓ) == σₐₓtrue
            @test SD.σₐₓ2sₐₓ(setup,σₐₓtrue) == -1
            @test SD.sₐₓ2Δσ(setup,sₐₓ) == Δσtrue
            @test SD.Δσ2sₐₓ(setup,Δσtrue) == sₐₓ
        end
        # pc = 1
        let sₐₓ = -1, pc = 2
            setup = TriaxialSetup(; geom, control=csr, pc) 
            σₐₓtrue = -3/2 -2
            Δσtrue = -3/2
            @test SD.sₐₓ2σₐₓ(setup,sₐₓ) == σₐₓtrue
            @test SD.σₐₓ2sₐₓ(setup,σₐₓtrue) == sₐₓ
            @test SD.sₐₓ2Δσ(setup,sₐₓ) == Δσtrue
            @test SD.Δσ2sₐₓ(setup,Δσtrue) == sₐₓ
        end
    end

end

@testset "KI" begin
    # TEST KI WITH INVARIANTS
    let m = Model(;setup=loading2D, cm=cmi)
        a = m.cm.damage.r.a
        μ = m.cm.damage.r.μ
        D₀ = m.cm.damage.r.D₀
        pc = m.setup.pc
        ψ = m.cm.damage.r.ψ

        let s = -μ*pc/(1-μ)#- μ*pc/(0.5*(sqrt(3) - μ))
            τ, σₘ = SD.stress_invariants(s, m)
            # at D = D0 and τ/σ=μ -> KI = O
            @test DSB.compute_KI(r,σₘ,τ,D₀) ≈ 0
            @test SD.get_KI_from_s(s, D₀, m) ≈ 0
        end
    end

    let m = Model(;setup=loading3D, cm=cmi)
        a = m.cm.damage.r.a
        μ = m.cm.damage.r.μ
        D₀ = m.cm.damage.r.D₀
        pc = m.setup.pc
        ψ = m.cm.damage.r.ψ

        let s = - μ*pc/(0.5*(sqrt(3) - μ))
            τ, σₘ = SD.stress_invariants(s, m)
            # at D = D0 and τ/σ=μ -> KI = O
            @test isapprox(DSB.compute_KI(r,σₘ,τ,D₀), 0, atol=1e-10)
            @test isapprox(SD.get_KI_from_s(s, D₀, m), 0, atol=1e-10)
        end
    end

    # TEST KI WITH PRINCIPAL STRESSES
    let m = Model(;setup=loading2D, cm=cmp)
        a = m.cm.damage.r.a
        μ = m.cm.damage.r.μ
        D₀ = m.cm.damage.r.D₀
        pc = m.setup.pc
        ψ = m.cm.damage.r.ψ
        S = (sqrt(1+μ^2) + μ) / (sqrt(1+μ^2) - μ)

        # TEST KI INDEPENDENCE ON SIGN CONVENTION

        # σ₁ & σ₃ negative
        let σ₁ = -S*pc
            σ₃ = -pc
            σs = SD.flaw_shear_stress(σ₁,σ₃,ψ)
            σn = SD.flaw_normal_stress(σ₁,σ₃,ψ)
            @info "friction coef cancelling KI ($(abs(σs/σn))) is not equal to internal friction ($μ)"

            # at D = D0 and σ₁=Sσ₃ -> KI = O
            @test SD.get_KI_from_s(SD.σₐₓ2sₐₓ(m,σ₁), D₀, m) == 0
            # at D > D0 and σ₁=Sσ₃ -> KI = - sign(σ3)*c3*σ₃
            @test SD.get_KI_from_s(SD.σₐₓ2sₐₓ(m,σ₁), D₀ + 0.5(1-D₀), m)  ≈  -sign(σ₁) * sqrt(π*a) * σ₃ * SD.compute_c3(r,D₀ + 0.5(1-D₀))
            @test SD.get_KI_from_s(SD.σₐₓ2sₐₓ(m,σ₁), D₀ + 0.8(1-D₀), m)  ≈  -sign(σ₁) * sqrt(π*a) * σ₃ * SD.compute_c3(r,D₀ + 0.8(1-D₀))
        end
    end

end

# m = Model(;setup=loading2D, cm=cmp)
#         a = m.cm.damage.r.a
#         μ = m.cm.damage.r.μ
#         D₀ = m.cm.damage.r.D₀
#         pc = m.setup.pc
#         ψ = m.cm.damage.r.ψ
# S_func(μ) = (sqrt(1+μ^2) + μ) / (sqrt(1+μ^2) - μ)
#         σ₁ = -S*pc
#         σ₃ = -pc
#         σs = SD.flaw_shear_stress(σ₁,σ₃,ψ)
#         σn = SD.flaw_normal_stress(σ₁,σ₃,ψ)
        
# μ_vec = 0.3:0.01:0.8
# σ1_vec = [-S_func(μi)*pc for μi in μ_vec]

# σs_vec = [SD.flaw_shear_stress(σ₁i,σ₃,ψ) for σ₁i in σ1_vec]
# σn_vec = [SD.flaw_normal_stress(σ₁i,σ₃,ψ) for σ₁i in σ1_vec]

# lw = 4
# f,ax,pl = lines(μ_vec,abs.(σs_vec./σn_vec),linewidth=lw)
# lines!(ax,μ_vec[1]..μ_vec[end],x->x,linestyle=:dash,linewidth=lw)
# ax.xlabel = "μ"
# ax.ylabel = "τ/σₙ at KI=0"
# ax.xlabelsize = ax.ylabelsize = 30
# DataInspector();