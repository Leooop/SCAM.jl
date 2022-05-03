"blue notebook p123 :
This is the derivation of KI with mohr-coulomb friction, allowing free parametrization of ψ
!!! wrong derivation !!!"
# function compute_N1(r::Rheology,D)
#     α = cosd(r.ψ)
#     1/(((D/r.D₀)^(1/3)-1)*α + r.β)^(3/2)
# end

# function compute_N2(r::Rheology,D)
#     α = cosd(r.ψ)
#     2*sqrt(α)*((D-r.D₀)/r.D₀)^(1/6)
# end

# function compute_N3(r::Rheology,D)
#     α = cosd(r.ψ)
#     r.D₀^(2/3)/(α^2*(1-D^(2/3)))
# end

# compute_N1N2N3(r,D) = (compute_N1(r,D), compute_N2(r,D), compute_N3(r,D))
    
# compute_M1(r) = 0.5*sind(r.ψ)*( r.μ*(1-cosd(2*r.ψ)) - sind(2*r.ψ) )
# compute_M3(r) = 0.5*sind(r.ψ)*( r.μ*(1+cosd(2*r.ψ)) + sind(2*r.ψ) )
# compute_M1M3(r) = (compute_M1(r), compute_M3(r))
# # eq 15 Bhat2012 (A1 : *c2*c3), Perol&Bhat2016 (A1 : ...*c2)*c3):
# # Perol&Bhat2016 is the corrected version, and the one implemented
# compute_A1(N1,N2,N3,M1) = (N1 + N2*N3)*M1
# compute_A3(N1,N2,N3,M3) = (N1 + N2*N3)*M3 + N2

# compute_A1A3(N1,N2,N3,M1,M3) = (compute_A1(N1,N2,N3,M1), compute_A3(N1,N2,N3,M3))
# compute_A1A3(r::Rheology,D) = compute_A1A3(compute_N1N2N3(r,D)...,compute_M1M3(r)...)

# compute_KI(r, σ₁, σ₃, A1, A3) = sqrt(r.a/π) * (A1*σ₁ + A3*σ₃)
# compute_KI(r, σ₁, σ₃, D) = compute_KI(r, σ₁, σ₃, compute_A1A3(r,D)...)


# Bhat 2011 with KI twisted to match extreme cases where l=0 and l>>1 :
function compute_c1(r,D)
    α = cosd(r.ψ)
    return 1/(π^2 * α^(3/2) * ((D/r.D₀)^(1/3) - 1 + r.β/α)^(3/2))
end

function compute_c2(r,D)
    α = cosd(r.ψ)
    return 2 * ((D/r.D₀)^(1/3) - 1)^(1/2) * (r.D₀^(2/3)/(1-D^(2/3))) / (π^2 * α^(3/2)) 
end

function compute_c3(r,D)
    α = cosd(r.ψ)
    return 2/π * sqrt(α) * ((D/r.D₀)^(1/3) - 1)^(1/2) 
end
compute_c1c2c3(r::Rheology,D) = (compute_c1(r,D), compute_c2(r,D), compute_c3(r,D))


compute_A1(r) = π * sqrt(r.β/3) * (sqrt(1+r.μ^2) - r.μ) #TODO: Check if + or - μ makes sense
compute_A3(r, A1) = A1 * (sqrt(1+r.μ^2) + r.μ) / (sqrt(1+r.μ^2) - r.μ)
compute_A3(r) = compute_A1(r) * (sqrt(1+r.μ^2) + r.μ) / (sqrt(1+r.μ^2) - r.μ)
compute_A1A3(r::Rheology) = (A1 = compute_A1(r) ; A3 = compute_A3(r, A1) ; return (A1, A3))

"principal stress formulation of KI from Bhat 2011"
compute_KI(r::Rheology, σ₁, σ₃, A1, A3, c1, c2, c3) = sqrt(π*r.a) * ((σ₃*A3 - σ₁*A1) * (c1 + c2) + σ₃*c3)
compute_KI(r::Rheology, σ₁, σ₃, D) = compute_KI(r, σ₁, σ₃, compute_A1A3(r)..., compute_c1c2c3(r,D)...)

compute_KI_no_interaction(r::Rheology, σ₁, σ₃, A1, A3, c1, c2, c3) = sqrt(π*r.a) * ((σ₃*A3 - σ₁*A1)*c1 + σ₃*c3)
compute_KI_no_interaction(r::Rheology, σ₁, σ₃, D) = compute_KI_no_interaction(r, σ₁, σ₃, compute_A1A3(r)..., compute_c1c2c3(r,D)...)
#################
# use those instead to drop the match with extreme cases, vizualizing the difference shows a massive underestime of KI,
# it should not be used.

# c1 when droping Beta in the wedging force term of KI
function compute_c1psi(r,D)
    α = cosd(r.ψ)
    return 1/(π^2 * α^(3/2) * (D/r.D₀)^(1/2))
end
compute_c1c2c3psi(r::Rheology,D) = (compute_c1psi(r,D), compute_c2(r,D), compute_c3(r,D))
# mohr coulomb formulation of A1 & A3 :
compute_A1psi(r) = 0.5*π * sind(r.ψ) * (r.μ*(cosd(2*r.ψ) - 1) + sind(2*r.ψ))
compute_A3psi(r) = 0.5*π * sind(r.ψ) * (r.μ*(cosd(2*r.ψ) + 1) + sind(2*r.ψ))
compute_A1A3psi(r::Rheology) = (compute_A1psi(r), compute_A3psi(r))
##################
