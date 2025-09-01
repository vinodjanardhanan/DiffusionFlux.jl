module DiffusionFlux

using IdealGas, TransportProperties, RxnHelperUtils, LinearAlgebra

export  D_Kn, KozeneyCarman, effective_coefficients!, flux_ficks!, flux_interface!, flux_dgm!, flux_porous_media_fick!
export flux_dgm_Dij_update!

"""
A structure to define the porous medis properties 
-   ϵ : porosity
-   τ : tortuosity
-   pore_dia : pore diameter 
-   part_dia : particle diameter
"""
struct Properties
    ϵ::Float64
    τ::Float64
    pore_dia::Float64
    part_dia::Float64
end
export Properties


"""
Non allocating array for calculation of fluxes 
-   D_kn_e : Array for Knudsen diffusion coeffcients
-   D_ij_e : Array for binary diffusion coeffcients
-   D_kl_DGM : Array space for use in DGM 
-   Dkm : Array for storing diffusion coeffcients of k in the rest of the mixture
"""
mutable struct WorkSpace
    D_kn_e::Array{Float64}
    D_ij_e::Matrix{Float64}    
    Dkl_DGM::Matrix{Float64}
    Dkm::Array{Float64}
    I_matrix::Matrix{Float64} # Identity matrix
    p_vec::Array{Float64} # pressure vector    
end

export WorkSpace


"""
Function to calculate the permeability using Kozeney Carman relationship
# Usage
    KozeneyCarman(pm::Properties)
-   pm: Struct of the type Properties    
"""
KozeneyCarman(pm::Properties) = pm.ϵ^3 * pm.part_dia^2/(72pm.τ *(1-pm.ϵ)^2)
    

"""
Function to calcuate the Knudsen diffusion coeffcients 
# Usage
    D_Kn(thermo_obj, pm::Properties, T)
-   thermo_obj : SpeciesThermProperties (Refer IdealGas)    
-   pm : Struct of the type Properties    
-   T : Temperature (K)
"""
D_Kn(molwts, pm::Properties, T) = (pm.pore_dia/3.0)*sqrt(8*IdealGas.R*T/π) * sqrt.(1.0 ./ molwts)
    

#=
Function to evaluate H matrix which is part of Dusty Gas Model
=#
function H_matrix!(H::Matrix, D_kn_e::Array{Float64,1}, D_ij_e::Matrix{Float64}, mole_fracs::Array{Float64,1})
    for k in eachindex(mole_fracs)
        sum = 0.0
        for l in eachindex(mole_fracs)
            H[k,l] = -mole_fracs[k]/D_ij_e[k,l]
            if k != l
                sum += mole_fracs[l]/D_ij_e[k,l]
            end
        end
        H[k,k] = (1/D_kn_e[k]) + sum       
    end            
end

#=
Function to evaluate Matrix of Diffusion coeffcients for Dusty Gas Model 
=#
function D_kl_DGM!(Dkl_DGM, D_kn_e::Array{Float64,1}, D_ij_e::Matrix{Float64}, mole_fracs::Array{Float64}, I_mat::Matrix{Float64})
    # n = length(mole_fracs)
    H_matrix!(Dkl_DGM, D_kn_e, D_ij_e, mole_fracs)   
    F = lu!(Dkl_DGM)     
    ldiv!(Dkl_DGM, F, I_mat)
end


"""
Function to calculate binary diffusion coeffcients and Knudsen diffusion coeffcients,
which are required for the evaluation of DGM fluxes
# Usage  
effective_coefficients!(dgm_obj::WorkSpace, pm::Properties, sp_trd, p::Float64, T::Float64 ,molwts::Array{Float64})
-   dgm_obj : A struct of the type WorkSpace which stores the diffusion coeffcients matrix
-   pm : Porous media properties (struct of the type Properties)
-   sp_trd : Array of Species transport data, which is obtained by calling create_transport_data of TransportProperties
-   p : pressure in (Pa)
-   T : Temperature (K)
-   molwts : molecular weights vector 
"""
function effective_coefficients!(dgm_obj::WorkSpace, pm::Properties, sp_trd, p::Float64, T::Float64, molwts::Array{Float64})
    # calculate the effective binary diffusion coeffcients    
    dgm_obj.D_ij_e = (pm.ϵ/pm.τ) * D_ij(sp_trd, T, p, molwts)
    
    # calculate the effective Knudsen diffusion coeffcients
    dgm_obj.D_kn_e = (pm.ϵ/pm.τ) * D_Kn(molwts, pm, T)

end

"""
Function to calculate the mass fluxes (kg/m^2-s) at the interface between two cells using DGM 
The fluxes at the boundary cells are not evaluated in this function. If There
are n cells, there will be n-1 internal cell faces 
# Usage 
flux_dgm!(jks::Array{Array{T,1},1}, C::Array{Array{T,1},1}, pm::Properties, dgm_obj::WorkSpace , molwts::Array{T}, Temp::T, δ::T)
-   jks : vector{vector} for storing the fluxes jks[n] stands for the flux at the interface between n and (n+1)th cell
-   C : vector{vector} concentrations in all cells 
-   pm : Porous media properties (struct of the type Properties)
-   dgm_obj : A struct of the type WorkSpace which stores the diffusion coeffcients matrix
-   sp_tr_data : Array of Species transport data, which is obtained by calling create_transport_data of TransportProperties
-   molwts : vector of molecular weights
-   Temp : Temperature (K)
-   δ : distance between cell centers 
"""
function flux_dgm!(jks::Array{Array{T,1},1}, C::Matrix{T}, pm::Properties, dgm_obj::WorkSpace, sp_tr_data, molwts::Array{T}, Temp::T, δ::T) where T
    # Permeability
    Bg = KozeneyCarman(pm)
    # number of cells 
    ncells = size(C, 2)
    # pressure in each cell 
    for i in 1:ncells
        dgm_obj.p_vec[i] = sum(@view C[:, i]) * IdealGas.R * Temp
    end
    
    #=
    Only an approximate value of viscosity is enough for these calculations 
    =#
    # mole_fracs = sum(map(Cs->Cs/sum(Cs) , C))/ncells # this is an approximate value averaged over all cells
    mole_fracs = vec(sum(C, dims=2)) ./ sum(C) # approximate molefractions 
    μ = viscosity(sp_tr_data,Temp,molwts,mole_fracs)    
    sp_flux = zeros(length(mole_fracs))
    for i in 1:ncells-1
        C_1  = @view C[:, i]
        C_2  = @view C[:, i+1]
        ∇C = (C_2- C_1 )/δ  # vector of concentration gradients 
        ∇p = (dgm_obj.p_vec[i+1]-dgm_obj.p_vec[i])/δ # pressure gradient in the cell 
        C_avg = 0.5*(C_1 + C_2 ) # Average concentration at the cell faces 
        mole_fracs = C_avg/sum(C_avg) # average mole fractions at the cell faces 
        D_kl_DGM!(dgm_obj.Dkl_DGM, dgm_obj.D_kn_e, dgm_obj.D_ij_e, mole_fracs, dgm_obj.I_matrix)
        
        @inbounds @simd for k in eachindex(mole_fracs)
            jk_1 = sum(dgm_obj.Dkl_DGM[k,:] .* ∇C)
            jk_2 = sum((dgm_obj.Dkl_DGM[k,:] .* C_avg) ./ dgm_obj.D_kn_e)*Bg*∇p/μ
            sp_flux[k] = jk_1 + jk_2
        end
        jks[i] = -sp_flux .* molwts
    end
                
end



"""
Function to calculate the mass fluxes (kg/m^2-s) at the interface between two cells using DGM 
The fluxes at the boundary cells are not evaluated in this function. If There
are n cells, there will be n-1 internal cell faces 
# Usage 
flux_dgm_Dij_update!(jks::Array{Array{T,1},1}, C::Array{Array{T,1},1}, pm::Properties, dgm_obj::WorkSpace , molwts::Array{T}, Temp::T, δ::T)
-   jks : vector{vector} for storing the fluxes jks[n] stands for the flux at the interface between n and (n+1)th cell
-   C : vector{vector} concentrations in all cells 
-   pm : Porous media properties (struct of the type Properties)
-   dgm_obj : A struct of the type WorkSpace which stores the diffusion coeffcients matrix
-   sp_tr_data : Array of Species transport data, which is obtained by calling create_transport_data of TransportProperties
-   molwts : vector of molecular weights
-   Temp : Temperature (K)
-   δ : distance between cell centers 
"""
# function flux_dgm_Dij_update!(jks::Array{Array{T,1},1}, C::Array{Array{T,1},1}, pm::Properties, dgm_obj::WorkSpace, sp_tr_data, molwts::Array{T}, Temp::T, δ::T) where T
#     # Permeability
#     Bg = KozeneyCarman(pm)
#     # number of cells 
#     ncells = length(C)
#     # pressure in each cell 
#     p_vec = map(Cs->sum(Cs)*IdealGas.R*Temp, C)
    
#     #=
#     Only an approximate value of viscosity is enough for these calculations 
#     =#
#     mole_fracs = sum(map(Cs->Cs/sum(Cs) , C))/ncells
#     μ = viscosity(sp_tr_data,Temp,molwts,mole_fracs)    

#     for i in 1:ncells-1
#         ∇C = (C[i+1]-C[i])/δ  # vector of concentration gradients 
#         dgm_obj.D_ij_e = (pm.ϵ/pm.τ) * D_ij(sp_tr_data, Temp, p_vec[i], molwts)        
#         ∇p = (p_vec[i+1]-p_vec[i])/δ # pressure gradient in the cell 
#         C_avg = 0.5*(C[i+1]+C[i]) # Average concentration at the cell faces 
#         mole_fracs = C_avg/sum(C_avg) # average mole fractions at the cell faces 
#         D_kl_DGM!(dgm_obj.Dkl_DGM, dgm_obj.D_kn_e, dgm_obj.D_ij_e, mole_fracs)
#         sp_flux = zeros(length(mole_fracs))
#         for k in eachindex(mole_fracs)
#             jk_1 = sum(dgm_obj.Dkl_DGM[k,:] .* ∇C)
#             jk_2 = sum((dgm_obj.Dkl_DGM[k,:] .* C_avg) ./ dgm_obj.D_kn_e)*Bg*∇p/μ
#             sp_flux[k] = jk_1 + jk_2
#         end
#         jks[i] = -sp_flux .* molwts
#     end
                
# end

"""
Function to calcuate the mass flux (kg/m^2-s) at the inteface between flow channel and porous media.
The pressure driven flux is assumed to be negligible at the inteface 
# Usage
interface_flux!(jks::Array{T,1},C_ch, C_pm, dgm_obj::WorkSpace, molwts::Array{T}, Temp::T, δ::T)
-   jks : Storage for species flux (1D Array)
-   C_ch : concentration in the flow channel
-   C_pm : concentration in the porous media
-   dgm_obj : Diffusion matrix for DGM flux calculations
-   molwts  : vector of molecular weights
-   Temp : Temperature (K)
-   δ   : distance between the cell centers 
"""
function flux_interface!(jks::Array{T,1},C_ch, C_pm, dgm_obj::WorkSpace, molwts::Array{T}, δ::T) where T
    ∇C = (C_ch - C_pm)/δ    
    C_avg = 0.5*(C_ch+C_pm)
    mole_fracs = C_avg/sum(C_avg) # average mole fractions at the cell faces 
    # n = length(C_ch)
    # Dkl_DGM = zeros(n,n)    
    D_kl_DGM!(dgm_obj.Dkl_DGM, dgm_obj.D_kn_e, dgm_obj.D_ij_e, mole_fracs, dgm_obj.I_matrix)
    for k in eachindex(C_ch)
        jks[k] = -sum(dgm_obj.Dkl_DGM[k,:] .* ∇C) * molwts[k]
    end

end


"""
Function to calculate the diffusion fluxes in Kg/m^2-s using Ficks law. The fluxes are corrected to 
ensure that the sum is zero 
# Usage
flux_ficks!(jks::Array{T}, Dkm::Array{T} ,C1::Array{T}, C2::Array{T}, δ::T)
-   jks : vector of fluxes 
-   C1 : concentration vector
-   C2 : concentration vector
-   δ  : distance between the cell centers  
"""
function flux_ficks!(jks::Array{T}, Dkm::Array{T} ,C1::Array{T}, C2::Array{T}, molwts, δ::T) where T
    @. jks = -Dkm * (C1-C2)/δ
    jks .*= molwts # convert to kg/m^2-s
    j_corr = sum(jks)
    C = 0.5*(C1+C2)
    mole_fracs = C/sum(C)
    mass_fracs = similar(mole_fracs)
    molefrac_to_massfrac!(mass_fracs, mole_fracs, molwts)
    @. jks -= mass_fracs*j_corr
end    



"""
Function to calculate the mass fluxes (kg/m^2-s) at the interface between two cells using modified Fickian model 
The fluxes at the boundary cells are not evaluated in this function. If There
are n cells, there will be n-1 internal cell faces 
# Usage 
flux_porous_media_fick!(jks::Array{Array{T,1},1}, C::Array{Array{T,1},1}, pm::Properties, sp_tr_data, ws::WorkSpace , molwts::Array{T}, Temp::T, δ::T)
-   jks : vector{vector} for storing the fluxes jks[n] stands for the flux at the interface between n and (n+1)th cell
-   C : vector{vector} concentrations in all cells 
-   pm : Porous media properties (struct of the type Properties)
-   sp_tr_data : Array of Species transport data, which is obtained by calling create_transport_data of TransportProperties
-   ws : A struct of the type WorkSpace which stores the diffusion coeffcients matrix
-   molwts : vector of molecular weights
-   Temp : Temperature (K)
-   δ : distance between cell centers 
"""
function flux_porous_media_fick!(jks::Array{Array{T,1},1}, C::Array{Array{T,1},1}, pm::Properties, sp_tr_data, ws::WorkSpace , molwts::Array{T}, Temp::T, δ::T) where T
    ncells = length(C)    
    p_vec = map(Cs->sum(Cs)*IdealGas.R*Temp, C)
    # average pressure 
    p = sum(p_vec)/ncells 
    ncell = length(C)
    # Permeability
    Bg = KozeneyCarman(pm)
    ws.D_kn_e = (pm.ϵ/pm.τ)*D_Kn(molwts, pm, Temp)
    # this is not effective binary diffusion coeffcients
    ws.D_ij_e = D_ij(sp_tr_data, Temp, p, molwts)
    
    # Average mole fractions in the whole porous media 
    mole_fracs = sum(map(Cs->Cs/sum(Cs) , C))/ncells
    # calculate the viscosity based on the average mole fractions 
    μ = viscosity(sp_tr_data,Temp,molwts,mole_fracs)        
    for i in 1:ncell-1
        # mole fractions at the inteface between i'th and i+1'th cell 
        mole_fracs = 0.5*( C[i]/sum(C[i]) + C[i+1]/sum(C[i+1]))
        ∇C = (C[i+1]-C[i])/δ
        D_km!(ws.Dkm, ws.D_ij_e, mole_fracs, molwts)
        ws.Dkm *= pm.ϵ/pm.τ
        De = (1.0 ./ ws.D_kn_e) + (1.0 ./ ws.Dkm)
        De = 1.0 ./ De
        jk1 = -De .* ∇C
        jk2 = - Bg*(p_vec[i+1]-p_vec[i])/δ/μ
        jks[i] = jk1 .+ jk2
        @. jks[i] *= molwts
    end
    
end


end
