var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = DiffusionFlux","category":"page"},{"location":"#DiffusionFlux","page":"Home","title":"DiffusionFlux","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for DiffusionFlux.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [DiffusionFlux]","category":"page"},{"location":"#DiffusionFlux.Properties","page":"Home","title":"DiffusionFlux.Properties","text":"A structure to define the porous medis properties \n\nϵ : porosity\nτ : tortuosity\npore_dia : pore diameter \npart_dia : particle diameter\n\n\n\n\n\n","category":"type"},{"location":"#DiffusionFlux.D_Kn-Tuple{Any, Properties, Any}","page":"Home","title":"DiffusionFlux.D_Kn","text":"Function to calcuate the Knudsen diffusion coeffcients \n\nUsage\n\nD_Kn(thermo_obj, pm::Properties, T)\n\nthermo_obj : SpeciesThermProperties (Refer IdealGas)    \npm : Struct of the type Properties    \nT : Temperature (K)\n\n\n\n\n\n","category":"method"},{"location":"#DiffusionFlux.KozeneyCarman-Tuple{Properties}","page":"Home","title":"DiffusionFlux.KozeneyCarman","text":"Function to calculate the permeability using Kozeney Carman relationship\n\nUsage\n\nKozeneyCarman(pm::Properties)\n\npm: Struct of the type Properties    \n\n\n\n\n\n","category":"method"},{"location":"#DiffusionFlux.dgm_objects!-Tuple{DGMObjects, Properties, Any, Float64, Float64, Array{Float64}}","page":"Home","title":"DiffusionFlux.dgm_objects!","text":"Function to calculate binary diffusion coeffcients and Knudsen diffusion coeffcients, which are required for the evaluation of DGM fluxes\n\nUsage\n\ndgmobjects!(dgmobj::DGMObjects, pm::Properties, sp_trd, p::Float64, T::Float64 ,molwts::Array{Float64})\n\ndgm_obj : A struct of the type DGMObjects which stores the diffusion coeffcients matrix\npm : Porous media properties (struct of the type Properties)\nsptrd : Array of Species transport data, which is obtained by calling createtransport_data of TransportProperties\np : pressure in (Pa)\nT : Temperature (K)\nmolwts : molecular weights vector \n\n\n\n\n\n","category":"method"},{"location":"#DiffusionFlux.flux_dgm!-Union{Tuple{T}, Tuple{Array{Vector{T}, 1}, Array{Vector{T}, 1}, Properties, DGMObjects, Any, Array{T}, T, T}} where T","page":"Home","title":"DiffusionFlux.flux_dgm!","text":"Function to calculate the mass fluxes (kg/m^2-s) at the interface between two cells using DGM  The fluxes at the boundary cells are not evaluated in this function. If There are n cells, there will be n-1 internal cell faces \n\nUsage\n\nfluxdgm!(jks::Array{Array{T,1},1}, C::Array{Array{T,1},1}, pm::Properties, dgmobj::DGMObjects , molwts::Array{T}, Temp::T, δ::T)\n\njks : vector{vector} for storing the fluxes jks[n] stands for the flux at the interface between n and (n+1)th cell\nC : vector{vector} concentrations in all cells \npm : Porous media properties (struct of the type Properties)\ndgm_obj : A struct of the type DGMObjects which stores the diffusion coeffcients matrix\nmolwts : vector of molecular weights\n\n\n\n\n\n","category":"method"},{"location":"#DiffusionFlux.flux_ficks!-Union{Tuple{T}, Tuple{Array{T}, Array{T}, Array{T}, Array{T}, Any, T}} where T","page":"Home","title":"DiffusionFlux.flux_ficks!","text":"Function to calculate the diffusion fluxes in Kg/m^2-s using Ficks law. The fluxes are corrected to  ensure that the sum is zero \n\nUsage\n\nflux_ficks!(jks::Array{T}, Dkm::Array{T} ,C1::Array{T}, C2::Array{T}, δ::T)\n\njks : vector of fluxes \nC1 : concentration vector\nC2 : concentration vector\nδ  : distance between the cell centers  \n\n\n\n\n\n","category":"method"},{"location":"#DiffusionFlux.interface_flux!-Union{Tuple{T}, Tuple{Vector{T}, Any, Any, DGMObjects, Array{T}, T}} where T","page":"Home","title":"DiffusionFlux.interface_flux!","text":"Function to calcuate the mass flux (kg/m^2-s) at the inteface between flow channel and porous media. The pressure driven flux is assumed to be negligible at the inteface \n\nUsage\n\ninterfaceflux!(jks::Array{T,1},Cch, Cpm, dgmobj::DGMObjects, molwts::Array{T}, Temp::T, δ::T)\n\njks : Storage for species flux (1D Array)\nC_ch : concentration in the flow channel\nC_pm : concentration in the porous media\ndgm_obj : Diffusion matrix for DGM flux calculations\nmolwts  : vector of molecular weights\nTemp : Temperature (K)\nδ   : distance between the cell centers \n\n\n\n\n\n","category":"method"}]
}
