using DiffusionFlux, IdealGas, TransportProperties
using Test

@testset "DiffusionFlux.jl" begin

    if Sys.isapple() || Sys.islinux()
        lib_dir = "lib/"
    elseif Sys.iswindows()
        lib_dir = "lib\\"
    end
    tr_file = joinpath(lib_dir, "transport.dat")
    therm_file = joinpath(lib_dir, "therm.dat")
    species = ["H2", "H2O", "CO", "CH4", "CO2", "O2", "N2"]
    sp_trd = create_transport_data(species, tr_file)
    thermo_obj = create_thermo(species, therm_file)
    # channel concentration
    ml1 = [0.12, 0.66, 0.0, 0.0, 0.22, 0.0, 0.0]
    # concentration within the electrode 
    ml2 = [0.1, 0.5, 0.4, 0.0, 0.0, 0.0, 0.0]
    T = 1073.15
    Cb = 1e5/IdealGas.R/T
    c1 = Cb * ml1 
    c2 = Cb * ml2

    @testset "Testing Fickinan diffusion " begin
        Dkm = similar(c1)
        mole_fracs = 0.5*(ml1+ml2)
        D_km!(Dkm, sp_trd, T, 1e5, thermo_obj.molwt, mole_fracs)
        δ=1e-4
        jks = similar(c1)
        flux_ficks!(jks, Dkm, c1, c2, thermo_obj.molwt, δ)
        @test sum(jks) == 0.0
    end

    @testset "Testing DGM Fluxes (interface) " begin
        pm = Properties(0.35, 3.5, 1e-6, 2.5e-6)
        n = length(c1)
        Dkn_e = Array{Float64,1}(undef, n)
        D_ij_e = Matrix{Float64}(undef, n, n)
        Dkl_DGM = Matrix{Float64}(undef, n, n)
        Dkm = Array{Float64,1}()
        dgm_objs = WorkSpace(Dkn_e, D_ij_e, Dkl_DGM, Dkm)
        δ=4.88e-5
        jks = similar(c1)
        effective_coefficients!(dgm_objs, pm, sp_trd, 1e5, T, thermo_obj.molwt)
        flux_interface!(jks,c1, c2, dgm_objs, thermo_obj.molwt, δ)        
        @test sum(jks) < 0.005
    end


    @testset "Testing DGM fluxes (porous media) " begin
        pm = Properties(0.35, 3.5, 1e-6, 2.5e-6)
        n = length(c1)
        Dkn_e = Array{Float64,1}(undef, n)
        D_ij_e = Matrix{Float64}(undef, n, n)        
        Dkl_DGM = Matrix{Float64}(undef, n, n)
        Dkm = Array{Float64,1}()
        dgm_objs = WorkSpace(Dkn_e, D_ij_e, Dkl_DGM, Dkm)
        δ=4.88e-5        
        effective_coefficients!(dgm_objs, pm, sp_trd, 1e5, T, thermo_obj.molwt)
        
        ncells = 10
        jks = Array{Array{Float64,1},1}(undef, ncells-1)
        conc = Array{Array{Float64,1},1}(undef, ncells)
        for i in 1:ncells
            conc[i] = c1
        end
        
        flux_dgm!(jks, conc, pm, dgm_objs, sp_trd,thermo_obj.molwt, T, δ)
        total = reduce(.+, reduce(.+, jks))
        @test total == 0

    end

    @testset "Testing modified Ficks model" begin
        pm = Properties(0.35, 3.5, 1e-6, 2.5e-6)
        n = length(c1)
        Dkn_e = Array{Float64,1}(undef, n)
        D_ij = Matrix{Float64}(undef, n, n)        
        Dkl_DGM = Matrix{Float64}(undef,1,1)
        Dkm = Array{Float64,1}(undef, n)
        diff_coeffs = WorkSpace(Dkn_e, D_ij, Dkl_DGM, Dkm)
        ncells = 10
        conc = Array{Array{Float64,1},1}(undef, ncells)        
        for i in 1:ncells
            conc[i] = c1
        end
        δ=4.88e-5        
        jks = Array{Array{Float64,1},1}(undef, ncells-1)
        flux_porous_media_fick!(jks, conc, pm, sp_trd, diff_coeffs , thermo_obj.molwt, T, δ)
        total = reduce(.+, reduce(.+, jks))
        @test total == 0
    end

end
