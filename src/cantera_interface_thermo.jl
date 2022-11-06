#####################################################
#=
Interface to cantera/src/clib/ct.cpp file thermo objects
This contains all information about the thermodynamic state and composition
Notes: 
1. Cantera treats temp and density as independent variables
    it will hold there constant UNLESS they are explicitly changed
    or pressure is directly changed


C. Pavan 2022-10-28
=#
####################################################

######## Base object & constructors ###########

struct thermo_base
    ind::Int
    nspec::Int
end

# number of species
function _num_spec(ph::Integer)
    nspec=ccall(Libdl.dlsym(lib,:thermo_nSpecies),
        Csize_t,(Cint,),
        ph)
    return nspec
end

# create the thermo base
# if full path is true, must give entire path to the input file
# otherwise, searches in usr/local/share/cantera/data/
function thermo_base(file::String; full_path::Bool=true)
    # load in the file
    sym=Libdl.dlsym(lib,:thermo_newFromFile)
    if !full_path
        file="/usr/local/share/cantera/data/" * file
    end
    therm_obj=ccall(sym,
        Cint, (Cstring, Cstring) ,
        file,"")
    nspec = _num_spec(therm_obj)
    return thermo_base(therm_obj,nspec)
end

############## Single thermo state setting #################

# Temperature (using this over-rides pressure)
function set_T(tbase::thermo_base, T::Float64)
    sym=Libdl.dlsym(lib,:thermo_setTemperature)
    ccall(sym,
        Cint,(Cint,Cdouble),
        tbase.ind,T)
    return nothing
end

# pressure  (using this over-rides density)
function set_P(tbase::thermo_base, P::Float64)
    sym=Libdl.dlsym(lib,:thermo_setPressure)
    ccall(sym,
        Cint,(Cint,Cdouble),
        tbase.ind,P)
    return nothing
end

# density (using this over-rides pressure)
function set_rho(tbase::thermo_base, ρ::Float64)
    sym=Libdl.dlsym(lib,:thermo_setDensity)
    ccall(sym,
        Cint,(Cint,Cdouble),
        tbase.ind,ρ)
    return nothing
end

########### Multi-thermo state setting #############

# set enthalpy + pressure
function set_HP(tbase::thermo_base,HP::Tuple{Float64,Float64})
    sym=Libdl.dlsym(lib,:thermo_set_HP)
    HP_vec=[HP[1],HP[2]]
    ccall(sym,
        Cint,(Cint,Ptr{Cdouble}),
        tbase.ind,HP_vec)
    return nothing
end

# set int energy + pressure
function set_EP(tbase::thermo_base,EP::Tuple{Float64,Float64})
    sym=Libdl.dlsym(lib,:thermo_set_UP)
    EP_vec=[EP[1],EP[2]]
    ccall(sym,
        Cint,(Cint,Ptr{Cdouble}),
        tbase.ind,EP_vec)
    return nothing
end

# set int energy + density
function set_ER(tbase::thermo_base,ER::Tuple{Float64,Float64})
    sym=Libdl.dlsym(lib,:thermo_set_UV)
    EV_vec=[ER[1],1/ER[2]]
    ccall(sym,
        Cint,(Cint,Ptr{Cdouble}),
        tbase.ind,EV_vec)
    return nothing
end

# set temp + density
function set_TR(tbase::thermo_base,TR::Tuple{Float64,Float64})
    set_T(tbase,TR[1])
    set_rho(tbase,TR[2])
    return nothing
end

# set density + pressure
function set_RP(tbase::thermo_base,RP::Tuple{Float64,Float64})
    sym=Libdl.dlsym(lib,:thermo_set_RP)
    RP_vec=[RP[1],RP[2]]
    ccall(sym,
        Cint,(Cint,Ptr{Cdouble}),
        tbase.ind,RP_vec)
    return nothing
end

########### Composition Setting #############

# species mole fractions
# By name, X has format "A:1,B:2,C:3"
function set_X(tbase::thermo_base, X::String)
    sym=Libdl.dlsym(lib,:thermo_setMoleFractionsByName)
    ccall(sym,
        Cint,(Cint,Cstring),
        tbase.ind,X)
    return nothing
end
# set all species
function set_X(tbase::thermo_base,X::Array{Float64,1},norm=true)
    # norm is a boolean flag for whether to normalize
    sym=Libdl.dlsym(lib,:thermo_setMoleFractions)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble},Cint),
        tbase.ind,tbase.nspec,X, norm)
    return nothing
end

# species mass fractions
# By name, Y has format "A:1,B:2,C:3"
function set_Y(tbase::thermo_base, Y::String)
    sym=Libdl.dlsym(lib,:thermo_setMassFractionsByName)
    ccall(sym,
        Cint,(Cint,Cstring),
        tbase.ind,Y)
    return nothing
end
# set all species
function set_Y(tbase::thermo_base,Y::Array{Float64,1},norm=true)
    # norm is a boolean flag for whether to normalize
    sym=Libdl.dlsym(lib,:thermo_setMassFractions)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble},Cint),
        tbase.ind,tbase.nspec,Y, norm)
    return nothing
end

######### Set full state ###############

function set_TPX(thermo::thermo_base,TPX::Tuple{Float64,Float64,Union{String,Array{Float64,1}}})
    set_T(thermo,TPX[1])
    set_X(thermo,TPX[3])
    set_P(thermo,TPX[2])
    return nothing
end
function set_TPY(thermo::thermo_base,TPY::Tuple{Float64,Float64,Union{String,Array{Float64,1}}})
    set_T(thermo,TPY[1])
    set_Y(thermo,TPY[3])
    set_P(thermo,TPY[2])
    return nothing
end
function set_HPX(thermo::thermo_base,HPX::Tuple{Float64,Float64,Union{String,Array{Float64,1}}})
    set_Y(thermo,HPX[3])
    set_HP(thermo,HPX[1:2])
    return nothing
end
function set_HPY(thermo::thermo_base,HPY::Tuple{Float64,Float64,Union{String,Array{Float64,1}}})
    set_Y(thermo,HPY[3])
    set_HP(thermo,HPY[1:2])
    return nothing
end

# set by internal energy (E), density (R) and mass fraction (Y)
function set_ERY(thermo::thermo_base,ERY::Tuple{Float64,Float64,Union{String,Array{Float64,1}}})
    set_Y(thermo,ERY[3])
    set_ER(thermo,ERY[1:2])
    return nothing
end

# set by internal energy (E), density (R) and mole fraction (X)
function set_ERX(thermo::thermo_base,ERX::Tuple{Float64,Float64,Union{String,Array{Float64,1}}})
    set_X(thermo,ERX[3])
    set_ER(thermo,ERX[1:2])
    return nothing
end


############### Get Thermo State #################3
# temperature
function get_T(tbase::thermo_base)
    T=ccall(Libdl.dlsym(lib,:thermo_temperature),
        Float64,(Cint,),
        tbase.ind)
    return T
end

# pressure
function get_P(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_pressure)
    P=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return P
end

# density
function get_rho(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_density)
    rho=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return rho
end

# cp
function get_cp(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_cp_mass)
    cp=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return cp
end

# cv
function get_cv(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_cv_mass)
    cv=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return cv
end


# molar cp
function get_cp_mole(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_cp_mole)
    cp=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return cp
end

# enthalpy
function get_h(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_enthalpy_mass)
    h=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return h
end

# get E
function get_e(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_intEnergy_mass)
    e=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return e
end

# mean molecular weight
function get_mean_MW(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_meanMolecularWeight)
    MW=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return MW
end

# species molecular weight
function get_MW(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_getMolecularWeights)
    MW=Array{Float64,1}(undef, tbase.nspec)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble}),
        tbase.ind, tbase.nspec,MW)
    return MW
end

# species cp
function get_spec_molar_cp(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_getCp_R)
    cp=Array{Float64,1}(undef, tbase.nspec)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble}),
        tbase.ind, tbase.nspec,cp)
    cp*=Ru
    return cp
end

# species enthalpies
function get_spec_molar_enthalpies(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_getEnthalpies_RT)
    h=Array{Float64,1}(undef, tbase.nspec)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble}),
        tbase.ind, tbase.nspec,h)
    h*=Ru*get_T(tbase)
    return h
end 


############ Get Composition ##################

# get all species
function get_X(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_getMoleFractions)
    X=Array{Float64,1}(undef,tbase.nspec)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble}),
        tbase.ind,tbase.nspec,X)
    return X
end

# get all species
function get_Y(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_getMassFractions)
    Y=Array{Float64,1}(undef,tbase.nspec)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble}),
        tbase.ind,tbase.nspec,Y)
    return Y
end

############### Misc ##########################3

# get species names -> used by gas object
function _get_spec_name(tbase::thermo_base,spec_num::Integer)
    # remember cpp indexes from 0 
    spec_num-=1   
    sym=Libdl.dlsym(lib,:thermo_getSpeciesName)
    nm=""
    len=ccall(sym,
        Cint,(Cint,Csize_t,Csize_t,Cstring),
        tbase.ind,spec_num,1,nm)
    if len!=0
        nm=string(0,base=10,pad=len)
        ccall(sym,
        Cint,(Cint,Csize_t,Csize_t,Cstring),
        tbase.ind,spec_num,len,nm)
    end 
    return nm[1:len-1]
end
