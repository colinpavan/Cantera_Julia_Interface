# base thermo object - stores info about the phase
# does not have any functions for accessing the c object
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
function thermo_base(file::String; full_path::Bool=false)
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

#### All the functions for accessing/setting thermo properties

# Temperature access
function set_T(tbase::thermo_base, T::Float64)
    sym=Libdl.dlsym(lib,:thermo_setTemperature)
    ccall(sym,
        Cint,(Cint,Cdouble),
        tbase.ind,T)
    return nothing
end

function get_T(tbase::thermo_base)
    T=ccall(Libdl.dlsym(lib,:thermo_temperature),
        Float64,(Cint,),
        tbase.ind)
    return T
end

# pressure access
function set_P(tbase::thermo_base, P::Float64)
    sym=Libdl.dlsym(lib,:thermo_setPressure)
    ccall(sym,
        Cint,(Cint,Cdouble),
        tbase.ind,P)
    return nothing
end

function get_P(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_pressure)
    P=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return P
end

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

# get all species
function get_X(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_getMoleFractions)
    X=Array{Float64,1}(undef,tbase.nspec)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble}),
        tbase.ind,tbase.nspec,X)
    return X
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

# get all species
function get_Y(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_getMassFractions)
    Y=Array{Float64,1}(undef,tbase.nspec)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble}),
        tbase.ind,tbase.nspec,Y)
    return Y
end

# get thermo properties
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

# Enthalpy access functions
# get enthalpy
function get_h(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_enthalpy_mass)
    h=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return h
end

# set enthalpy (with pressure)
function set_HP(tbase::thermo_base,HP::Tuple{Float64,Float64})
    sym=Libdl.dlsym(lib,:thermo_set_HP)
    HP_vec=[HP[1],HP[2]]
    ccall(sym,
        Cint,(Cint,Ptr{Cdouble}),
        tbase.ind,HP_vec)
    return nothing
end

# Internal Energy access functions
# get E
function get_e(tbase::thermo_base)
    sym=Libdl.dlsym(lib,:thermo_intEnergy_mass)
    e=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return e
end

# set internal energy (with pressure)
function set_EP(tbase::thermo_base,EP::Tuple{Float64,Float64})
    sym=Libdl.dlsym(lib,:thermo_set_UP)
    EP_vec=[EP[1],EP[2]]
    ccall(sym,
        Cint,(Cint,Ptr{Cdouble}),
        tbase.ind,EP_vec)
    return nothing
end

# set internal energy (with density)
function set_ER(tbase::thermo_base,ER::Tuple{Float64,Float64})
    sym=Libdl.dlsym(lib,:thermo_set_UV)
    EV_vec=[ER[1],1/ER[2]]
    ccall(sym,
        Cint,(Cint,Ptr{Cdouble}),
        tbase.ind,EV_vec)
    return nothing
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


# multi-parameter setters

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

# get species names
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

#=
# mutable thermo phase object
mutable struct thermo
    base::thermo_base
    set_T::Function
    get_T::Function
    set_P::Function
    get_P::Function
    set_X::Function
    get_X::Function
    set_Y::Function
    get_Y::Function
    set_TPX::Function
    set_TPY::Function
    rho::Function
    cp::Function
    cp_molar::Function
end


# Temperature access
function _set_T(tbase::thermo_base, T::Float64)
    sym=Libdl.dlsym(tbase.ct_lib,:thermo_setTemperature)
    ccall(sym,
        Cint,(Cint,Cdouble),
        tbase.ind,T)
    return nothing
end

function _get_T(tbase::thermo_base)
    sym=Libdl.dlsym(tbase.ct_lib,:thermo_temperature)
    T=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return T
end


# working with names and indices
 

# not working right now
function _get_spec_index(tbase::thermo_base,spec_name::String)
    sym=Libdl.dlsym(tbase.ct_lib,:thermo_getName)
    ind=ccall(sym,
        Cint,(Cint,Csize_t,Cstring),
        tbase.ind,length(spec_name),spec_name)
    return ind
end

# species mole fractions
# By name, X has format "A:1,B:2,C:3"
function _set_X(tbase::thermo_base, X::String)
    sym=Libdl.dlsym(tbase.ct_lib,:thermo_setMoleFractionsByName)
    ccall(sym,
        Cint,(Cint,Cstring),
        tbase.ind,X)
    return nothing
end
# set all species
function _set_X(tbase::thermo_base,X::Array{Float64,1},norm=true)
    # norm is a boolean flag for whether to normalize
    sym=Libdl.dlsym(tbase.ct_lib,:thermo_setMoleFractions)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble},Cint),
        tbase.ind,tbase.nspec,X, norm)
    return nothing
end

# get all species
function _get_X(tbase::thermo_base)
    sym=Libdl.dlsym(tbase.ct_lib,:thermo_getMoleFractions)
    X=Array{Float64,1}(undef,tbase.nspec)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble}),
        tbase.ind,tbase.nspec,X)
    return X
end

# species mass fractions
# By name, Y has format "A:1,B:2,C:3"
function _set_Y(tbase::thermo_base, Y::String)
    sym=Libdl.dlsym(tbase.ct_lib,:thermo_setMassFractionsByName)
    ccall(sym,
        Cint,(Cint,Cstring),
        tbase.ind,Y)
    return nothing
end
# set all species
function _set_Y(tbase::thermo_base,Y::Array{Float64,1},norm=true)
    # norm is a boolean flag for whether to normalize
    sym=Libdl.dlsym(tbase.ct_lib,:thermo_setMassFractions)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble},Cint),
        tbase.ind,tbase.nspec,Y, norm)
    return nothing
end

# get all species
function _get_Y(tbase::thermo_base)
    sym=Libdl.dlsym(tbase.ct_lib,:thermo_getMassFractions)
    X=Array{Float64,1}(undef,tbase.nspec)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble}),
        tbase.ind,tbase.nspec,Y)
    return Y
end

# get thermo properties
# density
function _get_rho(tbase::thermo_base)
    sym=Libdl.dlsym(tbase.ct_lib,:thermo_density)
    rho=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return rho
end

# cp
function _get_cp(tbase::thermo_base)
    sym=Libdl.dlsym(tbase.ct_lib,:thermo_cp_mass)
    cp=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return cp
end

# molar cp
function _get_cp_mole(tbase::thermo_base)
    sym=Libdl.dlsym(tbase.ct_lib,:thermo_cp_mole)
    cp=ccall(sym,
        Float64,(Cint,),
        tbase.ind)
    return cp
end


# create the mutable thermo object
# if full path is true, must give entire path to the input file
# otherwise, searches in usr/local/share/cantera/data/
function thermo(tbase::thermo_base)
    # load in the file
    # functions for acccessing temperature
    Tsetter(T)=_set_T(tbase,T)
    Tgetter()=_get_T(tbase)
    # functions for acccessing pressure
    Psetter(P)=_set_P(tbase,P)
    Pgetter()=_get_P(tbase)
    # functions for accessing mole fractions
    Xsetter(X::Union{String,Array{Float64,1}})=_set_X(tbase,X)
    Xgetter()=_get_X(tbase)
    # functions for accessing mass fractions
    Ysetter(Y::Union{String,Array{Float64,1}})=_set_Y(tbase,Y)
    Ygetter()=_get_Y(tbase)
    # state setting functions
    function TPX_set(TPX::Tuple)
        Tsetter(TPX[1])
        Xsetter(TPX[3])
        Psetter(TPX[2])
        return nothing
    end
    function TPY_set(TPY::Tuple)
        Tsetter(TPY[1])
        Ysetter(TPY[3])
        Psetter(TPY[2])
        return nothing
    end
    # thermo properties functions
    rho_getter()=_get_rho(tbase)
    cp_getter()=_get_cp(tbase)
    cp_mole_getter()=_get_cp_mole(tbase)
    return thermo(tbase,
        Tsetter,Tgetter,
        Psetter, Pgetter,
        Xsetter,Xgetter,
        Ysetter,Ygetter,
        TPX_set,TPY_set,
        rho_getter,
        cp_getter, cp_mole_getter
        )
end

=#