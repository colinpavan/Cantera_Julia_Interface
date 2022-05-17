# immutable base object
struct kin_base
    ind::Int
    nreact::Int
    nspec::Int
end

# TNumber of Reactions
function _nreact(rind::Integer)
    sym=Libdl.dlsym(lib,:kin_nReactions)
    nr=ccall(sym,
        Cint,(Cint,),
        rind)
    return nr
end
    
# create the kinetics base
# if full path is true, must give entire path to the input file
# otherwise, searches in usr/local/share/cantera/data/
function kin_base(file::String, thermo::thermo_base; full_path::Bool=false)
    # load in the file
    sym=Libdl.dlsym(lib,:kin_newFromFile)
    if !full_path
        file="/usr/local/share/cantera/data/" * file
    end
    kin_obj=ccall(sym,
        Cint, (Cstring, Cstring, Cint, Cint, Cint, Cint, Cint) ,
        file,"", thermo.ind,-1,-1,-1,-1)
    nreact = _nreact(kin_obj)
    return kin_base(kin_obj,
        nreact, thermo.nspec)
end

# kinetics access opjects

function net_rates_of_progress(kbase::kin_base)
    sym=Libdl.dlsym(lib,:kin_getNetRatesOfProgress)
    netROP=Array{Float64,1}(undef,kbase.nreact)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble}),
        kbase.ind, kbase.nreact,netROP)
    return netROP
end

function net_production_rates(kbase::kin_base)
    sym=Libdl.dlsym(lib,:kin_getNetProductionRates)
    netprod=Array{Float64,1}(undef,kbase.nspec)
    ccall(sym,
        Cint,(Cint,Csize_t,Ptr{Cdouble}),
        kbase.ind, kbase.nspec,netprod)
    return netprod
end

function _delta(kbase::kin_base, ind::Int)
    sym=Libdl.dlsym(lib,:kin_getDelta)
    δ=Array{Float64,1}(undef,kbase.nreact)
    ccall(sym,
        Cint,(Cint,Cint,Csize_t,Ptr{Cdouble}),
        kbase.ind, ind,kbase.nreact,δ)
    return δ
end

delta_enthalpy(kbase::kin_base)=_delta(kbase,0)
delta_gibbs(kbase::kin_base)=_delta(kbase,1)
delta_entropy(kbase::kin_base)=_delta(kbase,2)