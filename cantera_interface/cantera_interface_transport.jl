# immutable base object
struct trans_base
    ind::Int
    nspec::Int
end
# create the kinetics base
function trans_base(thermo::thermo_base; loglevel=0)
    sym=Libdl.dlsym(lib,:trans_new)
    trans_obj=ccall(sym,
        Cint, (Cstring, Cint,Cint) ,
        "Mix",thermo.ind,loglevel)
    return trans_base(trans_obj, thermo.nspec)
end

function thermal_cond(tbase::trans_base)
    sym=Libdl.dlsym(lib,:trans_thermalConductivity)
    k=ccall(sym,
        Cdouble, (Cint,),
        tbase.ind)
    return k
end

function mixture_diff_coeffs(tbase::trans_base)
    sym=Libdl.dlsym(lib,:trans_getMixDiffCoeffs)
    D=Array{Float64,1}(undef,tbase.nspec)
    k=ccall(sym,
        Cdouble, (Cint,Cint,Ptr{Cdouble}),
        tbase.ind, tbase.nspec,D)
    return D
end
