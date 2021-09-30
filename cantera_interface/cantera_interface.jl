module cantera
export one_atm

using Libdl
const one_atm=101325.0
const Ru=8314.46261815324
# this will find cantera provide libcantera_shared.so is on your system path
const lib=Libdl.dlopen("libcantera_shared.so")

include("cantera_interface_thermo.jl")
include("cantera_interface_kinetics.jl")
include("cantera_interface_transport.jl")
#include("cantera_interface_onedim.jl")


function ct_error_get()
    sym=Libdl.dlsym(lib,:ct_getCanteraError)
    tmp=""
    len=ccall(sym,
        Cint,(Cint, Cstring),
        0,tmp)
    tmp=String(['a' for i in 1:len])
    ccall(sym,
        Cint, (Cint,Cstring),
        len, tmp)
    print("$tmp\n")
    return nothing
end


# equivalent of python "SolutionArray" object
# I'm building functionality as I need it
struct solutionArray
    Nel::Int
    phase::Array{thermo_base,1}
    kin::Array{kin_base,1}
    trans::Array{trans_base,1}
end

function solutionArray(file::String, Nelement::Int,
    init_TPX::Tuple{Float64,Float64,Union{String,Float64}};
    full_path::Bool=false)
    phase=Array{thermo_base,1}(undef,Nelement)
    kin=Array{kin_base,1}(undef,Nelement)
    trans=Array{trans_base,1}(undef,Nelement)
    for i=1:Nelement
        phase[i]=thermo_base(file, full_path=full_path)
        kin[i]=kin_base(file,phase[i], full_path=full_path)
        trans[i]=trans_base(phase[i])
        set_TPX(phase[i],init_TPX)
    end
    return solutionArray(Nelement,phase,kin, trans)
end

function get_T(S::solutionArray)
    T=[get_T(c) for c in S.phase]
    return T
end

function get_X(S::solutionArray)
    X=Array{Float64,2}(undef,S.Nel,S.phase[1].nspec)
    for i=1:S.Nel
        X[i,:]=get_X(S.phase[i])
    end
    return X
end

function get_Y(S::solutionArray)
    Y=Array{Float64,2}(undef,S.Nel,S.phase[1].nspec)
    for i=1:S.Nel
        Y[i,:]=get_Y(S.phase[i])
    end
    return Y
end

function get_cp(S::solutionArray)
    cp=[get_cp(c) for c in S.phase]
    return cp
end

function get_cv(S::solutionArray)
    cv=[get_cv(c) for c in S.phase]
    return cv
end

function get_λ(S::solutionArray)
    λ=[thermal_cond(c) for c in S.trans]
    return λ
end

function get_rho(S::solutionArray)
    ρ=[get_rho(c) for c in S.phase]
    return ρ
end

function get_D_mix(S::solutionArray)
    D=Array{Float64,2}(undef,S.Nel,S.phase[1].nspec)
    for i=1:S.Nel
        D[i,:]=mixture_diff_coeffs(S.trans[i])
    end
    return D
end

function get_mean_MW(S::solutionArray)
    MW=[get_mean_MW(c) for c in S.phase]
    return MW
end

# this doesnt change with position, so just call once
function get_MW(S::solutionArray)
    return get_MW(S.phase[1])
end

function get_spec_molar_cp(S::solutionArray)
    cp=Array{Float64,2}(undef,S.Nel,S.phase[1].nspec)
    for i=1:S.Nel
        cp[i,:]=get_spec_molar_cp(S.phase[i])
    end
    return cp
end

function net_production_rates(S::solutionArray)
    wdot=Array{Float64,2}(undef,S.Nel,S.phase[1].nspec)
    for i=1:S.Nel
        wdot[i,:]=net_production_rates(S.kin[i])
    end
    return wdot
end

function hdot(S::solutionArray)
    hdot=[sum(net_rates_of_progress(k).*delta_enthalpy(k)) for k in S.kin]
    return hdot
end

function enthalpy(S::solutionArray)
    h=[get_h(c) for c in S.phase]
    return h
end


# this evaluates all the above and stores them in local variables
# makes it so only need to access cantera on update
# not on repeat calls
# intentionally does NOT do all variables - only the ones needed for a flame
struct solutionArray_local
    S::solutionArray
    T::Array{Float64,1}
    P::Float64
    Y::Array{Float64,2}
    X::Array{Float64,2}
    λ::Array{Float64,1}
    cp::Array{Float64,1}
    R::Array{Float64,1}
    rho::Array{Float64,1}
    spec_MW::Array{Float64,1}
    mean_MW::Array{Float64,1}
    partial_molar_cp::Array{Float64,2}
    D_mix::Array{Float64,2}
    net_production_rates::Array{Float64,2}
    hdot::Array{Float64,1}
end

function update_gas_array(SL::solutionArray_local; thermal_only=false)
    SL.T.=get_T(SL.S)
    SL.X.=get_X(SL.S)
    SL.Y.=get_Y(SL.S)
    SL.λ.=get_λ(SL.S)
    SL.cp.=get_cp(SL.S)
    SL.rho.=get_rho(SL.S)
    SL.mean_MW.=get_mean_MW(SL.S)
    SL.partial_molar_cp.=get_spec_molar_cp(SL.S)
    SL.D_mix.=get_D_mix(SL.S)
    if !thermal_only
        SL.net_production_rates.=net_production_rates(SL.S)
        SL.hdot.=hdot(SL.S)
        SL.R.=Ru./SL.mean_MW
    end
    return nothing
end

function update_gas_array_partial(SL::solutionArray_local, ind; thermal_only=false)
    for i in ind
        SL.T[i]=get_T(SL.S.phase[i])
        SL.X[i,:].=get_X(SL.S.phase[i])
        SL.Y[i,:].=get_Y(SL.S.phase[i])
        SL.λ[i]=thermal_cond(SL.S.trans[i])
        SL.cp[i]=get_cp(SL.S.phase[i])
        SL.rho[i]=get_rho(SL.S.phase[i])
        SL.mean_MW[i]=get_mean_MW(SL.S.phase[i])
        SL.partial_molar_cp[i,:].=get_spec_molar_cp(SL.S.phase[i])
        SL.D_mix[i,:].=mixture_diff_coeffs(SL.S.trans[i])
        if !thermal_only
            SL.net_production_rates[i,:].=net_production_rates(SL.S.kin[i])
            SL.hdot[i]=sum(net_rates_of_progress(SL.S.kin[i]).*
                delta_enthalpy(SL.S.kin[i]))
            SL.R[i]=Ru./SL.mean_MW[i]
        end
    end
    return nothing
end


function solutionArray_local(S::solutionArray)
    return solutionArray_local(
        S,
        get_T(S),
        get_P(S.phase[1]),
        get_X(S),
        get_Y(S),
        get_λ(S),
        get_cp(S),
        get_cp(S).-get_cv(S),
        get_rho(S),
        get_MW(S),
        get_mean_MW(S),
        get_spec_molar_cp(S),
        get_D_mix(S),
        net_production_rates(S),
        hdot(S)
    )
end
end
