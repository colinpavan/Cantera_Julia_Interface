# equivalent of python "SolutionArray" object
# I'm building functionality as I need it
mutable struct solutionArray
    Nel::Int
    gas::gas
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

function solutionArray(file::String, Nel::Int,
    init_TPX::Tuple{Float64,Float64,Union{String,Float64}}=nothing;
    full_path::Bool=false)
    if isnothing(init_TPX)
        g=gas(file, full_path=full_path)
    else
        g=gas(file,init_TPX, full_path=full_path)
    end
    Nspec=g.Nspec
    T=Array{Float64,1}(undef,Nel)
    P=0.0
    Y=Array{Float64,2}(undef,Nel,Nspec)
    X=Array{Float64,2}(undef,Nel,Nspec)
    λ=Array{Float64,1}(undef,Nel)
    cp=Array{Float64,1}(undef,Nel)
    R=Array{Float64,1}(undef,Nel)
    rho=Array{Float64,1}(undef,Nel)
    spec_MW=Array{Float64,1}(undef,Nel)
    mean_MW=Array{Float64,1}(undef,Nel)
    partial_molar_cp=Array{Float64,2}(undef,Nel,Nspec)
    D_mix=Array{Float64,2}(undef,Nel,Nspec)
    net_production_rates=Array{Float64,2}(undef,Nel,Nspec)
    hdot=Array{Float64,1}(undef,Nel)    
    S=solutionArray(Nel, g,T,P,Y,X,λ,cp,R,rho,
    spec_MW,mean_MW,partial_molar_cp,D_mix,
    net_production_rates, hdot)
    for i=1:Nel
        _update_gas_array_single(S,i,get_TPX(g))
    end
    return S
end

function _update_gas_array_single(S::solutionArray, ind::Int, 
    TPX::Tuple{Float64,Float64,Union{String,Array{Float64,1}}}; thermal_only=false)
    set_TPX(S.gas,TPX)
    _get_thermo_props!(S,ind)
    _get_trans_props!(S,ind)
    if !thermal_only
        _get_kin_props!(S,ind)
    end
end
function _update_gas_array_single(S::solutionArray, ind::Int, 
    T::Float64; thermal_only=false)
    set_T(S.gas,T)
    _get_thermo_props!(S,ind)
    _get_trans_props!(S,ind)
    if !thermal_only
        _get_kin_props!(S,ind)
    end
end
function _update_gas_array_single(S::solutionArray, ind::Int, 
    X::Array{Float64,1}; thermal_only=false)
    set_X(S.gas,X)
    _get_thermo_props!(S,ind)
    _get_trans_props!(S,ind)
    if !thermal_only
        _get_kin_props!(S,ind)
    end
end

function _get_thermo_props!(S::solutionArray,ind)
    S.T[ind]=get_T(S.gas)
    S.X[ind,:].=get_X(S.gas)
    S.Y[ind,:].=get_Y(S.gas)
    S.cp[ind]=get_cp(S.gas)
    S.rho[ind]=get_rho(S.gas)
    S.mean_MW[ind]=get_mean_MW(S.gas)
    S.partial_molar_cp[ind,:].=get_spec_molar_cp(S.gas) 
    S.R[ind]=Ru./S.mean_MW[ind]
end
function _get_trans_props!(S::solutionArray,ind)
    S.λ[ind]=get_λ(S.gas)
    S.D_mix[ind,:].=mixture_diff_coeffs(S.gas.trans)
end
function _get_kin_props!(S::solutionArray,ind)
    S.net_production_rates[ind,:].=net_production_rates(S.gas.kin)
    S.hdot[ind]=sum(net_rates_of_progress(S.gas.kin).*
    delta_enthalpy(S.gas.kin))
end



function update_gas_array_partial(S::solutionArray,ind::Union{Array{Int,1},Int},
    TPX::Union{Array{Tuple,1},Tuple{Float64,Float64,Union{String,Float64}}})
    if typeof(TPX)==tuple
        _update_gas_array_single(S,ind,TPX)
    else
        for i=1:length(Ind)
            _update_gas_array_single(S,ind[i],TPX[i])
        end
    end
end

function update_gas_array(S::solutionArray,TPX::Array{Tuple,1})
    for i=1:S.Nel
        _update_gas_array_single(S,i,TPX)
    end
end

function set_T(S::solutionArray,ind::Int,T::Float64;thermal_only::Bool=false)
    _update_gas_array_single(S,ind,T, thermal_only=thermal_only)
    return nothing    
end
function set_T(S::solutionArray,ind::Array{Int,1},T::Array{Float64,1};thermal_only::Bool=false)
    for i=1:length(ind)
        _update_gas_array_single(S,ind[i],T[ind[i]], thermal_only=thermal_only)
    end
    return nothing    
end
function set_T(S::solutionArray,T::Array{Float64,1};thermal_only::Bool=false)
    for i=1:S.Nel
        _update_gas_array_single(S,i,T[i], thermal_only=thermal_only)
    end
    return nothing    
end

function set_X(S::solutionArray,ind::Int,X::Array{Float64,1})
    _update_gas_array_single(S,ind,X)
    return nothing    
end
function set_X(S::solutionArray,ind::Array{Int,1},X::Array{Float64,2})
    for i=1:length(ind)
        _update_gas_array_single(S,ind[i],X[ind[i],:])
    end
    return nothing    
end
function set_X(S::solutionArray,X::Array{Float64,2})
    for i=1:S.Nel
        _update_gas_array_single(S,i,X[i,:])
    end
    return nothing    
end

function change_size!(S::solutionArray,Nnew::Int)
    if Nnew<S.Nel
        S.phase=S.phase[1:Nnew]
        S.kin=S.kin[1:Nnew]
        S.trans=S.trans[1:Nnew]
    elseif  Nnew>S.Nel
        add=Nnew-S.Nel
        phase=Array{thermo_base,1}(undef,add)
        kin=Array{kin_base,1}(undef,add)
        trans=Array{trans_base,1}(undef,add)
        for i=1:add
            phase[i]=thermo_base(S.file)
            kin[i]=kin_base(S.file, phase[i])
            trans[i]=trans_base(phase[i])
        end
        append!(S.phase,phase)
        append!(S.kin,kin)
        append!(S.trans,trans)
    end
    S.Nel=Nnew
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

function change_size!(SL::solutionArray_local, Nnew::Int)
    change_size!(SL.S, Nnew)
    SL=solutionArray_local(SL.S)
end