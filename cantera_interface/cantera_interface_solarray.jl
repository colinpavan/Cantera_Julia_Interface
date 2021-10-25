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
    fields=fieldnames(typeof(S))
    for k in 1:length(fields)
        tmp=getfield(S,k)
        if isa(tmp,Array)
            if size(tmp,2)==1
                tmp_new=Array{Float64,1}(undef,Nnew)
            else
                tmp_new=Array{Float64,2}(undef,Nnew,size(tmp,2))
            end
            for j=1:Nnew
                tmp_new[j,:].=tmp[mod(j-1,S.Nel)+1,:]
            end
            setfield!(S,fields[k],tmp_new)
        end
    end
    S.Nel=Nnew
end