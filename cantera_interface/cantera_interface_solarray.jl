# equivalent of python "SolutionArray" object
# I'm building functionality as I need it
mutable struct solutionArray
    Nel::Int
    gas::gas
    T::Array{Float64,1}
    P::Array{Float64,1}
    Y::Array{Float64,2}
    X::Array{Float64,2}
    λ::Array{Float64,1}
    cp::Array{Float64,1}
    R::Array{Float64,1}
    rho::Array{Float64,1}
    h::Array{Float64,1}
    e::Array{Float64,1}
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
    P=ones(Nel)*init_TPX[2]
    Y=Array{Float64,2}(undef,Nel,Nspec)
    X=Array{Float64,2}(undef,Nel,Nspec)
    λ=Array{Float64,1}(undef,Nel)
    cp=Array{Float64,1}(undef,Nel)
    R=Array{Float64,1}(undef,Nel)
    rho=Array{Float64,1}(undef,Nel)
    h=Array{Float64,1}(undef,Nel)
    e=Array{Float64,1}(undef,Nel)
    spec_MW=get_MW(g.phase)
    mean_MW=Array{Float64,1}(undef,Nel)
    partial_molar_cp=Array{Float64,2}(undef,Nel,Nspec)
    D_mix=Array{Float64,2}(undef,Nel,Nspec)
    net_production_rates=Array{Float64,2}(undef,Nel,Nspec)
    hdot=Array{Float64,1}(undef,Nel)    
    S=solutionArray(Nel, g,T,P,Y,X,λ,cp,R,rho,h,e,
    spec_MW,mean_MW,partial_molar_cp,D_mix,
    net_production_rates, hdot)
    for i=1:Nel
        _fetch_properties(S,i)
    end
    return S
end

function _fetch_properties(S::solutionArray,ind::Int,thermal_only::Bool=false)
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
    S.h[ind]=get_h(S.gas)
    S.e[ind]=get_e(S.gas)
    S.P[ind]=get_P(S.gas)
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

function set_T(S::solutionArray,ind::Int,
    T::Float64;thermal_only::Bool=false)
    set_T(S.gas,T)
    _fetch_properties(S,ind,thermal_only)
    return nothing    
end

function set_T(S::solutionArray,ind::Array{Int,1},
    T::Array{Float64,1};thermal_only::Bool=false)
    for i=1:length(ind)
        set_T(S.gas,T[ind[i]])
        _fetch_properties(S,ind[i],thermal_only)
    end
    return nothing    
end

function set_T(S::solutionArray,
    T::Array{Float64,1};thermal_only::Bool=false)
    for i=1:S.Nel
        set_T(S.gas,T[i])
        _fetch_properties(S,i,thermal_only)
    end
    return nothing    
end

function set_X(S::solutionArray,ind::Int,
    X::Array{Float64,1};thermal_only::Bool=false)
    set_X(S.gas,X)
    _fetch_properties(S,ind,thermal_only)
    return nothing    
end

function set_X(S::solutionArray,ind::Array{Int,1},
    X::Array{Float64,2};thermal_only::Bool=false)
    for i=1:length(ind)
        set_X(S.gas,X[ind[i],:])
        _fetch_properties(S,ind[i],thermal_only)
    end
    return nothing    
end

function set_X(S::solutionArray,X::Array{Float64,2};
    thermal_only::Bool=false)
    for i=1:S.Nel
        set_X(S.gas,X[i,:])
        _fetch_properties(S,i,thermal_only)
    end
end

function set_Y(S::solutionArray,ind::Int,Y::Array{Float64,1};
    thermal_only::Bool=false)
    set_Y(S.gas,Y)
    _fetch_properties(S,ind,thermal_only)
    return nothing    
end

function set_Y(S::solutionArray,ind::Array{Int,1},Y::Array{Float64,2};
    thermal_only::Bool=false)
    for i=1:length(ind)
        set_Y(S.gas,Y[ind[i],:])
        _fetch_properties(S,ind[i],thermal_only)
    end
    return nothing    
end

function set_Y(S::solutionArray,Y::Array{Float64,2};
    thermal_only::Bool=false)
    for i=1:S.Nel
        set_Y(S.gas,Y[i,:])
        _fetch_properties(S,i,thermal_only)
    end
end

function set_TPX(S::solutionArray,ind::Int,
    TPX::Tuple{Float64,Float64,Union{String,Array{Float64,1}}};thermal_only::Bool=false)
    set_TPX(S.gas,TPX)
    _fetch_properties(S,ind,thermal_only)
    return nothing    
end

function set_TPX(S::solutionArray,ind::Array{Int,1},
    TPX::Array{Tuple{Float64,Float64,Array{Float64,1}},1};thermal_only::Bool=false)
    for i=1:length(ind)
        set_TPX(S.gas,TPX[ind[i]])
        _fetch_properties(S,ind[i],thermal_only)
    end
    return nothing    
end

function set_TPX(S::solutionArray,
    TPX::Array{Tuple{Float64,Float64,Array{Float64,1}},1};thermal_only::Bool=false)
    for i=1:S.Nel
        set_TPX(S.gas,TPX[i])
        _fetch_properties(S,i,thermal_only)
    end
end

function set_TPY(S::solutionArray,ind::Int,
    TPY::Tuple{Float64,Float64,Union{String,Array{Float64,1}}};thermal_only::Bool=false)
    set_TPY(S.gas,TPY)
    _fetch_properties(S,ind,thermal_only)
    return nothing    
end

function set_TPY(S::solutionArray,ind::Array{Int,1},
    TPY::Array{Tuple{Float64,Float64,Array{Float64,1}},1};thermal_only::Bool=false)
    for i=1:length(ind)
        set_TPY(S.gas,TPY[ind[i]])
        _fetch_properties(S,ind[i],thermal_only)
    end
    return nothing    
end

function set_TPY(S::solutionArray,
    TPY::Array{Tuple{Float64,Float64,Array{Float64,1}},1};thermal_only::Bool=false)
    for i=1:S.Nel
        set_TPY(S.gas,TPY[i])
        _fetch_properties(S,i,thermal_only)
    end
end

# property setting via HPY
# single index
function set_HPY(S::solutionArray,ind::Int,
    HPY::Tuple{Float64,Float64,Union{String,Array{Float64,1}}};thermal_only::Bool=false)
    set_HPY(S.gas,HPY)
    _fetch_properties(S,ind,thermal_only)
    return nothing    
end

# multi index
function set_HPY(S::solutionArray,ind::Array{Int,1},
    HPY::Array{Tuple{Float64,Float64,Array{Float64,1}},1};thermal_only::Bool=false)
    for i=1:length(ind)
        set_HPY(S.gas,HPY[ind[i]])
        _fetch_properties(S,ind[i],thermal_only)
    end
    return nothing    
end

# all elements
function set_HPY(S::solutionArray,
    HPY::Array{Tuple{Float64,Float64,Array{Float64,1}},1};thermal_only::Bool=false)
    for i=1:S.Nel
        set_HPY(S.gas,HPY[i])
        _fetch_properties(S,i,thermal_only)
    end
end

# property setting via HPX
# single index
function set_HPX(S::solutionArray,ind::Int,
    HPX::Tuple{Float64,Float64,Union{String,Array{Float64,1}}};thermal_only::Bool=false)
    set_HPY(S.gas,HPX)
    _fetch_properties(S,ind,thermal_only)
    return nothing    
end

# multi index
function set_HPX(S::solutionArray,ind::Array{Int,1},
    HPX::Array{Tuple{Float64,Float64,Array{Float64,1}},1};thermal_only::Bool=false)
    for i=1:length(ind)
        set_HPX(S.gas,HPX[ind[i]])
        _fetch_properties(S,ind[i],thermal_only)
    end
    return nothing    
end

# all elements
function set_HPX(S::solutionArray,
    HPX::Array{Tuple{Float64,Float64,Array{Float64,1}},1};thermal_only::Bool=false)
    for i=1:S.Nel
        set_HPX(S.gas,HPX[i])
        _fetch_properties(S,i,thermal_only)
    end
end

# property setting via ERY
# single index
function set_ERY(S::solutionArray,ind::Int,
    ERY::Tuple{Float64,Float64,Union{String,Array{Float64,1}}};thermal_only::Bool=false)
    set_ERY(S.gas,ERY)
    _fetch_properties(S,ind,thermal_only)
    return nothing    
end

# multi index
function set_ERY(S::solutionArray,ind::Array{Int,1},
    ERY::Array{Tuple{Float64,Float64,Array{Float64,1}},1};thermal_only::Bool=false)
    for i=1:length(ind)
        set_ERY(S.gas,ERY[ind[i]])
        _fetch_properties(S,ind[i],thermal_only)
    end
    return nothing    
end

# all elements
function set_ERY(S::solutionArray,
    ERY::Array{Tuple{Float64,Float64,Array{Float64,1}},1};thermal_only::Bool=false)
    for i=1:S.Nel
        set_HPY(S.gas,ERY[i])
        _fetch_properties(S,i,thermal_only)
    end
end

# property setting via ERX
# single index
function set_ERX(S::solutionArray,ind::Int,
    ERX::Tuple{Float64,Float64,Union{String,Array{Float64,1}}};thermal_only::Bool=false)
    set_ERX(S.gas,ERX)
    _fetch_properties(S,ind,thermal_only)
    return nothing    
end

# multi index
function set_ERX(S::solutionArray,ind::Array{Int,1},
    ERX::Array{Tuple{Float64,Float64,Array{Float64,1}},1};thermal_only::Bool=false)
    for i=1:length(ind)
        set_ERX(S.gas,ERX[ind[i]])
        _fetch_properties(S,ind[i],thermal_only)
    end
    return nothing    
end

# all elements
function set_ERX(S::solutionArray,
    ERX::Array{Tuple{Float64,Float64,Array{Float64,1}},1};thermal_only::Bool=false)
    for i=1:S.Nel
        set_HPY(S.gas,ERX[i])
        _fetch_properties(S,i,thermal_only)
    end
end

function change_size!(S::solutionArray,Nnew::Int)
    fields=fieldnames(typeof(S))
    for k in 1:length(fields)
        tmp=getfield(S,k)
        if isa(tmp,Array) && size(tmp,1)==S.Nel
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
    return nothing
end