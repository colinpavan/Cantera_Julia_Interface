# equivalent of python "SolutionArray" object
# I'm building functionality as I need it
struct gas
    Nspec::Int
    spec_names::Array{String,1}
    phase::thermo_base
    kin::kin_base
    trans::trans_base
end

# case for not initializing TPX
function gas(file::String; full_path::Bool=false)
    phase=thermo_base(file,full_path=full_path)
    kin=kin_base(file,phase, full_path=full_path)
    trans=trans_base(phase)
    Nspec=phase.nspec
    spec_names=[_get_spec_name(phase,i) for i=1:Nspec]
    return gas(Nspec,spec_names,phase,kin,trans)
end

# case for initializing TPX
function gas(file::String, init_TPX::Tuple{Float64,Float64,Union{String,Float64}};
    full_path::Bool=false)
    phase=thermo_base(file,full_path=full_path)
    kin=kin_base(file,phase, full_path=full_path)
    trans=trans_base(phase)
    set_TPX(phase,init_TPX)
    Nspec=phase.nspec
    spec_names=[_get_spec_name(phase,i) for i=1:Nspec]
    return gas(Nspec,spec_names,phase,kin,trans)
end

function get_T(G::gas)
    return get_T(G.phase)
end
function set_T(G::gas,T)
    set_T(G.phase,T)
end

function get_P(G::gas)
    return get_P(G.phase)
end
function set_P(G::gas,T)
    set_P(G.phase,T)
end

function get_X(G::gas)
    return get_X(G.phase)
end

function set_X(G::gas,X)
    set_X(G.phase,X)
end

function get_Y(G::gas)
    return get_Y(G.phase)
end

function get_cp(G::gas)
    return get_cp(G.phase)
end

function get_cv(G::gas)
    return get_cv(G.phase)
end

function get_λ(G::gas)
    return thermal_cond(G.trans)
end

function get_rho(G::gas)
    return get_rho(G.phase)
end

function get_D_mix(G::gas)
    return mixture_diff_coeffs(G.trans)
end

function get_mean_MW(G::gas)
    return get_mean_MW(G.phase)
end

function get_MW(G::gas)
    return get_MW(G.phase)
end

function get_spec_molar_cp(G::gas)
    return get_spec_molar_cp(G.phase)
end

function net_production_rates(G::gas)
    return net_production_rates(G.kin)
end

function hdot(G::gas)
    return sum(net_rates_of_progress(G.kin).*delta_enthalpy(G.kin))
end

function enthalpy(G::gas)
    return get_h(G.phase)
end


# this evaluates all the above and stores them in local variables
# makes it so only need to access cantera on update
# not on repeat calls
# intentionally does NOT do all variables - only the ones needed for a flame
mutable struct gas_local
    G::gas
    T::Float64
    P::Float64
    Y::Array{Float64,1}
    X::Array{Float64,1}
    λ::Float64
    cp::Float64
    R::Float64
    rho::Float64
    spec_MW::Array{Float64,1}
    mean_MW::Float64
    partial_molar_cp::Array{Float64,1}
    D_mix::Array{Float64,1}
    net_production_rates::Array{Float64,1}
    hdot::Float64
end

function update_gas(GL::gas_local; thermal_only=false)
    GL.T=get_T(GL.G)
    GL.X.=get_X(GL.G)
    GL.Y.=get_Y(GL.G)
    GL.λ=get_λ(GL.G)
    GL.cp=get_cp(GL.G)
    GL.rho=get_rho(GL.G)
    GL.mean_MW=get_mean_MW(GL.G)
    GL.partial_molar_cp.=get_spec_molar_cp(GL.G)
    GL.D_mix.=get_D_mix(GL.G)
    if !thermal_only
        GL.net_production_rates=net_production_rates(GL.G)
        GL.hdot=hdot(GL.G)
        GL.R=Ru./GL.mean_MW
    end
    return nothing
end

function gas_local(G::gas)
    return gas_local(
        G,
        get_T(G),
        get_P(G),
        get_X(G),
        get_Y(G),
        get_λ(G),
        get_cp(G),
        get_cp(G).-get_cv(G),
        get_rho(G),
        get_MW(G),
        get_mean_MW(G),
        get_spec_molar_cp(G),
        get_D_mix(G),
        net_production_rates(G),
        hdot(G)
    )
end