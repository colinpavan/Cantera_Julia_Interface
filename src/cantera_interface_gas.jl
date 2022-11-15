#####################################################
#=
Interface to cantera/src/clib/ct.cpp file
Calls 3 separate functions for dealing with thermo, kinetics and transport

C. Pavan 2022-10-28
=#
####################################################


####### Structure definition and initialization ######

struct gas
    Nspec::Int
    spec_names::Array{String,1}
    phase::thermo_base
    kin::kin_base
    trans::trans_base
end

# case for not initializing TPX
function gas(file::String; full_path::Bool=true)
    phase=thermo_base(file,full_path=full_path)
    kin=kin_base(file,phase, full_path=full_path)
    trans=trans_base(phase)
    Nspec=phase.nspec
    spec_names=[_get_spec_name(phase,i) for i=1:Nspec]
    return gas(Nspec,spec_names,phase,kin,trans)
end

# case for initializing TPX
function gas(file::String, init_TPX::Tuple{Float64,Float64,Union{String,Float64}};
    full_path::Bool=true)
    phase=thermo_base(file,full_path=full_path)
    kin=kin_base(file,phase, full_path=full_path)
    trans=trans_base(phase)
    set_TPX(phase,init_TPX)
    Nspec=phase.nspec
    spec_names=[_get_spec_name(phase,i) for i=1:Nspec]
    return gas(Nspec,spec_names,phase,kin,trans)
end

############### Thermo setter functions #################
# temperature
function set_T(G::gas,T::Float64)
    set_T(G.phase,T)
end

# pressure
function set_P(G::gas,P::Float64)
    set_P(G.phase,P)
end

# pressure
function set_rho(G::gas,rho::Float64)
    set_P(G.phase,rho)
end

############# Composition setter functions ###########3
# by mole fraction
function set_X(G::gas,X::AbstractVector{Float64}, norm::Bool=true)
    set_X(G.phase,X, norm)
end

############## Multi-thermo setter functions ############

function set_HP(G::gas,HP::Tuple{Float64,Float64})
    return set_HP(G.phase, HP)
end

function set_TR(G::gas,TR::Tuple{Float64,Float64})
    return set_TR(G.phase, TR)
end

function set_TPX(G::gas,TPX::Tuple{Float64,Float64,Union{String,Array{Float64,1}}})
    set_TPX(G.phase,TPX)
end

function set_TPY(G::gas,TPY::Tuple{Float64,Float64,Union{String,Array{Float64,1}}})
    set_TPY(G.phase,TPY)
end

function set_HPY(G::gas,HPY::Tuple{Float64,Float64,Union{String,Array{Float64,1}}})
    set_HPY(G.phase,HPY)
end
function set_HPX(G::gas,HPX::Tuple{Float64,Float64,Union{String,Array{Float64,1}}})
    set_HPY(G.phase,HPX)
end

function set_ERY(G::gas,ERY::Tuple{Float64,Float64,Union{String,Array{Float64,1}}})
    set_ERY(G.phase,ERY)
end
function set_ERX(G::gas,ERX::Tuple{Float64,Float64,Union{String,Array{Float64,1}}})
    set_ERX(G.phase,ERX)
end

############### thermo getter functions #########3
# temperature
function get_T(G::gas)
    return get_T(G.phase)
end

# pressure
function get_P(G::gas)
    return get_P(G.phase)
end

# cp
function get_cp(G::gas)
    return get_cp(G.phase)
end

# cv
function get_cv(G::gas)
    return get_cv(G.phase)
end

# density
function get_rho(G::gas)
    return get_rho(G.phase)
end

# mean molecular weight
function get_mean_MW(G::gas)
    return get_mean_MW(G.phase)
end

# species molecular weights
function get_MW(G::gas)
    return get_MW(G.phase)
end

# species molar cp
function get_spec_molar_cp(G::gas)
    return get_spec_molar_cp(G.phase)
end

# species molar enthalpy
function get_spec_molar_enthalpies(G::gas)
    return get_spec_molar_enthalpies(G.phase)
end

# enthalpy
function enthalpy(G::gas)
    return get_h(G.phase)
end

# enthalpy (accessible by a different name)
function get_h(G::gas)
    return get_h(G.phase)
end

# internal energy
function get_e(G::gas)
    return get_e(G.phase)
end

function get_spec_heat_ratio(G::gas)
    return get_cp(G.phase)/get_cv(G.phase)
end

function get_n(G::gas)
   return get_X(G.phase)*(get_P(G.phase)/(1.380649e-23*get_T(G.phase))) 
end

function get_ncm3(G::gas)
    return get_X(G.phase)*(get_P(G.phase)/(1.380649e-17*get_T(G.phase))) 
 end
 

##########3 Composition getter functions ############
# by mole fraction
function get_X(G::gas)
    return get_X(G.phase)
end

#by mass fraction
function get_Y(G::gas)
    return get_Y(G.phase)
end

######### Transport getter functions ###########
# thermal conductivity
function get_λ(G::gas)
    return thermal_cond(G.trans)
end

# mixture averaged diffusion coefficiences
function get_D_mix(G::gas)
    return mixture_diff_coeffs(G.trans)
end

########3 Kinetics getter functions ###########

# net production rates
function net_production_rates(G::gas)
    return net_production_rates(G.kin)
end

# total enthalpy release
function hdot(G::gas)
    return sum(net_rates_of_progress(G.kin).*delta_enthalpy(G.kin))
end

########## multi-param thermo getters ##############33

function get_TPX(G)
    return (get_T(G),get_P(G),get_X(G))
end
function get_TPY(G)
    return (get_T(G),get_P(G),get_Y(G))
end

function get_ERX(G)
    return (get_e(G),get_rho(G),get_X(G))
end
function get_ERY(G)
    return (get_e(G),get_rho(G),get_Y(G))
end



# # this evaluates all the above and stores them in local variables
# # makes it so only need to access cantera on update
# # not on repeat calls
# # intentionally does NOT do all variables - only the ones needed for a flame
# mutable struct gas_local
#     G::gas
#     T::Float64
#     P::Float64
#     Y::Array{Float64,1}
#     X::Array{Float64,1}
#     λ::Float64
#     cp::Float64
#     R::Float64
#     rho::Float64
#     spec_MW::Array{Float64,1}
#     mean_MW::Float64
#     partial_molar_cp::Array{Float64,1}
#     D_mix::Array{Float64,1}
#     net_production_rates::Array{Float64,1}
#     hdot::Float64
# end

# function update_gas(GL::gas_local; thermal_only=false)
#     GL.T=get_T(GL.G)
#     GL.X.=get_X(GL.G)
#     GL.Y.=get_Y(GL.G)
#     GL.λ=get_λ(GL.G)
#     GL.cp=get_cp(GL.G)
#     GL.rho=get_rho(GL.G)
#     GL.mean_MW=get_mean_MW(GL.G)
#     GL.partial_molar_cp.=get_spec_molar_cp(GL.G)
#     GL.D_mix.=get_D_mix(GL.G)
#     if !thermal_only
#         GL.net_production_rates=net_production_rates(GL.G)
#         GL.hdot=hdot(GL.G)
#         GL.R=Ru./GL.mean_MW
#     end
#     return nothing
# end

# function gas_local(G::gas)
#     return gas_local(
#         G,
#         get_T(G),
#         get_P(G),
#         get_X(G),
#         get_Y(G),
#         get_λ(G),
#         get_cp(G),
#         get_cp(G).-get_cv(G),
#         get_rho(G),
#         get_MW(G),
#         get_mean_MW(G),
#         get_spec_molar_cp(G),
#         get_D_mix(G),
#         net_production_rates(G),
#         hdot(G)
#     )
# end