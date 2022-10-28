cd(@__DIR__)
using Pkg
Pkg.activate("../.")
push!(Base.DL_LOAD_PATH,"/usr/local/lib/") # change this to wherever libcantera_shared.so is on your machine)
using Revise
include("../cantera_interface.jl")
using .cantera
ct=cantera


gas=ct.gas("../../../data/Comb_Mech/gri30.yaml")


function get_thermo_state(g::ct.gas)
    return [ct.get_rho(g),ct.get_T(g), ct.get_P(g), ct.get_h(g), ct.get_e(g)]
end

ct.set_TPX(gas, (300.0,ct.one_atm,"N2:1"))
s1=get_thermo_state(gas)

ct.set_X(gas, "O:1")
s2=get_thermo_state(gas)

ct.set_P(gas,50e3)
s3=get_thermo_state(gas)

ct.set_X(gas, "O2:1")
s4=get_thermo_state(gas)
