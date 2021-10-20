cd(@__DIR__)
using Pkg
Pkg.activate(".")
# to use this, the library libcantera_shared.so must be on your system search path
ct=include("./cantera_interface/cantera_interface.jl")

# gas array with direct cantera interface
gas1=ct.gas("gri30.yaml")
# same as above with initialization of TPX
gas2=ct.gas("gri30.yaml",(300.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))
# this stores a local copy (inside julia) for faster access
gas2_loc=ct.gas_local(gas2)
# change the temperature
ct.set_T(gas2,500.0)
# note this isn't changed in the local copy
gas2_loc.T
# update the local copy
ct.update_gas(gas2_loc)
# now the local copy is updated
gas2_loc.T
gas2_loc.X