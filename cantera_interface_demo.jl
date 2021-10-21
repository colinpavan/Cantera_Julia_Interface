cd(@__DIR__)
using Pkg
Pkg.activate(".")
# to use this, the library libcantera_shared.so must be on your system search path
ct=include("./cantera_interface/cantera_interface.jl")

# gas with direct cantera interface
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

# gas array - usefull for 1D simulations
Nel=50
gas_array=ct.solutionArray("gri30.yaml",Nel,(300.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))
# local copy for fast access
gas_array_local=ct.solutionArray_local(gas_array)
[ct.set_T(gas_array.phase[i], 300.0+i) for i=1:Nel]
ct.update_gas_array(gas_array_local)
# compare timing of local vs library copies - 2 order of magnitude difference!
using BenchmarkTools
@btime gas_array_local.T
@btime ct.get_T(gas_array)
# if you are going to access a property more than once
# it will always be more efficient to store it in the local object

# change size of gas array
ct.change_size!(gas_array,60)
ct.change_size!(gas_array,40)
ct.change_size!(gas_array_local, 70)