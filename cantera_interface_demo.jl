cd(@__DIR__)
using Pkg
Pkg.activate(".")
# to use this, the library libcantera_shared.so must be on your system search path
ct=include("./cantera_interface/cantera_interface.jl")

# gas with direct cantera interface
gas1=ct.gas("gri30.yaml")
# same as above with initialization of TPX
gas2=ct.gas("gri30.yaml",(1500.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))
# testing 0D simulation

r=ct.reactor(gas2,"ConstPressureReactor")
ct.reactor_sim_full(r,1e-1)
ct.get_TPY(gas2)
ct.ct_error_get()

using BenchmarkTools
gas2=ct.gas("gri30.yaml",(1500.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))
r=ct.reactor(gas2,"IdealGasConstPressureReactor")
@btime ct.reactor_sim_full(r,1e-1)
gas2=ct.gas("gri30.yaml",(1500.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))
r2=ct.reactor(gas2,"ConstPressureReactor")
@btime ct.reactor_sim_full(r2,1e-1)
gas2=ct.gas("gri30.yaml",(1500.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))
r3=ct.reactor(gas2,"IdealGasConstPressureReactor")
RN=ct.reactor_sim_initialize(r3)
tmax=1e-1
@btime ct.reactor_sim_advance(RN,tmax)
using PyCall
gas2=ct.gas("gri30.yaml",(1500.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))
ctp=pyimport("cantera")
gas=ctp.Solution("gri30.yaml")
r3=ctp.IdealGasConstPressureReactor()
gas.TPX=ct.get_T(gas2),ct.get_P(gas2),ct.get_X(gas2)
r3.insert(gas)
sim=ctp.ReactorNet([r3])
@btime sim.advance(1e-1)
tmp=sim.get_state()

# compare timing of local vs library copies - 2 order of magnitude difference!
using BenchmarkTools
T=collect(range(300,350,length=Nel))
gas_array=ct.solutionArray("gri30.yaml",Nel,(300.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))
@btime gas_array=ct.solutionArray("gri30.yaml",Nel,(300.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))
@btime ct.set_T(gas_array,T)
@btime ct.set_T(gas_array,T, thermal_only=true)
@btime ct.set_T(gas_array,5,400.0)
@btime ct.set_T(gas_array,5,400.0, thermal_only=true)
@btime ct.set_TPX(gas_array,10,(400.0, ct.one_atm,"CH4:1,O2:2.5,N2:7.52"))
@btime ct.set_TPX(gas_array,10,(400.0, ct.one_atm,"CH4:1,O2:2.5,N2:7.52"), thermal_only=true)
X=gas_array.X[10,:]
@btime ct.set_TPX(gas_array,10,(400.0, ct.one_atm,X))
@btime ct.set_TPX(gas_array,10,(400.0, ct.one_atm,X), thermal_only=true)

ct.change_size!(gas_array,200)
ct.change_size!(gas_array,1000)
