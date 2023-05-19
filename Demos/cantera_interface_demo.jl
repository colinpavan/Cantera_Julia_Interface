#######################################################
#=
Demo of how to use Julia Cantera interface
C. Pavan 2023-05-10


=#
#######################################################

# Load Cantera module

cd(@__DIR__)
using Pkg
Pkg.activate(".")
push!(Base.DL_LOAD_PATH,"/usr/local/lib/") # change this to wherever libcantera_shared.so is on your machine)
include("../cantera_interface.jl")
using .cantera
ct=cantera

# general utility for pulling c++ errors
ct.ct_error_get() # check for errors

# Add in a gas

# gas with direct cantera interface
gas1=ct.gas("gri30.yaml", full_path=false) # if full_path = false, looks where cantera stores its files normally
# same as above with initialization of TPX
gas2=ct.gas("gri30.yaml",(1500.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))


# Get/Set State Demos

# Set by TPX (Temp/Pressure/Mole Fraction)
ct.set_TPX(gas1,(800.0,101.325e3,"N2:0.79,O2:0.21"))
# set by ERY (Internal energy/Density/Mass Fraction)
ct.set_ERY(gas2, (3e5, 1.0, "O2:1"))
# set temp only
ct.set_T(gas1,900.0)

#=
Note on setting order (if setting thermo states individually):
Cantera treats temperature and density as independent variables.
The only variables you can set in isolation are temp, density (rho) and pressure
Changing pressure changes density at constant temperature
Changing h and e can only be done in while also specifying P or rho
=#

ct.get_T(gas2)
ct.get_rho(gas2)
ct.set_P(gas2,2e5)
ct.get_T(gas2)
ct.get_rho(gas2)

ct.set_HP(gas2,(7e5,50e3))
ct.get_T(gas2)
ct.get_rho(gas2)

# mole/mass fractions can be specified either as a string (above) or as a vector
# get the names of the species in the order they are stored
specs=gas2.spec_names
X=ct.get_X(gas1)
# find CH4 index and change the molar concentration of that species
CH4_ind=findfirst("CH4" .== specs)
X[CH4_ind]=0.05
ct.set_X(gas2,X)


############# Transport properties #############

λ=ct.get_λ(gas2) # thermal cond.
D=ct.get_D_mix(gas2) # mixture averaged diffusivity coefficinets for each species


########### Kinetics properties #############

wdot=ct.net_production_rates(gas2) # production rates of each species`
hdot=ct.hdot(gas2) # net enthalpy release (per unit mass)



########### 0D Reactor #############
# compare timing to calling through python
using BenchmarkTools
@btime begin
gas2=ct.gas("gri30.yaml",(1500.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))
r=ct.reactor(gas2,"ConstPressureReactor")
tmax=1e-1
ct.reactor_sim_full(r,1e-1)
global Tfinal=ct.get_T()
end

# python version of same thing -> only works if you have python installed
using PyCall
gas2=ct.gas("gri30.yaml",(1500.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))
ctp=pyimport("cantera")
@btime begin
gasp=ctp.Solution("gri30.yaml")
rp=ctp.ConstPressureReactor()
gasp.TPX=ct.get_T(gas2),ct.get_P(gas2),ct.get_X(gas2)
rp.insert(gas)
sim=ctp.ReactorNet([r3])
sim.advance(1e-1)
tmp=sim.get_state()
end

# check same results


# "Solution Array" object holds state locally for an arbitray number of points
Nel=51
T=collect(range(300,350,length=Nel))
gas_array=ct.solutionArray("gri30.yaml",Nel,(300.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))
# compare access times - faster to use locally saved one
@btime gas_array.T[1] # value saved in julia
@btime ct.get_T(gas2) # get value from c++ library

# change size of array
ct.change_size!(gas_array,200)

# set state
ct.set_T(gas_array,1, 500.0)
# optional keyword 'thermal-only' avoids recalculating/storing the kinetics (slow)
# these are slower than directly communication with the gas object (i.e. ct.set_T) because
# it both sets the state and then loads it into julia
# useful if you will be using many of the state variables
@btime ct.set_T(gas_array, 5,600.0)
@btime ct.set_T(gas_array, 5,600.0, thermal_only=true)

