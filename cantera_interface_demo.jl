cd(@__DIR__)
using Pkg
Pkg.activate(".")

if (:cantera in names(Main))
    using Libdl
    Libdl.dlclose(cantera.lib)
end
ct=include("./cantera_interface/cantera_interface.jl")

# gas array with direct cantera interface
gas1=ct.gas("gri30.yaml")
gas2=ct.gas("gri30.yaml",(300.0,cantera.one_atm,"CH4:1,O2:2,N2:7.52"))