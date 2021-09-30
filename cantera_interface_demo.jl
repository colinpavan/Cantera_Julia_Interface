cd(@__DIR__)
using Pkg
Pkg.activate(".")

if !(:cantera in names(Main))
    include("./cantera_interface/cantera_interface.jl")
    using .cantera
end

# gas array with direct cantera interface
gas=cantera.thermo_base("gri30.yaml")
cantera.set_T(gas,300.0)
