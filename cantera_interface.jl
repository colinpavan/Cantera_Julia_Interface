module cantera
export one_atm, Ru

using Libdl
const one_atm=101325.0
const Ru=8314.46261815324
const kb=1.380649e-23
# this will find cantera provide libcantera_shared.so is on your system path
# if statement ensures library only opened once
if @isdefined lib
    Libdl.dlclose(lib)
end

const lib=Libdl.dlopen("libcantera_shared.so")

include("src/cantera_interface_thermo.jl")
include("src/cantera_interface_kinetics.jl")
include("src/cantera_interface_transport.jl")
include("src/cantera_interface_gas.jl")
include("src/cantera_interface_solarray.jl")
include("src/cantera_interface_zeroD.jl")


function ct_error_get()
    sym=Libdl.dlsym(lib,:ct_getCanteraError)
    tmp=""
    len=ccall(sym,
        Cint,(Cint, Cstring),
        0,tmp)
    tmp=String(['a' for i in 1:len])
    ccall(sym,
        Cint, (Cint,Cstring),
        len, tmp)
    print("$tmp\n")
    return nothing
end





end