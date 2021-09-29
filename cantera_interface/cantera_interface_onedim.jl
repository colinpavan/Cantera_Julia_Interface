# immutable base object
struct flame_base
    ct_lib::Ptr
    doms::Array{Int32,1}
end

#=
this clears all existing domains
    sym=Libdl.dlsym(ct.lib,:ct_clearOneDim)
    ccall(sym,
        Cint,(),
        )
=#
# initialization of inlet
function inlet_init(ct::ct_lib, inlet::Int32, phase::thermo)
    sym=Libdl.dlsym(ct.lib,:bdry_setTemperature)
    ccall(sym,
        Cint, (Cint, Cdouble),
        inlet, phase.get_T())
    return nothing
end

# initialization of stflow
function _set_flow_P(ct::ct_lib, flow_dom::Int32, phase::thermo)
    sym=Libdl.dlsym(ct.lib,:stflow_setPressure)
    ccall(sym,
        Cint, (Cint, Cdouble),
        flow_dom, phase.get_P())
    return nothing
end

function _set_grid_default(ct::ct_lib, flow_dom::Int32)
    sym=Libdl.dlsym(ct.lib,:domain_setupGrid)
    grid=collect(range(0,0.1,length=6))
    ccall(sym,
        Cint, (Cint, Csize_t,Ptr{Float64}),
        flow_dom, length(grid),grid)
    return nothing
end

function set_grid(flbase::flame_base, grid::Array{Float64,1})
    sym=Libdl.dlsym(flbase.ct_lib,:domain_setupGrid)
    ccall(sym,
        Cint, (Cint, Csize_t,Ptr{Float64}),
        flbase.doms[2], length(grid),grid)
    return nothing
end

# initialization of outlet
function _outlet_position(ct::ct_lib, outlet::Int32)
    sym=Libdl.dlsym(ct.lib,:domain_setupGrid)
    grid=[0.1]
    ccall(sym,
        Cint, (Cint, Csize_t,Ptr{Float64}),
        outlet, 1 ,grid)
    return nothing
end

# create the oneD base
function flame_base(ct::ct_lib,phase::thermo,
    kin::kin_base,trans::trans_base; itype=2)
    # itype=1 is axisymmetric, 2 is freeflame
    # first create input domain
    doms=Array{Int32,1}(undef, 3)
    sym=Libdl.dlsym(ct.lib,:inlet_new)
    doms[1]=ccall(sym,
        Cint, (),
    )
    inlet_init(ct,doms[1],phase)

    # now create flow domain
    sym=Libdl.dlsym(ct.lib,:stflow_new)
    doms[2]=ccall(sym,
        Cint, (Cint, Cint,Cint, Cint) ,
        phase.base.ind, kin.ind,trans.ind,itype)
    _set_flow_P(ct,doms[2],phase)
    _set_grid_default(ct, doms[2])

    # finally create the outlet
    sym=Libdl.dlsym(ct.lib,:outlet_new)
    doms[3]=ccall(sym,
        Cint, (),
    )
    _outlet_position(ct, doms[3])

    return flame_base(ct.lib, doms)
end

# simulation object

struct Sim1D_base
    ind::Int
    ct_lib::Ptr
    doms::Array{Int32,1}
end

mutable struct Sim1D
    base::Sim1D_base
    solve::Function
end

function _sim1D_solve(sim::Sim1D_base,loglevel, refine_grid)
    sym=Libdl.dlsym(sim.ct_lib,:sim1D_solve)
    solved=ccall(sym,
        Cint, (Cint,Cint, Cint),
        sim.ind, loglevel,refine_grid)
    return solved
end

function Sim1D_base(flbase::flame_base)
    sym=Libdl.dlsym(flbase.ct_lib,:sim1D_new)
    sim_obj=ccall(sym,
        Cint, (Csize_t, Ptr{Cint}) ,
        3,flbase.doms)
    return Sim1D_base(sim_obj,flbase.ct_lib,
        flbase.doms)
end

function Sim1D(sim::Sim1D_base)
    solver(;loglevel=1,refine_grid=true)=_sim1D_solve(sim,loglevel, refine_grid)
    return Sim1D(sim,
        solver)
end