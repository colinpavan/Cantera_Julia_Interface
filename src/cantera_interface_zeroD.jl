struct reactor
    ind::Int
    gas::gas
end

function reactor(G::gas,type::String="ConstPressureReactor")
    # load in the file
    sym=Libdl.dlsym(lib,:reactor_new)
    r_obj=ccall(sym,
        Cint, (Cstring,) ,
        type)
    ccall(Libdl.dlsym(lib,:reactor_setThermoMgr),
        Cint,(Cint,Cint),r_obj,G.phase.ind)
    ccall(Libdl.dlsym(lib,:reactor_setKineticsMgr),
        Cint,(Cint,Cint),r_obj,G.kin.ind)
    return reactor(r_obj,G)
end

function delete_reactor(r::reactor)
    sym=Libdl.dlsym(lib,:reactor_del)
    ccall(sym,Cint,(Cint,),r.ind)
    r=nothing
    return nothing
end

function reactor_sim_full(r::reactor,tmax::Float64)
    ind=ccall(Libdl.dlsym(lib,:reactornet_new),
        Cint,(),)
    ccall(Libdl.dlsym(lib,:reactornet_addreactor),
        Cint,(Cint,Cint),ind,r.ind)
    ccall(Libdl.dlsym(lib,:reactornet_setInitialTime),
        Cint,(Cint,Cdouble),ind,0.0)
    ccall(Libdl.dlsym(lib,:reactornet_advance),
        Cint,(Cint,Cdouble),ind,tmax)
    ccall(Libdl.dlsym(lib,:reactornet_del),
        Cint,(Cint,),ind)
    return nothing
end

struct reactor_net
    ind::Int
    r::reactor
end

function reactor_sim_initialize(r::reactor)
    ind=ccall(Libdl.dlsym(lib,:reactornet_new),
        Cint,(),)
    ccall(Libdl.dlsym(lib,:reactornet_addreactor),
        Cint,(Cint,Cint),ind,r.ind)
    return reactor_net(ind,r)
end

function reactor_sim_advance(RN::reactor_net,tmax::Float64)
    # for some reason setting the time and then
    # immediately advancing is slow
    # better to just continue from last point
    #ccall(Libdl.dlsym(lib,:reactornet_setInitialTime),
    #    Cint,(Cint,Cdouble),RN.ind,0.0)
    t0=ccall(Libdl.dlsym(lib,:reactornet_time),
        Cint,(Cint,),RN.ind)
    ccall(Libdl.dlsym(lib,:reactornet_advance),
        Cint,(Cint,Cdouble),RN.ind,t0+tmax)
end

function reactor_sim_advance2(RN::reactor_net,tmax::Float64)
    ccall(Libdl.dlsym(lib,:reactornet_advance),
        Cint,(Cint,Cdouble),RN.ind,tmax)
end


function reactor_sim_delete(RN::reactor_net)
    ccall(Libdl.dlsym(lib,:reactornet_del),
        Cint,(Cint,),RN.ind)
    return nothing
end