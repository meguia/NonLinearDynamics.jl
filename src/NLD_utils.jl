using Plots
using DifferentialEquations
using ForwardDiff
using IntervalRootFinding
using StaticArrays
using Distances

# utilities

inbox(x,y,xlims,ylims) = (xlims[1]<x<xlims[2]) & (ylims[1]<y<ylims[2])

function realplot!(p1,x,y;plotops...)    
    idx = isreal.(x).*isreal.(y)
    plot!(p1,real(x[idx]),real(y[idx]);plotops...)
end    

function realplot!(x,y;plotops...)    
    idx = isreal.(x).*isreal.(y)
    plot!(real(x[idx]),real(y[idx]);plotops...)
end 

function realplot(x,y;plotops...)    
    p1 = plot()
    realplot!(p1,x,y;plotops...) 
end    


# flows 1D

function dualplot_flow1D(f,x,x0,p,sol;
    size=(900,300),plotops...)

    (ymin,ymax)=extrema(sol.u)    
    (fmin,fmax)=extrema(f.(x,(p,),0.0))
    p1 = plot([x[1],x[end]],[0,0],label="",color=:gray)
    if (x[1]*x[end]<0)
        plot!(p1,[0,0],[fmin,fmax],label="",color=:gray,linestyle=:dash)
    end 
    plot!(p1,x,f.(x,(p,),0.0),label="f(x)",xlabel="x",ylabel="f(x)",color=:blue)
    plot!(p1,sol.u,sol.u*0,label="x",color=:red)
    scatter!(p1,[x0],[0],label="x0",color=:green)
    p2 = plot(sol,label="x(t)",ylabel="x",ylim=(0.9*ymin,1.1*ymax))
    plot(p1,p2,layout=(1,2);size=size,fmt=:png,plotops...)
end    

function flow1D(f::Function,x0::Float64,tmax::Float64,p;
    xlims=[-1.0,1.0],plotops...)

    xrange = xlims[2]-xlims[1]
    x=xlims[1]:xrange/100:xlims[2]
    prob = ODEProblem(f,x0,(0,tmax),p)
    sol = solve(prob)
    dualplot_flow1D(f,x,x0,p,sol;plotops...)
end 

function flow1D(f::Function,x0::Float64,tmax::Float64,p,condition::Function;
    xlims=[-1.0,1.0],plotops...)

    xrange = xlims[2]-xlims[1]
    x=xlims[1]:xrange/100:xlims[2]
    prob = ODEProblem(f,x0,(0,tmax),p)
    condition_terminate(u,t,integrator) = condition(u)
    affect!(integrator) = terminate!(integrator)
    sol = solve(prob,callback=DiscreteCallback(condition_terminate,affect!))
    dualplot_flow1D(f,x,x0,p,sol;plotops...)
end

function flow1D(f::Function,x0::Float64,tmax::Float64,p,tperturb::Float64,Aperturb::Float64,condition::Function;
    xlims=[-1.0,1.0],plotops...)

    xrange = xlims[2]-xlims[1]
    x=xlims[1]:xrange/100:xlims[2]
    prob = ODEProblem(f,x0,(0,tmax),p)
    condition_terminate(u,t,integrator) = condition(u)
    affect1!(integrator) = terminate!(integrator)
    affect2!(integrator) = integrator.u += Aperturb
    cb1 = DiscreteCallback(condition_terminate,affect1!)
    cb2 = PresetTimeCallback([tperturb],affect2!)
    sol = solve(prob,callback=CallbackSet(cb1,cb2))
    dualplot_flow1D(f,x,x0,p,sol;plotops...)
end    

function potential1D(V::Function,x0::Float64,tmax::Float64,p;
    xlims=[-1.0,1.0],size=(900,300),plotops...)

    xrange = xlims[2]-xlims[1]
    x=xlims[1]:xrange/100:xlims[2]
    f(x,p,t) = -ForwardDiff.derivative(x -> V(x,p), x)
    prob = ODEProblem(f,x0,(0,tmax),p)
    sol = solve(prob)
    p1 = plot(x,V.(x,(p,)),label="V(x)",xlabel="x",ylabel="V")
    scatter!(p1,[x0],[V(x0,p)],label="x0")
    plot!(p1,sol,V.(sol.u,(p,)),color=:black,linewidth=2,label="")
    p2 = plot(sol,label="x(t)",ylabel="x",ylim=ylims)
    plot(p1,p2,layout=(1,2);size=size,fmt=:png,plotops...)
end    

function potential1D(V::Function,x0::Float64,tmax::Float64,p,tperturb::Float64,Aperturb::Float64,condition::Function;
    xlims=[-1.0,1.0],size=(900,300),plotops...)

    xrange = xlims[2]-xlims[1]
    x=xlims[1]:xrange/100:xlims[2]
    f(x,p,t) = -ForwardDiff.derivative(x -> V(x,p), x)
    prob = ODEProblem(f,x0,(0,tmax),p)
    condition_terminate(u,t,integrator) = condition(u)
    affect1!(integrator) = terminate!(integrator)
    affect2!(integrator) = integrator.u += Aperturb
    cb1 = DiscreteCallback(condition_terminate,affect1!)
    cb2 = PresetTimeCallback([tperturb],affect2!)
    sol = solve(prob,callback=CallbackSet(cb1,cb2))
    p1 = plot(x,V.(x,(p,)),label="V(x)",xlabel="x",ylabel="V")
    scatter!(p1,[x0],[V(x0,p)],label="x0")
    plot!(p1,sol,V.(sol.u,(p,)),color=:black,linewidth=2,label="")
    p2 = plot(sol,label="x(t)",ylabel="x",ylim=ylims)
    plot(p1,p2,layout=(1,2);size=size,fmt=:png,plotops...)
end    

# 2D

function myquiver!(p1,x, y, u, v; 
    scale=0.3,arrowscale=0.07, color=:black, linealpha=1,plotops...)
    
    function arrow0!(p1, x, y, u, v, as, lc, la)
        nuv = sqrt(u^2 + v^2)
        v1, v2 = [u;v] / nuv,  [-v;u] / nuv
        v4 = (3*v1 + v2)/3.1623 
        v5 = v4 - 2*(v4'*v2)*v2
        v4, v5 = as*nuv*v4, as*nuv*v5
        plot!(p1,[x,x+u], [y,y+v], lc=lc,la=la)
        plot!(p1,[x+u,x+u-v5[1]], [y+v,y+v-v5[2]], lc=lc, la=la)
        plot!(p1,[x+u,x+u-v4[1]], [y+v,y+v-v4[2]], lc=lc, la=la)
    end    
    scatter!(p1,x, y, mc=color, ms=1.5, ratio=1, ma=0.5;legend=false,plotops...)
    for (x,y,u,v) in zip(x,y,u,v)
        arrow0!(p1,x,y,u*scale,v*scale,arrowscale,color,linealpha)
    end
    p1
end

function myquiver(x, y, u, v; 
    scale=0.3,arrowscale=0.07, color=:black, linealpha=1,plotops...)

    p1 = plot(xlabel="x",ylabel="y",fmt=:png)
    myquiver!(p1,x, y, u, v; scale=scale,arrowscale=arrowscale, color=color, linealpha=linealpha,plotops...)
end

function flow2d_grid(f,p,xlims,ylims,npts)
    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    dx = xrange/(npts-1)
    dy = yrange/(npts-1)
    xgrid = xlims[1]:dx:xlims[2]
    ygrid = ylims[1]:dy:ylims[2]
    u = [[x;y] for y in ygrid, x in xgrid]
    Z = u*0
    f.(Z,u,(p,),0)
    dx,dy,xgrid,ygrid,Z
end    

function flow2d_vectorfield!(p1,f::Function,p;
    xlims=[-1.0,1.0],ylims=[-1.0,1.0],npts=21,scale=1.0,plotops...)

    dx,dy,xgrid,ygrid,Z = flow2d_grid(f,p,xlims,ylims,npts)
    u = vec(getindex.(Z,1))
    v = vec(getindex.(Z,2))
    x = vec(repeat(xgrid',length(ygrid)))
    y = repeat(ygrid,length(xgrid))
    maxuv = max(maximum(abs.(u)),maximum(abs.(v)))
    delta = min(dx,dy)
    scale = scale*delta/maxuv
    myquiver!(p1,x,y,u,v;scale=scale,arrowscale=0.5,color=:blue,linealpha=0.5,plotops...)
end    

function flow2d_vectorfield(f::Function,p;
    xlims=[-1.0,1.0],ylims=[-1.0,1.0],npts=21,scale=1.0,size=(600,600),plotops...)

    p1 = plot(xlabel="x",ylabel="y",xlims=xlims,ylims=ylims,fmt=:png,legend=:none,size=size)
    flow2d_vectorfield!(p1,f,p;xlims=xlims,ylims=ylims,npts=npts,scale=scale,plotops...)
end    


function flow2d_vectorfield(f::Function,u0::Vector{Float64},tmax::Float64,p;
    xlims=[-1.0,1.0],ylims=[-1.0,1.0],npts=21,scale=1.0,size=(600,600),plotops...)

    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    p1 = flow2d_vectorfield(f,p;xlims=xlims,ylims=ylims,npts=npts,scale=scale,size=size,plotops...)
    condition(u,t,integrator) = (u[1]*u[1]+u[2]*u[2]) > max(xrange*xrange,yrange*yrange)
    affect!(integrator) = terminate!(integrator)
    sol = solve(ODEProblem(f,u0,(0.0,tmax),p),callback=DiscreteCallback(condition,affect!))
    plot!(p1,sol,vars=(1,2),c=:black,arrow=true,xlims=xlims,ylims=ylims)
end    

function flow2d_vectorfield(f::Function,u0_array::Vector{Vector{Float64}},tmax::Float64,p;
    xlims=[-1.0,1.0],ylims=[-1.0,1.0],npts=21,scale=1.0,size=(600,600),plotops...)

    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    p1 = flow2d_vectorfield(f,p;xlims=xlims,ylims=ylims,npts=npts,scale=scale,size=size,plotops...)
    prob = ODEProblem(f,u0_array[1],(0.0,tmax),p)
    ensamble_prob = EnsembleProblem(prob,prob_func=(prob,i,repeat;u0=u0_array)->(remake(prob,u0=u0[i])))
    condition(u,t,integrator) = (u[1]*u[1]+u[2]*u[2]) > max(xrange*xrange,yrange*yrange)
    affect!(integrator) = terminate!(integrator)
    sol = solve(ensamble_prob,EnsembleThreads(),trajectories=length(u0_array),callback=DiscreteCallback(condition,affect!))
    plot!(p1,sol,vars=(1,2),arrows=true,c=:black,linewidth=0.5,xlims=xlims,ylims=ylims)
end    

function flow2d_nullclines(f::Function,p;
    xlims=[-1.0,1.0],ylims=[-1.0,1.0],npts=21,regions=true,vectorfield=false,size=(700,500),plotops...)
    #=
    =#
    _,_,xgrid,ygrid,Z = flow2d_grid(f,p,xlims,ylims,npts)
    p1 = plot(;xlabel="x",ylabel="y",legend=:none,fmt=:png,size=size,plotops...)
    if regions
        contourf!(p1,xgrid,ygrid,getindex.(Z,1),levels=[-1000,0,1000],alpha=0.8, c=[:red,:green])
        contourf!(p1,xgrid,ygrid,getindex.(Z,2),levels=[-1000,0,1000],alpha=0.4, c=[:blue,:yellow])
        contour!(p1,xgrid,ygrid,getindex.(Z,1),levels=[0],c=:red)
        contour!(p1,xgrid,ygrid,getindex.(Z,2),levels=[0],c=:green)
    else
        contour!(p1,xgrid,ygrid,getindex.(Z,1),levels=[0],c=:red,alpha=0.5)
        contour!(p1,xgrid,ygrid,getindex.(Z,2),levels=[0],c=:green,alpha=0.5)
    end
    if vectorfield
        flow2d_vectorfield!(p1,f,p;xlims=0.9*xlims,ylims=0.9*ylims)
    end  
    p1
end

function flow2d_nullclines(f::Function,u0::Vector{Float64},tmax::Float64,p;
        xlims=[-1.0,1.0],ylims=[-1.0,1.0],npts=21,regions=true,vectorfield=false,size=(700,500),plotops...)
    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    p1 = flow2d_nullclines(f,p;xlims=xlims,ylims=ylims,npts=npts,regions=regions,vectorfield=vectorfield,size=size,plotops...)
    condition(u,t,integrator) = (u[1]*u[1]+u[2]*u[2]) > max(xrange*xrange,yrange*yrange)
    affect!(integrator) = terminate!(integrator)
    sol = solve(ODEProblem(f,u0,(0.0,tmax),p),callback=DiscreteCallback(condition,affect!))
    plot!(p1,sol,vars=(1,2),c=:black,arrow=true,xlims=xlims,ylims=ylims)
end    

function flow2d_nullclines(f::Function,u0_array::Vector{Vector{Float64}},tmax::Float64,p;
    xlims=[-1.0,1.0],ylims=[-1.0,1.0],npts=21,regions=true,vectorfield=false,size=(700,500),plotops...)
    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    #u0_arr = vec([[xlims[1]+i*xrange/Ngrid,ylims[1]+j*yrange/Ngrid] for i=0:Ngrid, j=0:Ngrid])
    p1 = flow2d_nullclines(f,p;xlims=xlims,ylims=ylims,npts=npts,regions=regions,vectorfield=vectorfield,size=size,plotops...)
    prob = ODEProblem(f,u0_array[1],(0.0,tmax),p)
    ensamble_prob = EnsembleProblem(prob,prob_func=(prob,i,repeat;u0=u0_array)->(remake(prob,u0=u0[i])))
    condition(u,t,integrator) = (u[1]*u[1]+u[2]*u[2]) > max(xrange*xrange,yrange*yrange)
    affect!(integrator) = terminate!(integrator)
    sol = solve(ensamble_prob,EnsembleThreads(),trajectories=length(u0_array),callback=DiscreteCallback(condition,affect!))
    plot!(p1,sol,vars=(1,2),arrows=true,c=:black,linewidth=0.5,xlims=xlims,ylims=ylims)
end    

function flow2d_animated(f::Function,p,N::Int64,dt::Float64;
    Ngrid=10,fps=15,xlims=[-1.0,1.0],ylims=[-1.0,1.0],size=(400,400),nullclines=false,fname="")

    if isempty(fname)
        fname = string(Symbol(f))
        if fname[end] == '!'
            fname = chop(fname)
        end  
    end      
    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    u0_arr = vec([[xlims[1]+i*xrange/Ngrid,ylims[1]+j*yrange/Ngrid] for i=0:Ngrid, j=0:Ngrid])
    prob = ODEProblem(f,u0_arr[1],(0.0,N*dt),p)
    ensamble_prob = EnsembleProblem(prob,prob_func=(prob,i,repeat;u0=u0_arr)->(remake(prob,u0=u0[i])))
    sol = solve(ensamble_prob,EnsembleThreads(),trajectories=length(u0_arr))
    x = reduce(hcat,[sol[n](0:dt:N*dt,idxs=1).u for n=1:length(sol)])
    y = reduce(hcat,[sol[n](0:dt:N*dt,idxs=2).u for n=1:length(sol)])
    if (nullclines)
        p1 = flow2d_nullclines(f,p;xlims=xlims,ylims=ylims,size=size) 
    else    
        p1 = plot([x[1:2,:]],[y[1:2,:]],color=:black,xlims=xlims,ylims=ylims,legend=false,size=size)
    end  
    # make animation 
    anim = @animate for n=3:N
        plot!(p1,x[n-1:n,:],y[n-1:n,:],color=:black)
    end
    print("Saving animation ...")
    gif(anim, fname * "_fps_" * string(fps) * ".gif", fps=fps)
end       

function classification_linear(A::Matrix{Float64};
    Ngrid=5,tmax=2.0,xlims=[-1.2,1.2],ylims=[-1.2,1.2],circular=true)

    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    if circular
        Npts = Ngrid*Ngrid
        u0_arr = u0_arr = vec([[cos(i*2pi/Npts);sin(i*2pi/Npts)] for i=1:Npts])
    else
        u0_arr = vec([[xlims[1]+i*xrange/Ngrid,ylims[1]+j*yrange/Ngrid] for i=0:Ngrid, j=0:Ngrid])
    end    
    prob = ODEProblem((u,p,t) -> A*u, [0;0], (0,tmax))
    ensamble_prob = EnsembleProblem(prob,prob_func=(prob,i,repeat;u0=u0_arr)->(remake(prob,u0=u0[i])))
    sol = solve(ensamble_prob,EnsembleThreads(),trajectories=length(u0_arr))
    p2 = plot(sol,vars=(1,2),xlims=xlims,ylims=ylims,arrows=true,c=:black,linewidth=0.5,fmt=:png)
    maxtr = 2.0
    maxdet = 1.3
    trv = -maxtr:maxtr/30:maxtr
    trdet = abs2.(trv)/4
    p1 = plot(trv,trdet,xaxis=("Tr(A)",(-maxtr,maxtr)),yaxis=("Det(A)",(-maxdet,maxdet)),legend=false,fmt = :png)
    plot!(p1,[-maxtr,maxtr],[0,0])
    plot!(p1,[0,0],[0,maxdet])
    txt=["Saddle","Stable\nNode","Stable\nFocus","Unstable\nFocus","Unstable\nNode"]
    annotate!([(0,-0.5*maxdet,(txt[1],12,:gray))])
    annotate!([(-0.8*maxtr,0.2*maxdet,(txt[2],10,:gray)),(-0.5*maxtr,0.7*maxdet,(txt[3],10,:gray))])
    annotate!([(0.8*maxtr,0.2*maxdet,(txt[5],10,:gray)),(0.5*maxtr,0.7*maxdet,(txt[4],10,:gray))])
    atxt = string(A[1,1]) * " " * string(A[1,2]) * "\n" * string(A[2,1]) * " " * string(A[2,2])
    scatter!(p1,[tr(A)],[det(A)],c=:red)
    annotate!([(tr(A),det(A)-0.1*maxdet,(atxt,9,:red,:center))])
    plot(p2,p1,layout=(1,2),size=(900,400))
end

function flow2d_manifolds!(p1,f,f_jac,u0_array,p;
    tmax=30,delta=0.001,repulsor=false,xlims=[-1.0,1.0],ylims=[-1.0,1.0],size=(700,500),plotops...)

    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    condition(u,t,integrator) = (u[1]*u[1]+u[2]*u[2]) > max(xrange*xrange,yrange*yrange)
    affect!(integrator) = terminate!(integrator)
    for u0 in u0_array
        A = f_jac(u0,p)
        if det(A)<0
            av = eigen(A)
            for n=1:2
                if real(av.values[n])>0
                    u1 = u0+delta*av.vectors[:,n]
                    sol = solve(ODEProblem(f,u1,(0.0,tmax),p),callback=DiscreteCallback(condition,affect!))
                    plot!(p1,sol,vars=(1,2),c=:red,label="Wu")
                    u1 = u0-delta*av.vectors[:,n]
                    sol = solve(ODEProblem(f,u1,(0.0,tmax),p),callback=DiscreteCallback(condition,affect!))
                    plot!(p1,sol,vars=(1,2),c=:red,label="")
                else
                    u1 = u0+delta*av.vectors[:,n]
                    sol = solve(ODEProblem(f,u1,(0.0,-tmax),p),callback=DiscreteCallback(condition,affect!))
                    plot!(p1,sol,vars=(1,2),c=:blue,label="Ws")
                    u1 = u0-delta*av.vectors[:,n]
                    sol = solve(ODEProblem(f,u1,(0.0,-tmax),p),callback=DiscreteCallback(condition,affect!))
                    plot!(p1,sol,vars=(1,2),c=:blue,label="")
                end  
            end
        else
            if tr(A)<0
                scatter!(p1,[u0[1]],[u0[2]],c=:black,label="A")
            else
                scatter!(p1,[u0[1]],[u0[2]],c=:red,label="R")
                if repulsor
                    u1 = u0.+[delta;delta]
                    sol = solve(ODEProblem(f,u1,(0.0,3*tmax),p),callback=DiscreteCallback(condition,affect!))
                    plot!(p1,sol,vars=(1,2),c=:purple,alpha=0.5,label="")
                end    
            end    
        end    
    end
    plot(p1;legend=false,xlims=xlims,ylims=ylims,size=size,plotops...)
end

function flow2d_manifolds(f::Function,f_jac::Function,u0_array::Vector{Vector{Float64}},p::Vector{Float64};
    tmax=30,delta=0.001,repulsor=false,xlims=[-1.0,1.0],ylims=[-1.0,1.0],size=(700,500),plotops...)

    p1=plot(xlabel="x",ylabel="y",fmt=:png)
    flow2d_manifolds!(p1,f,f_jac,u0_array,p;tmax=tmax,delta=delta,repulsor=repulsor,xlims=xlims,ylims=ylims,size=size,plotops...)
end

function phase_portrait(f::Function,p::Vector{Float64};
    tmax=50,delta=0.001,xlims=[-1.0,1.0],ylims=[-1.0,1.0],size=(700,500),plotops...)
    # Find fixedpoints in interval 
    function fsv((x,y))
        du = f(similar([x,y]),[x,y],p,0)
        return SVector(du[1],du[2])
    end    
    X = xlims[1]..xlims[2]
    Y = ylims[1]..ylims[2]
    rts = roots(fsv, X × Y)
    u0_arr=[[mid(rt.interval[1]);mid(rt.interval[2])] for rt in rts]
    p1 = flow2d_nullclines(f,p;xlims=xlims,ylims=ylims,npts=50,regions=false)
    f_jac(u0,p) = ForwardDiff.jacobian(x -> f(similar(x),x,p,0), u0)
    flow2d_manifolds!(p1,f,f_jac,u0_arr,p;tmax=tmax,delta=delta,repulsor=true,xlims=xlims,ylims=ylims,size=size,plotops...)
end    
    
function attractor_basin(f::Function,p::Vector{Float64},attractors::Vector{Vector{Float64}},maxdist::Float64;
    delta=0.1,tmax=1000.0,xlims=(-1.0,1.0),ylims=(-1.0,1.0),size=(900,600),plotops...)

    natt = length(attractors)
    clist = [:black,:red,:blue,:yellow,:green,:purple,:cyan,:orange]
    if (natt>7) 
        error("maximum number of attractors is 7")
    end    
    distance2D(x1,x2) = sqrt((x1[1]-x2[1])^2+(x1[2]-x2[2])^2)
    x = xlims[1]:delta:xlims[2]
    y = ylims[1]:delta:ylims[2]
    nx = length(x)
    ny = length(y)
    u0_list=[[xi,yi,0.0] for xi in x for yi in y]
    prob = ODEProblem(f,u0_list[1],(0,tmax),p)
    ensamble_prob = EnsembleProblem(prob,prob_func=(prob,i,repeat;u0=u0_list)->(remake(prob,u0=u0_list[i])))
    sol = solve(ensamble_prob,EnsembleThreads(),trajectories=length(u0_list),save_everystep=false)
    println("solved")
    M = zeros(Int8,nx,ny)
    for (n,pt) in enumerate(sol)
        for (m,at) in enumerate(attractors)
            if (distance2D(pt.u[end][1:2],at)<maxdist)
                M[(n-1)÷ny+1,mod(n-1,ny)+1] = m
            end
        end    
    end
    M[1,1]=0
    println("plotting...")
    contourf(x,y,M',c=clist[1:natt+1],linewidth=0,legend=false,xlabel="x",ylabel="y";size=size,fmt=:png,plotops...)    
end    

function flow2d_forced(f::Function,u0::Vector{Float64},p::Vector{Float64},period::Float64; 
    tcycles=0, npts=100,ncycles=10,size=(900,400),xlims=false,ylims=false,plotops...)
    # Assume that the 3rd variable (z) is periodic and plot x,y as a function of z modulo period
    # skip trans period of the forcing
    tmax = min(ncycles*period,100000)
    sol = solve(ODEProblem(f,u0,(0.0,tmax),p))
    ncycles = Int(floor(tmax/period))
    p1 = plot(xlabel="x",ylabel="y",zlabel="ϕ/2π",legend=false,fmt = :png)
    p2 = plot(xlabel="x",ylabel="y",legend=false,fmt = :png)
    for j = tcycles+1:ncycles
        ts = range((j-1)*period,stop=j*period,length=npts+1)
        plot!(p1,sol(ts,idxs=1),sol(ts,idxs=2),(0:npts)/npts,color=j)
        scatter!(p1,[sol(ts[end],idxs=1)],[sol(ts[end],idxs=2)],[1],markersize=2,color=j+1)
        scatter!(p1,[sol(ts[1],idxs=1)],[sol(ts[1],idxs=2)],[0],markersize=2,color=j,framestyle=:box)
        plot!(p2,sol(ts,idxs=1),sol(ts,idxs=2),color=j,linealpha=0.3)
        scatter!(p2,[sol(ts[1],idxs=1)],[sol(ts[1],idxs=2)],markersize=2,color=j)
    end
    scatter!(p2,[sol(ncycles*period,idxs=1)],[sol(ncycles*period,idxs=2)],markersize=2,color=ncycles+1)    
    if xlims isa Tuple
        xlims!(p1,xlims)
        xlims!(p2,xlims)
    end    
    if ylims isa Tuple
        ylims!(p1,ylims)
        ylims!(p2,ylims)
    end    
    plot(p1,p2,layout=(1,2);size=size,plotops...)
end        

function poincare_forced!(p1,f, u0, p, period; 
    tcycles=0, ncycles=10,msize=0.5,col=:blue,xlims=false,ylims=false,plotops...)
    # Assume that the 3rd variable (z) is periodic and plot x,y as a function of z modulo period
    # skip trans period of the forcing
    trans = solve(ODEProblem(f,u0,(0.0,tcycles*period),p))
    u0=trans.u[end]
    if period>1000
        tmax=min(ncycles*period,100000)
    else
        tmax=ncycles*period
    end    
    tsave = period:period:ncycles*period
    sol = solve(ODEProblem(f,u0,(0.0,tmax),p),saveat=tsave)
    scatter!(p1,sol,vars=(1,2);legend=false,markersize=msize,color=col,markerstrokewidth=0,plotops...)
    if (xlims)
        (xmin,xmax) = extrema(getindex.(vcat(trans.u,sol.u),1))
        xrange=xmax-xmin
        xlims!(p1,xmin-0.1*xrange,xmax+0.1*xrange)
    end
    if (ylims)    
        (ymin,ymax) = extrema(getindex.(vcat(trans.u,sol.u),2))
        yrange=ymax-ymin
        ylims!(p1,ymin-0.1*yrange,ymax+0.1*yrange)
    end    
    p1
end   

function poincare_forced(f::Function,u0::Vector{Float64},p::Vector{Float64},period::Float64;
    tcycles=0, ncycles=10,size=(500,500),msize=0.5,col=:blue,xlims=false,ylims=false,plotops...)

    p1=plot(xlabel="x",ylabel="y",size=size,fmt=:png)
    poincare_forced!(p1,f, u0, p, period; tcycles=tcycles,ncycles=ncycles,msize=msize,col=col,xlims=xlims,ylims=ylims,plotops...)
end    

function poincare_forced_zoom(f::Function,u0::Vector{Float64},p::Vector{Float64},period::Float64;
    npts=1000,maxiter=1000,size=(600,600),xlims=[-1,1],ylims=[-1,1],plotops...)
    # Assume that the 3rd variable (z) is periodic and plot x,y as a function of z modulo period
    # skip trans period of the forcing
    tcycles = 30
    trans = solve(ODEProblem(f,u0,(0.0,tcycles*period),p))
    u0=trans.u[end]
    p1 = plot(;xlabel="x",ylabel="y",legend=false,size=size,fmt=:png,plotops...)
    ncycles = 1000
    tmax=ncycles*period
    kpts = 0
    kiter = 0
    while (kpts<npts) & (kiter<maxiter)
        tsave = period:period:ncycles*period
        sol = solve(ODEProblem(f,u0,(0.0,tmax),p),saveat=tsave)
        xlist = getindex.(sol.u,1)
        ylist = getindex.(sol.u,2)
        idx=inbox.(xlist,ylist,Ref(xlims),Ref(ylims))
        scatter!(getindex.(sol.u[idx],1),getindex.(sol.u[idx],2),color=:blue,markersize=0.5,xlims=xlims,ylims=ylims)
        u0 = sol.u[end]
        kpts += sum(idx)
        kiter += 1
    end
    println(kpts)
    println(kiter)
    p1
end     

function recurrence_plot(f::Function,u0::Vector{Float64},p::Vector{Float64},period::Float64;
    dd=0.002,steps=10,tcycles=0,npts=300,ncycles=10,size=(900,450),plotops...)

    trans = solve(ODEProblem(f,u0,(0.0,tcycles*period),p))
    u0=trans.u[end]
    if period>1000
        tmax=min(ncycles*period,10000)
    else
        tmax=ncycles*period
    end  
    sol = solve(ODEProblem(f,u0,(0,tmax),p));
    ts = range(0,tmax,length=npts)
    xy=hcat(sol(ts,idxs=1).u,sol(ts,idxs=2).u,mod.(sol(ts,idxs=2).u,period))
    dist = pairwise(Euclidean(),xy,dims=1)
    dst = floor.(dist/dd)/steps
    dst[dst.>steps] .= steps
    p1 = plot(sol,vars=(1,2),xlabel="x",ylabel="y")
    p2 = heatmap(ts,ts,dst,c=cgrad([:blue,:white]),legend=false,size=size)
    plot(p1,p2,layout=(1,2);size=size,fmt=:png,plotops...)
end    

function saddle_orbit2D(f::Function,u0::Vector{Float64},p::Vector{Float64},period::Float64;
    λ=0.001,maxiter=10000, disttol=1e-9, inftol=10)
    # funcion que mapea un ciclo (extrender a N ciclos)
    function map_period(f,u0,p,period)
        sol = solve(ODEProblem(f,[u0[1],u0[2],0.0],(0,period),p),Tsit5(),reltol=1e-12,save_everystep=false)
        return sol.u[end][1:2]
    end   
    # Funcion que aplica el metodo Schmelcher & Diakonos
    function Sk(f,prus,p,period,Λ)
        us = map_period(f,prus,p,period)
        return prus, prus+Λ*(us-prus)
    end
    # base de matrices del metodo para 2D
    Λ_base = λ*[[1 0; 0 1],[1 0; 0 -1],[-1 0; 0 1]]
    converged = false
    fp = u0
    # iteracion sobre los tres Λ base
    for Λ in Λ_base
        prus = u0
        # iteracion del metodo
        for i in 1:maxiter
            prus, us = Sk(f,prus, p,period,Λ)
            norm(us) > inftol && break
            if norm(prus - us) < disttol
                fp = us
                converged = true
                break
            end
            prus = us
        end
    end    
    return fp,converged
end


function saddle_manifolds_forced(f::Function,f_jac::Function,us::Vector{Float64},p::Vector{Float64},period::Float64;
    ncycles=[10,3],npts=300,delta=0.01,xlims=(-1.0,1.0),ylims=(-1.0,1.0),size=(900,600),plotops...)
    # Asume que us es la orbita saddle en el plano 
    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    condition(u,t,integrator) = (u[1]*u[1]+u[2]*u[2]) > max(xrange*xrange,yrange*yrange)
    affect!(integrator) = terminate!(integrator)
    t0_array = 0:period/npts:period
    sol = solve(ODEProblem(f,[us[1],us[2],0.0],(0,period),p),Tsit5(),reltol=1e-12,saveat=t0_array)
    u0_array=sol.u
    p1 = plot(;xlabel="x",ylabel="y",legend=false,size=size,fmt=:png,plotops...)
    for (u0,t0) in zip(u0_array,t0_array)
        A = f_jac(u0,p)
        av = eigen(A)
        for n=1:3
            u1 = u0+delta*av.vectors[:,n]
            u2 = u0-delta*av.vectors[:,n]
            if real(av.values[n])>0
                tmax=ncycles[1]*period
                tsave = period:period:ncycles[1]*period
                sol = solve(ODEProblem(f,u1,(t0,tmax),p),callback=DiscreteCallback(condition,affect!),saveat=tsave)
                scatter!(p1,sol,vars=(1,2),c=:red,markerstrokewidth=0,markersize=1)
                sol = solve(ODEProblem(f,u2,(t0,tmax),p),callback=DiscreteCallback(condition,affect!),saveat=tsave)
                scatter!(p1,sol,vars=(1,2),c=:red,markerstrokewidth=0,markersize=1)
            elseif real(av.values[n])<0
                tmax=ncycles[2]*period
                tsave = 0:-period:-ncycles[2]*period
                sol = solve(ODEProblem(f,u1,(t0,-tmax),p),callback=DiscreteCallback(condition,affect!),saveat=tsave)
                scatter!(p1,sol,vars=(1,2),c=:blue,markerstrokewidth=0,markersize=1)
                sol = solve(ODEProblem(f,u2,(t0,-tmax),p),callback=DiscreteCallback(condition,affect!),saveat=tsave)
                scatter!(p1,sol,vars=(1,2),c=:blue,markerstrokewidth=0,markersize=1)
            end 
        end
    end    
    xlims!(p1,xlims)
    ylims!(p1,ylims)
end    


function butterfly(f::Function,u0::Vector{Float64},p::Vector{Float64}; 
    delta=1e-12,size=(900,400),xlims=false,ylims=false,plotops...)
    
    u1=u0
    u1[1]=u[1]+delta
    sol0 = solve(ODEProblem(f,u0,(0.0,tmax),p))
    sol1 = solve(ODEProblem(f,u1,(0.0,tmax),p))
    if length(u0) == 2
        p1 = plot(sol0,vars=(1,2),c=:red)
        plot!(p1,sol1,vars=(1,2),c=:blue)
    elseif length(u0) == 3        
        p1 = plot(sol0,vars=(1,2,3),c=:red)
        plot!(p1,sol1,vars=(1,2,3),c=:blue)  
    end    
    dd = 
    plot(p1,p2,layout=(1,2);size=size,plotops...)
end

function butterfly(f::Function,u0::Vector{Float64},p::Vector{Float64},tmax::Float64; 
    dim=3,dt=0.001,delta=1e-12,size=(900,400),xlims=false,ylims=false,plotops...)    
    
    u1=copy(u0)
    u1[1]=u1[1]+delta
    ts = range(0.0,stop=tmax,step=dt)
    sol0 = solve(ODEProblem(f,u0,(0.0,tmax),p))
    sol1 = solve(ODEProblem(f,u1,(0.0,tmax),p))
    if dim == 2
        p1 = plot(sol0,vars=(1,2),c=:red,xlabel="x",ylabel="y",label="u0")
        plot!(p1,sol1,vars=(1,2),c=:blue,label="u1")
    elseif dim == 3
        p1 = plot(sol0,vars=(1,2,3),c=:red,xlabel="x",ylabel="y",label="u0")
        plot!(p1,sol1,vars=(1,2,3),c=:blue,label="u1")  
    end    
    x0 = sol0(ts,idxs=1)
    x1 = sol1(ts,idxs=1)
    p2 = plot(ts,log10.(norm.(x0-x1)),xlabel="t",ylabel="log10(d)",label="")
    plot(p1,p2,layout=(1,2);size=size,plotops...)
end        