#11/19/20: Upload to Github

#Community consumption model ccm
#Organic matter consumption by microbial community
#ccm.jl by EJZ. Begin Apr 9, 2020

#5/14/20: update to include penalty as parameter input,
    #fix dCexdt for output as Ctot,
    #have flex printing of CM
    #copy dIC and bIC to preserve for output
#This is the SWITCH code, with the "regular" one commented out

using SparseArrays, LinearAlgebra
using Statistics
using NCDatasets
using Dates
using Printf

#make sum drop dimension automatically:
sumd(x,dims) = dropdims(sum(x,dims=dims),dims=dims)

#set subnormals to zero to avoid slowdown
#this way, don't have to cap b and d at a lower limit
set_zero_subnormals(true)

#make a Composite Type using "struct" to organize parameters
struct Prms
    fsave::String
    T::Int64
    Pt::Float64
    nd::Int64
    nb::Int64
    CM::Array{Bool,2}
    pen::Array{Float64,1}
    PA::Array{Float64,1}
    Cw::Array{Float64,1}
    D::Float64
    vmaxB::Array{Float64,1}
    km::Array{Float64,1}
    y::Array{Float64,1}
    mq::Array{Float64,1}
    mlin::Array{Float64,1}
    dIC::Array{Float64,1}
    bIC::Array{Float64,1}
    nrec::Float64
end

println("Loading function: runccm")
function runccm(prms)

    @printf("nd = %5.0f, nb = %5.0f, Pt =  %2.0f, T = %5.0f \n", prms.nd, prms.nb, prms.Pt, prms.T)
    open("jlog.txt","w") do f
        write(f,@sprintf("nd = %5.0f, nb = %5.0f, Pt =  %2.0f, T = %5.0f \n", prms.nd, prms.nb, prms.Pt, prms.T))
    end

    #Starting time and saving location
    tst = now()
    fsaven = string(prms.fsave,"_", Dates.format(tst, "yyyymmdd"), ".nc")
    #fsaven = string(prms.fsave,".nc")
    if isfile(fsaven) #if i do more than 1 in 1 day
        fsaven = string(prms.fsave,"_",Dates.format(tst, "yyyymmdd_HHMM"), ".nc")
    end
    println("Starting time: ", tst)
    println("Will be saved as: ", fsaven)

    #################################################################################
    #Prep for integration

    Cs = sparse(prms.CM)
    (II, JJ, _) = findnz(Cs) #this gives the list of all "on" indices in the matrix 

    println("CM[1:10,1:10]")
    for i=1:min(prms.nd,10)
        println(prms.CM[i,1:min(prms.nb,10)])
    end

    dt = 1e-3
    nt = Int(prms.T/dt)
    println("nt = ",nt)
    trec = nt÷prms.nrec #frequency of recording

    #set up empty matrices
    nrec1 = Int(prms.nrec+1) #bc i added time 0
    d = Array{Float64,2}(undef, prms.nd, nrec1)
    b = Array{Float64,2}(undef, prms.nb, nrec1)
    Cwall = Array{Float64,2}(undef, prms.nd, nrec1)
    PAall = Array{Float64,2}(undef, prms.nb, nrec1)
    Ctot = Array{Float64,1}(undef, nrec1) 
    timet = Array{Float64,1}(undef, nrec1) 
    vall = Array{Float64,3}(undef, prms.nd, prms.nb, nrec1)

    #################################################################################
    #Initial conditions IC

    dtemp = copy(prms.dIC) 

    btemp = copy(prms.bIC) 

    Cextemp = 0.

    #write IC as timet = 0
    timet[1] = 0
    d[:,1] = dtemp
    b[:,1] = btemp
    Cwall[:,1] .= NaN
    #Cwall[:,1]=fill!(zeros(prms.nd),NaN)
    PAall[:,1] .= NaN
    #PAall[:,1]=fill!(zeros(prms.nb),NaN)
    Ctot[1] = sum(dtemp)+sum(btemp)+Cextemp
    vall[:,:,1] .= NaN
    #################################################################################
    #Time loop

    dbdt1 = Array{Float64,1}(undef, prms.nb)
    dbdt2 = Array{Float64,1}(undef, prms.nb)
    dbdt3 = Array{Float64,1}(undef, prms.nb)
    dbdt4 = Array{Float64,1}(undef, prms.nb)
    dddt1 = Array{Float64,1}(undef, prms.nd)
    dddt2 = Array{Float64,1}(undef, prms.nd)
    dddt3 = Array{Float64,1}(undef, prms.nd)
    dddt4 = Array{Float64,1}(undef, prms.nd)
    b1 = Array{Float64,1}(undef, prms.nb)
    b2 = Array{Float64,1}(undef, prms.nb)
    b3 = Array{Float64,1}(undef, prms.nb)
    d1 = Array{Float64,1}(undef, prms.nd)
    d2 = Array{Float64,1}(undef, prms.nd)
    d3 = Array{Float64,1}(undef, prms.nd)
    #for switch:
    v = Array{Float64,2}(undef, prms.nd, prms.nb)
            
    @time for t = 1:nt #nt

        #@time begin #output at each timestep
        
        #print("v:")
        #print(v)
            
        #randomized supply according to Cw set above:
        Cwtemp = [rand() <= prms.Cw[i] for i = 1:prms.nd]
        Pd = prms.Pt.*Cwtemp/prms.nd

        #presence-absence of pops:
        PAtemp = [rand() <= prms.PA[i] for i = 1:prms.nb]

        dCexdt1 = ecof!(btemp, dtemp, Pd, PAtemp, II, JJ, prms, dbdt1, dddt1, v)

        b1 .= btemp .+ dt/2 .* dbdt1
        d1 .= dtemp .+ dt/2 .* dddt1

        dCexdt2 = ecof!(b1, d1, Pd, PAtemp, II, JJ, prms, dbdt2, dddt2, v)

        b2 .= btemp .+ dt/2 .* dbdt2
        d2 .= dtemp .+ dt/2 .* dddt2

        dCexdt3 = ecof!(b2, d2, Pd, PAtemp, II, JJ, prms, dbdt3, dddt3, v)

        b3 .= btemp .+ dt .* dbdt3
        d3 .= dtemp .+ dt .* dddt3

        dCexdt4 = ecof!(b3, d3, Pd, PAtemp, II, JJ, prms, dbdt4, dddt4, v)
       
        btemp .+= (dbdt1 .+ 2 .* dbdt2 .+ 2 .* dbdt3 .+ dbdt4).*(dt/6)
        dtemp .+= (dddt1 .+ 2 .* dddt2 .+ 2 .* dddt3 .+ dddt4).*(dt/6)
        Cextemp += dt/6 * (dCexdt1 + 2dCexdt2 + 2dCexdt3 + dCexdt4)

        dtemp[dtemp .< 0] .= 0
            
        #print every 100 timesteps:
        if mod(t, 100)==0
            @printf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", t*dt, prms.T, t*dt/prms.T*100, now())
        end
        
        if mod(t, trec)==0

            @printf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", t*dt, prms.T, t*dt/prms.T*100, now())
            open("jlog.txt","a") do f
                write(f,@sprintf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", t*dt, prms.T, t*dt/prms.T*100, now()))
            end
            j = Int(t÷trec + 1)
            timet[j]=t.*dt
            b[:,j]=btemp
            d[:,j]=dtemp
            Cwall[:,j]=Cwtemp
            PAall[:,j]=PAtemp
            Ctot[j] = sum(btemp)+sum(dtemp)+Cextemp
            #for switch:
            vall[:,:,j]=v
        end
    
        #end #end recording timestep time

    end #end time loop

    # #calculate uptake and v just for last timepoint: use vmaxB, not vmaxtemp
    # #NOT for switching
    vB = zeros(prms.nd,prms.nb) 
    v = zeros(prms.nd,prms.nb) 
    uptake = zeros(prms.nd,prms.nb) 
    # @inbounds for n = axes(II, 1)
    #     vB[II[n],JJ[n]] = prms.vmaxB[II[n]] * prms.pen[JJ[n]] * dtemp[II[n]]/(dtemp[II[n]]+prms.km[II[n]])
    #     v[II[n],JJ[n]] = prms.vmaxB[II[n]]* prms.pen[JJ[n]] * dtemp[II[n]]/(dtemp[II[n]]+prms.km[II[n]]) * prms.PA[JJ[n]]
    #     uptake[II[n],JJ[n]] = v[II[n],JJ[n]] * btemp[JJ[n]]
    # end
    @inbounds for m = 1:prms.nb
        maxu,uind = findmax(prms.CM[:,m] .* prms.vmaxB .* prms.pen[m] .* dtemp./(dtemp+prms.km) .* prms.y * PA[m])
        v[uind,m] = maxu/prms.y[uind]
        uptake[uind,m] = v[uind,m]*btemp[m] 
    end 

    tfn = now() #finish time
    
    #Write files
    savetoNC(fsaven, d, b, vB, v, uptake, timet, Ctot, Cwall, PAall, tst, tfn, prms)
    
    #return d, b, timet, Ctot, v, prms #, dCexdt1, dddt1
    return d, b, timet, Ctot, vall, prms #, dCexdt1, dddt1

end

#################################################################################
#Ecosystem equations and integration

println("Loading function: ecof")

function ecof!(b, d, Pd, PAtemp, II, JJ, prms, dbdt, dddt, v)

    fill!(v,0)

    d[d .< 0] .= 0
    
    dddt .= Pd - d.*prms.D
    
    dbdt .= -b.*(prms.mq.*b .+ prms.mlin .+ prms.D)
    
    dCexdt = sum(b.*(prms.mq.*b .+ prms.mlin))

    #uptake and respiration:
    #@inbounds for n = axes(II, 1)
#    for n = axes(II, 1)
        #n from 1 to n of interactions,
        #II[n] values are from 1 to nd (IF there is an interaction at i=nd!)
        #JJ[n] values are from 1 to nb
        
        #uptake = b[JJ[n]] * prms.vmaxB[II[n]] * prms.pen[JJ[n]] * d[II[n]]/(d[II[n]]+prms.km[II[n]]) * PAtemp[JJ[n]]
   
        #dddt[II[n]] += - uptake 
      
        #dbdt[JJ[n]] += uptake*prms.y[II[n]] 
       
        #dCexdt += uptake*(1 - prms.y[II[n]])
        #println(dCexdt)
    #end

    #switching:
    @inbounds for m = 1:prms.nb
        
        #find maximum growth rate and i index on all possible substrates:
        #right now, this DOESN'T scale with # of interactions!
        maxu,uind = findmax(prms.CM[:,m] .* prms.vmaxB .* prms.pen[m] .* d./(d+prms.km) .* prms.y * PAtemp[m])
        uptake = maxu/prms.y[uind]*b[m]

        #try:
        #every 10 timesteps (a few hours for the cell)
        #take out the constant factors and collapse
        #do manual find max: go through II and JJ and find max myself: might have to make new structure where new JJ is for each column

        #doing both of these is slower:
        #maxu = maximum(CM[:,m] .* prms.vmaxB .* d./(d+prms.km) .* prms.y * PAtemp[m])
        #uind = argmax(CM[:,m] .* prms.vmaxB .* d./(d+prms.km) .* prms.y * PAtemp[m])

        #this is about the same, slightly slower:
        #ups = CM[:,m] .* prms.vmaxB .* d./(d+prms.km) .* prms.y * PAtemp[m]
        #maxu = maximum(ups)
        #uind = argmax(ups)
       
        #v[uind,m] = maxu/prms.y[uind] 
        v[uind,m] = uptake 
        
        dddt[uind] += - uptake 
      
        dbdt[m] += uptake*prms.y[uind]

        #dCexdt += maxu*(1. /prms.y[uind] - 1.)
        dCexdt += uptake*(1 - prms.y[uind]) 

    end
      
    return dCexdt

end

#################################################################################

println("Loading function: savetoNC")
function savetoNC(fsaven, d, b, vB, v, uptake, timet, Ctot, Cwall, PAall, tst, tfn, prms)

    println("Saving to: ",fsaven)
    
    #if isfile(fsaven)
    #    f = NCDataset(fsaven, "a") #a for annotate -- this COPIES OVER the previous
    #else
    f = NCDataset(fsaven, "c") #c for create
    #end

    nrec1 = Int(prms.nrec+1) #bc i added time 0

    defDim(f,"nb",prms.nb)
    defDim(f,"nd",prms.nd)
    defDim(f,"nrec",nrec1)
   
    f.attrib["title"] = "Community consumption model i/o"
    f.attrib["Start time"] = string(tst)
    f.attrib["End time"] = string(tfn)
    f.attrib["Run time"] = string(tfn - tst) 

    w = defVar(f,"d",Float64,("nd","nrec"))
    w[:,:] = d
    w.attrib["units"] = "mmol/m3 C OM"
    
    w = defVar(f,"b",Float64,("nb","nrec"))
    w[:,:] = b
    w.attrib["units"] = "mmol/m3 C biomass"
    
    w = defVar(f,"vB",Float64,("nd","nb"))
    w[:,:,:] = vB
    w.attrib["units"] = "per d; specific uptake matrix before PA"
    
    w = defVar(f,"v",Float64,("nd","nb"))
    w[:,:,:] = v
    w.attrib["units"] = "per d; specific uptake matrix with PA impact"
    
    w = defVar(f,"uptake",Float64,("nd","nb"))
    w[:,:,:] = uptake
    w.attrib["units"] = "mmol/m3 C per d; uptake matrix"
    
    w = defVar(f,"dIC",Float64,("nd",))
    w[:,:] = prms.dIC
    w.attrib["units"] = "mmol/m3 C OM"
    
    w = defVar(f,"bIC",Float64,("nb",))
    w[:,:] = prms.bIC
    w.attrib["units"] = "mmol/m3 C biomass"
    
    w = defVar(f, "Ctot", Float64, ("nrec",))
    w[:] = Ctot
    w.attrib["units"] = "mmol/m3 total C"
    
    w = defVar(f,"Cw",Float64,("nd",))
    w[:,:] = prms.Cw 
    w.attrib["units"] = "Ind C supply weight: probability"
    
    w = defVar(f,"Cwall",Float64,("nd","nrec"))
    w[:,:] = Cwall 
    w.attrib["units"] = "Ind C supply weight: over time"
    
    w = defVar(f,"PA",Float64,("nb",))
    w[:,:] = prms.PA 
    w.attrib["units"] = "Ind C supply weight: over time"
    
    w = defVar(f,"PAall",Float64,("nb","nrec"))
    w[:,:] = PAall 
    w.attrib["units"] = "Ind C supply weight: over time"
    
    w = defVar(f,"CM",Float64,("nd","nb"))
    w[:,:] = prms.CM
    w.attrib["units"] = "Community consumption Matrix"
    
    w = defVar(f,"pen",Float64,("nb",))
    w[:,:] = prms.pen
    w.attrib["units"] = "penalty"
    
    w = defVar(f, "timet", Float64, ("nrec",))
    w[:] = timet
    w.attrib["units"] = "days"
    
    w = defVar(f, "Pt", Float64, ())
    w[:] = prms.Pt
    w.attrib["units"] = "Total C supply"

    w = defVar(f, "vmaxB", Float64, ("nd",))
    w[:] = prms.vmaxB
    w.attrib["units"] = "per d; max uptake rate"
    
    w = defVar(f, "km", Float64, ("nd",))
    w[:] = prms.km
    w.attrib["units"] = "mmol/m3; half-sat"
    
    w = defVar(f, "y", Float64, ("nd",))
    w[:] = prms.y
    w.attrib["units"] = "mol B/mol C; yield"
    
    w = defVar(f, "mq", Float64, ("nb",))
    w[:] = prms.mq 
    w.attrib["units"] = "m3/mmol/d; quadr mort"
    
    w = defVar(f, "mlin", Float64, ("nb",))
    w[:] = prms.mlin 
    w.attrib["units"] = "per d; lin mort"
    
    w = defVar(f, "D", Float64, ())
    w[:] = prms.D
    w.attrib["units"] = "Dilution rate"

    w = defVar(f, "T", Float64, ())
    w[:] = prms.T
    w.attrib["units"] = "d; total time"

    close(f)
end
