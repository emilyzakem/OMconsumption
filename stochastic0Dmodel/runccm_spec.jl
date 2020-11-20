#11/19/20: Upload to Github

#Apr 19, 2020 by ejz
#For input into function runccm() within ccm.jl
#Parameters are assembled into type Prms (using "struct"), which is defined in ccm.jl 
#runccm then requires prms as input, i.e. a specific instance of the type Prms

#output file name (to which date and .nc will be added): 

fsave = "runccm_spec_test"

##############################################
#Parameters

using Statistics

T = 365*10 #days -- run time
nd = 1000 #n OM pools
nb = 1000 #n microbial pops
Pt = 1 #mmol C/m3 total OM supply

nrec = 1 # of timepoints to record (will be nrec + 1 with t=0)

#Consumption matrix
#spec:
CM = [i == j for i = 1:nd, j=1:nb] #diagonal only (specialists)
#specgen:
#CM = [i >= j for i = 1:nd, j=1:nb] #fills the lower diagonal (spec to gen gradient)
#random:
#CM = rand(nd,nb) .> 0.8 #randomly fills above specified value

nup = sum(CM,dims=1)[1,:]

#penalty for generalists:
pen = 1 ./ nup

#presence absence for each population (sets overall probability):
PA = ones(nb) #length nb
#PA = rand(nb)

#supply weight for each pool (sets overall probability)
Cw = rand(nd) #length nd 

#Dilution rate for additional sink (for D and B)
D = 0 #1/d 

#max uptake rate (1/d)
xhi = log10(1e2)
xlo = log10(1e-2) 
ρmaxB = 10 .^ ((rand(nd).-0.5).*(xhi-xlo).+mean([xhi xlo]))

#affinity (m3/mmol/d)
xhi = log10(100)
xlo = log10(1)
aff = 10 .^ ((rand(nd).-0.5).*(xhi-xlo).+mean([xhi xlo]))
km = ρmaxB./aff #half-sat constant (mmol/m3)

#yield, mortality
y = rand(nd)*0.5 #mol B/mol OM
mq = rand(nb)*1. #quadratic mort: m3/mmol/d
mlin = rand(nb)*1e-2 #linear mort: 1/d
#mlin = zeros(nb)
##############################################
#Initial conditions:

xhi = log10(10)
xlo = log10(0.1)
dIC = 10 .^ ((rand(nd).-0.5).*(xhi-xlo).+mean([xhi xlo]))

bIC = ones(nb)*0.1

##############################################
#load model including Prms struct:
include("ccm.jl")

##############################################
#assemble into type Prms (defined in ccm.jl)
params = Prms(fsave,T,Pt,nd,nb,CM,pen,PA,Cw,D,ρmaxB,km,y,mq,mlin,dIC,bIC,nrec);

##############################################
#run model and save into fsave_*date*.nc:
runccm(params);




