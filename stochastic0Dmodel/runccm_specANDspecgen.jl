#11/19/20: Upload to Github

#Apr 19, 2020 by ejz
#For input into function runccm() within ccm.jl
#Parameters are assembled into type Prms (using "struct"), which is defined in ccm.jl 
#runccm then requires prms as input, i.e. a specific instance of the type Prms

#output file name (to which date and .nc will be added): 

fsave = "runccm_specANDspecgen_test"

##############################################
#Parameters

using Statistics
using StatsBase

T = 365*10 #days -- run time
nd = 1000 #n OM pools
nb = nd*2 #n microbial pops
Pt = 0.1 #mmol C/m3 total OM supply

nrec = 1000 

#Consumption matrix


#FIRST HALF: Spec:
CMone = zeros(nd,nb÷2)
#CM = [i == j for i = 1:nd, j=1:nb÷2] #diagonal only (specialists)
#CM[i == j for i = 1:nd, j=1:nb÷2]
for j = 1:nd
    CMone[j,j] = 1
end

#SECOND HALF: random
##############################################################################
#METHOD 1: assign nup first

#1.1. nup = # of substrates consumed by each pop
#lin space:
nupr = rand(1:nd,nb÷2) # n substrates consumed by each, nb long

# #log space: 
# xhi = log10(nd)
# xlo = log10(1)
# nupr = 10 .^ ((rand(nb÷2).-0.5).*(xhi-xlo).+mean([xhi xlo]))
# nupr = round.(Int,nupr) #round to integers:

#nupr[:] .= nd/10 #to test only impact of # consumers 
#now for $ieach B, assign its substrates randomly: rand(1:nd,nup[j]) but UNIQUELY
#uniquely using sample: sample(1:nd,nb,replace=false) from StatsBase

#food=[sample(1:nd,nupr[j], replace = false) for j=1:nb] #an array of indices

#1.2. assign a likelihood of consumption for each substrate
w = [1/nd:1/nd:1; ] #ordered
wv = StatsBase.ProbabilityWeights(w)
#sample(1:nd,wv,10,replace=false) #over range 1:nd, with weight wv, give 10 samples, all unique

food=[sample(1:nd,wv,nupr[j], replace = false) for j=1:nb÷2] #an array of indices

CMtwo = zeros(nd,nb÷2)
for j = 1:nb÷2
    CMtwo[food[j],j] .= 1
end

##############################################################################
#METHOD 2: assign cms (# of consumers) first

# #1.1. nup = # of substrates consumed by each pop
# #lin space:
# #cmsr = rand(1:nb÷2,nd)
# 
# #log sapce:
# xhi = log10(nb÷2)
# xlo = log10(1)
# cmsr = 10 .^ ((rand(nd).-0.5).*(xhi-xlo).+mean([xhi xlo]))
# cmsr = round.(Int,cmsr) #round to integers
# 
# #1.2. assign a likelihood of generality for each population (how many substrates does it eat?)
# #w = [1/nd:1/nd:1; ] #ordered
# w = [1:nb÷2; ] #ordered -- same as above
# #w[:] .= 1
# #w = w .^ 4 #if this is increased, then things are a bit more evenly distributed (makes some MUCH less likely)
# wv = StatsBase.ProbabilityWeights(w)
# 
# #1.3. each is a list of the consumers for each substrate
# consumers=[sample(1:nb÷2,wv,cmsr[j], replace = false) for j=1:nd] #an array of indices
# 
# CMtwo = zeros(nd,nb÷2)
# for j = 1:nd
#     CMtwo[j,consumers[j]] .= 1
# end

#####################################
#concatenate

CM = hcat(CMone,CMtwo)

#####################################

#convert to boolean:
CM = convert(Array{Bool}, CM .== 1)

nup = sum(CM,dims=1)[1,:]

#penalty for generalists:
#pen = 1 ./ nup
pen = 1 ./ nup.^2
#pen = 1 ./ nup.^0.5

#presence absence for each population (sets overall probability):
#PA = ones(nb) #length nb
PA = rand(nb)

#supply weight for each pool (sets overall probability)
#Cw = ones(nd) #length nd 
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
mq = rand(nb)*1. .+ 0.1 #quadratic mort: m3/mmol/d
mlin = rand(nb)*1e-2 #linear mort: 1/d
#mq = ones(nb)*1. #quadratic mort: m3/mmol/d
#mlin = ones(nb)*1e-2 #linear mort: 1/d

##############################################
#Initial conditions:

xhi = log10(10)
xlo = log10(0.1)
dIC = 10 .^ ((rand(nd).-0.5).*(xhi-xlo).+mean([xhi xlo]))
#dIC = ones(nd)

bIC = ones(nb)*0.1

##############################################
#load model including Prms struct:
include("ccm.jl")
#include("ccm_switch.jl")

##############################################
#assemble into type Prms (defined in ccm.jl)
params = Prms(fsave,T,Pt,nd,nb,CM,pen,PA,Cw,D,ρmaxB,km,y,mq,mlin,dIC,bIC,nrec);

##############################################
#run model and save into fsave_*date*.nc:
runccm(params);



