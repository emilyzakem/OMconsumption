!11/19/20: Fortran code by Emily Zakem for Zakem et al. "A unified theory for organic matter accumulation"

!Note 1: The steady state files included (*_fSS.txt) are of the model solution illustrated in the manuscript, which was spun up for 6000 years (so these could be plotted directly to reproduce the model results).
!HOWEVER, the initial condition script below resets the initial DOC and DON concentrations to low values, so that the model will need to run for another 6000 years to reproduce the figures here. 

!Note 2: As written, this model has options for 3 different P types and 3 zooplankton types, but here, only two P (p1 and p2) are used. (No explicit zooplankton grazers.)

PROGRAM EZM

IMPLICIT NONE

INTEGER,PARAMETER :: ndays=3e6, &!5e5, &
	ndistr=25, & !25 max, & !size of VmaxA, the length of the distribution
    Hp=2000, & !change this scale below, too!
	dz=5, &  
	nz=Hp/dz
REAL*8,PARAMETER :: dt=0.02D0, & !0.02 for 5m res, 0.05 for 10mres
	!ndays=0.05D0, & !for fractions of days
	H=2D3, & !!!!!change these scales,too!!!!!
	!
	!PHYSICAL params:
	euz=25D0, & !m 
	mlz=20D0, & !m mixed layer depth 
	o2sat=.2121D0, & !mol/m3 from calc_oxsat(25+273,35) in matlab. WOCE clim-avg surf T at 10S, E. Pac.
	Kgast=3D-5, & !m/s
	!deep oxygen relaxation
	o2satdeep=0.2D0, & !mol/m3, avg (for ~7 C) and 35
	t_o2relax=1D-2, & !1/day -from 0.01 to 0.1 
	kappazmin=1D-4, & !m2/s -value for most of the deep ocean (higher at top and bottom)
	kappazmax=1D-2, & !m2/s -value for most of the deep ocean (higher at top and bottom)
	Ws=10D0, & !m/day !was 10- 3 from Anderson 2007
	!Lateral N transfer:
    latnConc=0D0/365D0*1D-3, & !from Hickman 2010: 0.2 uM N/yr total
    !
    !Stoichiometries:
    Rb_CN = 6.6, &
    Rz_CN = 5, &
    Rp_CN = 6.6, &
	!BACTERIA METABOLISMS: 
    !
	!Bo: Aerobic Heterotroph
    yd_bo=0.2D0, & !mol cells/mol Detritus - aerobic bacteria yield
	yo_bo=yd_bo/20D0*4D0/(1D0-yd_bo), & !mol cells/mol O2 -where cells have 1mol N 
	!enh4_bo=(1D0/yd_bo-1D0), & !production of ammonia per mol cells produced
	ec=(1D0/yd_bo-1D0), & !production of ammonia per mol cells produced
	pd_max= 1D0, & !1/day - normalized max detritus uptake to match nitrate max uptake
	kd= 0.1D-3, & !org matter uptake half sat; used to be 5
	po_coef=2329100D0, & !m3/mol/day -
	alphaI=1.15D0, & !inverse of x% efficiency
	!alphaI=1D0, & !inverse of x% efficiency
    !
    !Litchman uptake of N for B and P (not nitrif for now)
    pn_maxL = 8.8D0, &!mol DIN/mol biomass N/day (VmaxN/QminN for diam = 0.5 Litchman)
    knL_nox = 0.081D-3, &!mol NOx/m3 (Kn Litchman)
    knL_nh4 = 0.041D-3, &!as in Ward, knL_nox/2
    !
	!Bnh4: Ammonia oxidizer
    ynh4_bnh4=1D0/112D0, &
    yo_bnh4=1D0/162D0, &
    eno2_bnh4=(1D0/ynh4_bnh4-1D0), & 
    pn_max= 40.6D0, & !was 30 before TempFun. 20D0, & !mol NO3 uptake/mol cellular N/day - normalized max nitrate uptake; see 2/16 notes
	kn= 0.133D-3, & !mol/m3 - DIN uptake half-sat; Litchman
	!
	!Bno2: Nitrite oxidizer
    yno2_bno2=1D0/334D0, &
    yo_bno2=1D0/162D0, &
	eno3_bno2=(1D0/yno2_bno2-1D0), & 
	pn_max_noo= pn_max/2.1544D0, & !NOB 5x higher quota !250D0, & !20D0, & !mol NO3 uptake/mol cellular N/day 
	kn_noo= kn*2.1544D0, & !NOB: 5x higher Q PLUS 10x lower affinity!. 0.5D-3, & !mol/m3 - DIN uptake half-sat; Litchman	
	!
	!P: phytoplankton
	Iinmax=1400D0, & !W/m2 
    !prochl:
    umaxp3=0.515D0, & !1/d -Ward 2014: 1*V^-0.15 for d=0.6 for Pro. max growth rate for p at 20C
	!one DIN pool:
	!knh4_effp3= 0.0018D-3, & !mol/m3 Litchman scaling, with aV^b with a for knh4 from Ward 2014 
    knh4_effp3= 0.0036D-3, & !mol/m3 Litchman scaling, with aV^b with a for knh4 from Ward 2014 
	knox_effp3= 0.0036D-3, & !mol/m3 Litchman
    !Diatom:
    umaxp2=3D0, & 
	!one DIN pool:
    !knh4_effp2= 0.164D-3, & !mol/m3 Litchman scaling, with aV^b with a for knh4 from Ward 2014 
    knh4_effp2= 0.327D-3, & !mol/m3 Litchman scaling, with aV^b with a for knh4 from Ward 2014 
    knox_effp2= 0.327D-3, &
    !
    amminhib=4.6*1e3*1D0, & !m3/mol NH4 (from 1/uM nh4 in follows 2007)
    chl2cmax=0.2D0, & !mg Chl/mmol C from Dutk 2015 
    chl2cmin=0.02D0, & !min Chl/mmol C 
    phimax=40D0, & !mmol C/mol photons/Ein  -quantum yield (mol C/mol photons)
    !a_chl=0.02D0, & !m2/mg chl a -absorption parameter by chlorophyll 
    a_chlp2=0.01D0, & !diatom
    a_chlp3=0.04D0, & !ll pro m2/mg chl a -absorption parameter (in paper as m2/mg chla) avg over all wavelengths 
    a_chlD=0D0,&!.04D3, & !m2/g chl a -for light attenutation: chlorophyll plus CDOM
	convI=2.77D18/6.02D23*86400D0, & !from W/m2 to Ein/m2/d 2.77e18[quanta/W/s]*1/6.02e23[Einstein/quanta]. from MBARI: http://www3.mbari.org/bog/nopp/par.html
    !
    !Zooplankton:
    gmax=0D0, & !1/d
    kg=1D-3, & !half-sat for grazing on pp
    gam=0.5D0, & !growth yield for zoo 
	mz=0.7D3, & !for quadratic mortality, so that mz*Z*Z = 1*0.1*0.1 and is like linear mz=0.1/day.
    !
	!GRAZINGandMORTALITy:
	!mlin=0D0, & !linear mortality for b and p
	fmlin = 0D0, &
    mlinA=pn_max*ynh4_bnh4*fmlin, &
    mlinN=pn_max*yno2_bno2*fmlin, &
    mlinB=pd_max*yd_bo*fmlin, & !mlinBall below
    mlinp2=umaxp2*fmlin, &
    mlinp3=umaxp3*fmlin, &
    mquad=1D3, & !quadratic mortality for b and p
    mortf=0D0, & !fraction of mort to DOM vs POM
    !maintf=0.5D0, &!fraction of loss to maintenance (respired to inorganics) vs mortality (to DOM/POM)
    !
	!Oxygen ratio for pp production and zoo consumption:
    RredO=1D0/(yo_bo*(1D0/yo_bo-1D0)), & !implied ratio of export production to balance surface oxygen (see 3/27/15 notes)
    !Temperature:
    TempAeArr = -4D3, &
    TemprefArr = 293.15D0, &    
    Tkel = 273.15D0, &
    TempCoeffArr = 0.8D0

    INTEGER :: t,mlboxes,j,i,jc,nt,startSS,ind,recordDaily,dailycycle,dommodel
REAL*8 :: zm(nz) = (/(j,j=0+dz/2,Hp-dz/2, dz)/)
REAL*8 :: z(nz+1) = (/(j,j=0,Hp, dz)/)
REAL*8 :: koverh, distn, adv, diff, cputime,dayint,dayint2,Iin,Cghost
REAL*8,DIMENSION(:),ALLOCATABLE :: time
REAL*8,DIMENSION(:,:),ALLOCATABLE :: sumall
REAL*8,DIMENSION(nz+1) :: w,wd,Kz,KzO
REAL*8,DIMENSION(nz+4) :: eqmask,inmask,Iz,latnout,latn, & 
	u_bo,u_bo_pa,u_bnh4,u_bno2,u_p1,u_p2,u_p3,u_p3nolim,bt,pt,btsq,ptsq,g,g2,g3, &
    Biot,Chlt,po,pd,pdom,pnh4b,pno2b,nlimtot,Temp,Q10r,TempFun, &  !these are replaced each time, don't need k and A,B,C terms.
	limnh4,limno2,limno3,inhibnh4, & !explicit limits to ease equations later
	limnh4b,limno2b,limnh4p3,limno2p3,limno3p3, & !different n lims for the b
    nh4,no2,no3,ntot,o,bo,bo_pa,bnh4,bno2,d,dom,zoo,zoo2,zoo3,p1,xp1,p2,xp2,p3,xp3, &
	knh4A,knh4B,knh4C,knh4D,nh4A,nh4B,nh4C, &
	kno2A,kno2B,kno2C,kno2D,no2A,no2B,no2C, &
	kno3A,kno3B,kno3C,kno3D,no3A,no3B,no3C, &
	koA,koB,koC,koD,oA,oB,oC, &
	kboA,kboB,kboC,kboD,boA,boB,boC, &
	kbo_paA,kbo_paB,kbo_paC,kbo_paD,bo_paA,bo_paB,bo_paC, &
	kbnh4A,kbnh4B,kbnh4C,kbnh4D,bnh4A,bnh4B,bnh4C, &
	kbno2A,kbno2B,kbno2C,kbno2D,bno2A,bno2B,bno2C, &
	kdA,kdB,kdC,kdD,dA,dB,dC,kdomA,kdomB,kdomC,kdomD,domA,domB,domC, &
	poc,kpocA,kpocB,kpocC,kpocD,pocA,pocB,pocC, &
	kzooA,kzooB,kzooC,kzooD,zooA,zooB,zooC, &
	kzoo2A,kzoo2B,kzoo2C,kzoo2D,zoo2A,zoo2B,zoo2C, &
	kzoo3A,kzoo3B,kzoo3C,kzoo3D,zoo3A,zoo3B,zoo3C, &
	kp1A,kp1B,kp1C,kp1D,p1A,p1B,p1C, &
	kxp1A,kxp1B,kxp1C,kxp1D,xp1A,xp1B,xp1C, &
	kp2A,kp2B,kp2C,kp2D,p2A,p2B,p2C, &
	kxp2A,kxp2B,kxp2C,kxp2D,xp2A,xp2B,xp2C, &
	kp3A,kp3B,kp3C,kp3D,p3A,p3B,p3C, &
	kxp3A,kxp3B,kxp3C,kxp3D,xp3A,xp3B,xp3C, &
    no3uptakeP,no2emitP, &
    Q1,Q2,Q3,Vnh4,Vno2,Vno3, & !for quota model (as well as xp variants above: p1 is like XQ (biomass), xp is like X)
    PC1,PCmax1,PC2,PCmax2,PC3,PCmax3,a_Ip2,a_Ip3,chl2c,chl2c_p1,chl2c_p2,chl2c_p3, & !for Chl:C model
    mortall,mortallC, &
    pNH4max,u_Norg,Rup, &
    pPON,pPOC,pNH4_bo,eDIC_bo,eNH4_bo, & !multiple hetB and POC/detritus in N form here
    dic,kdicA,kdicB,kdicC,kdicD,dicA,dicB,dicC !doc
REAL*8,DIMENSION(ndistr) :: Vmaxall,Vmaxall_w,yd_ball,yo_ball,kdall,mlinBall,burial !multiple hetB and POC/detritus in N form here
REAL*8,DIMENSION(nz+4,ndistr) :: boall,dall, &
    pDON,pDOC,u_boall,pNH4,eDIC,eNH4, & !multiple hetB and POC/detritus in N form here
    kboallA,kboallB,kboallC,kboallD,boallA,boallB,boallC, &
    kdallA,kdallB,kdallC,kdallD,dallA,dallB,dallC, &
    docall,kdocallA,kdocallB,kdocallC,kdocallD,docallA,docallB,docallC !doc

startSS=1
dayint=100D0 !for recording each timestep for resolved 1D for movies
dayint2=100D0
dailycycle=0 !for time-variant model
recordDaily=0 !for recording each dt for movies/2D


kdA(:)=0D0
kdB(:)=0D0
kdC(:)=0D0
kdD(:)=0D0
!
kpocA(:)=0D0
kpocB(:)=0D0
kpocC(:)=0D0
kpocD(:)=0D0
!
kdomA(:)=0D0
kdomB(:)=0D0
kdomC(:)=0D0
kdomD(:)=0D0
!
knh4A(:)=0D0
knh4B(:)=0D0
knh4C(:)=0D0
knh4D(:)=0D0
kno2A(:)=0D0
kno2B(:)=0D0
kno2C(:)=0D0
kno2D(:)=0D0
kno3A(:)=0D0
kno3B(:)=0D0
kno3C(:)=0D0
kno3D(:)=0D0
!
kboA(:)=0D0
kboB(:)=0D0
kboC(:)=0D0
kboD(:)=0D0
!
kbnh4A(:)=0D0
kbnh4B(:)=0D0
kbnh4C(:)=0D0
kbnh4D(:)=0D0
kbno2A(:)=0D0
kbno2B(:)=0D0
kbno2C(:)=0D0
kbno2D(:)=0D0
koA(:)=0D0
koB(:)=0D0
koC(:)=0D0
koD(:)=0D0
!
kzooA(:)=0D0
kzooB(:)=0D0
kzooC(:)=0D0
kzooD(:)=0D0
kzoo2A(:)=0D0
kzoo2B(:)=0D0
kzoo2C(:)=0D0
kzoo2D(:)=0D0
kzoo3A(:)=0D0
kzoo3B(:)=0D0
kzoo3C(:)=0D0
kzoo3D(:)=0D0
!
kp1A(:)=0D0
kp1B(:)=0D0
kp1C(:)=0D0
kp1D(:)=0D0
!
kp2A(:)=0D0
kp2B(:)=0D0
kp2C(:)=0D0
kp2D(:)=0D0
!
kp3A(:)=0D0
kp3B(:)=0D0
kp3C(:)=0D0
kp3D(:)=0D0
!
kboallA(:,:)=0D0
kboallB(:,:)=0D0
kboallC(:,:)=0D0
kboallD(:,:)=0D0
!
kdallA(:,:)=0D0
kdallB(:,:)=0D0
kdallC(:,:)=0D0
kdallD(:,:)=0D0
!
kdocallA(:,:)=0D0
kdocallB(:,:)=0D0
kdocallC(:,:)=0D0
kdocallD(:,:)=0D0
!
kdicA(:)=0D0
kdicB(:)=0D0
kdicC(:)=0D0
kdicD(:)=0D0
!

print*,'Run for total n of days:';print*,ndays
print*,'1D version'
print*,'nz is:'; print*,nz

nt=ndays/dt

print*,'Number of days:'
print*,ndays
print*,'Number of timesteps:'
print*,nt

!input the range of Vmax for each of the OM pools
!the SLOWER Vmaxes are first. the range is from 0.01 to 100
OPEN(UNIT=3,FILE='VmaxA_values_n25s2.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
    read(3,*) (Vmaxall(J),J=1,ndistr)    
    close(3)

!input the weight of production for each OM pool (based on lognormal distribution)
OPEN(UNIT=3,FILE='VmaxA_weights_n25s2.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
    read(3,*) (Vmaxall_w(J),J=1,ndistr)    
    close(3)

kdall=Vmaxall/10D3
!kdall=kd !increasing affinity with lability
yd_ball(:) = yd_bo
yo_ball=yd_ball/20D0*4D0/(1D0-yd_ball) !mol cells/mol O2 -where cells have 1mol N 

!mlinBall=min(0.1D0,Vmaxall*yd_ball*fmlin)
mlinBall=Vmaxall*yd_ball*fmlin

!now write for import into python:
OPEN(UNIT=5,FILE='VmaxA_fromF.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (Vmaxall(J),J=1,ndistr)
CLOSE(5)
OPEN(UNIT=5,FILE='VmaxA_w_fromF.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (Vmaxall_w(J),J=1,ndistr)
CLOSE(5)
OPEN(UNIT=5,FILE='kdall_fromF.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (kdall(J),J=1,ndistr)
CLOSE(5)
OPEN(UNIT=5,FILE='yd_ball_fromF.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (yd_ball(J),J=1,ndistr)
CLOSE(5)

print*,'ndistr:'
print*,ndistr
print*,'Vmaxall:'
print*,Vmaxall
print*,'Vmaxall weights:'
print*,Vmaxall_w
print*,'kdall:'
print*,kdall
print*,'yd_ball:'
print*,yd_ball
print*,'yo_ball:'
print*,yo_ball

mlboxes=100D0/dz !discrete n of boxes in the mixed layer, close to 100m total sum
koverh=Kgast/100D0/mlboxes *3600D0*24D0 !gas transfer coefficient for each of the n boxes comprising the ml

ALLOCATE(time(nt))
ind=ndays/dayint2
ALLOCATE(sumall(ind,20))

wd(:)=0D0!Ws
wd(1)=0D0
!if the following is commented out, then D will leave system
wd(nz+1)=0D0 !if this is not commented, the bottom box accumulates D

!1D model no advection:
	w(:)=0D0
	wd=wd+w; !vertical velocity combination for detritus

Temp(:)=0D0 !not sure if it matters but just in case

Kz=(kappazmax*exp(-z/mlz)+kappazmin)*3600D0*24D0 !larger at bottom boundary, too

Temp(3:nz+2)=12D0*exp(-zm/150D0)+12D0*exp(-zm/500D0)+2D0
TempFun = TempCoeffArr*exp(TempAeArr*(1D0/(Temp+Tkel)-1D0/TemprefArr))


KzO=Kz
KzO(1)=0D0
!for an open boundary: a sink for oxygen:
KzO(nz+1)=(kappazmin+1D-2*exp((H-H)/100D0))*3600D0*24D0 !the diffusion that would be on the bottom boundary, used for oxygen fixed conc./fixed flux in

!5/2016: for a closed boundary: a fixed o2
KzO(nz+1)=0D0

Kz(1)=0D0
Kz(nz+1)=0D0

!EJZ oct 2017: i turned this on for the first time
!Kz(nz)=0D0 !making last box accumulate D

eqmask(:)=0D0
eqmask(3:mlboxes+2)=1D0; !mask for air-sea equilibration

inmask(:)=0D0
inmask(3:nz+2)=1D0

!import previous steady state as IC:
!SSfiles

if (startSS.eq.1) then

		OPEN(UNIT=3,FILE='nh4_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (nh4(J),J=1,nz+4)
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='no2_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (no2(J),J=1,nz+4)
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='no3_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (no3(J),J=1,nz+4)
		CLOSE(3)

		OPEN(UNIT=3,FILE='d_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (d(J),J=1,nz+4)
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='dom_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (d(J),J=1,nz+4)
		CLOSE(3)

		OPEN(UNIT=3,FILE='o_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (o(J),J=1,nz+4)
		CLOSE(3)

		OPEN(UNIT=3,FILE='bo_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (bo(J),J=1,nz+4)
		CLOSE(3)

		OPEN(UNIT=3,FILE='bnh4_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (bnh4(J),J=1,nz+4)
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='bno2_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (bno2(J),J=1,nz+4)
		CLOSE(3)

		OPEN(UNIT=3,FILE='zoo_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (zoo(J),J=1,nz+4)
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='zoo2_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (zoo2(J),J=1,nz+4)
		CLOSE(3)
		
		OPEN(UNIT=3,FILE='zoo3_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (zoo3(J),J=1,nz+4)
		CLOSE(3)
		
        OPEN(UNIT=3,FILE='p1_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (p1(J),J=1,nz+4)
		CLOSE(3)

        OPEN(UNIT=3,FILE='p2_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (p2(J),J=1,nz+4)
		CLOSE(3)

        OPEN(UNIT=3,FILE='p3_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (p3(J),J=1,nz+4)
		CLOSE(3)

		OPEN(UNIT=3,FILE='boall_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
        DO I=1,ndistr
		read(3,*) (boall(J,I),J=1,nz+4)
        END DO
		CLOSE(3)

		OPEN(UNIT=3,FILE='dall_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
        DO I=1,ndistr
		read(3,*) (dall(J,I),J=1,nz+4)
        END DO
		CLOSE(3)

		OPEN(UNIT=3,FILE='docall_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
        DO I=1,ndistr
		read(3,*) (docall(J,I),J=1,nz+4)
        END DO
		CLOSE(3)

		OPEN(UNIT=3,FILE='dic_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (dic(J),J=1,nz+4)
		CLOSE(3)

		OPEN(UNIT=3,FILE='poc_fSS.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
		read(3,*) (poc(J),J=1,nz+4)
		CLOSE(3)

else
		
!NEW Initial Conditions
        !!Initial Conditions from simple distributions (may break now bc of steep gradients):
		nh4(:)=.01D-3
		no2(:)=.1D-3
		bo(:)=inmask*1e-6*1D0
		bnh4(:)=inmask*1e-6*1D0
		bno2(:)=inmask*1e-6*1D0
            p1(:)=inmask*1e-5*1D0
		    p2(:)=inmask*1e-5*1D0
            p3(:)=inmask*1e-5*1D0

		o(:)=inmask*0.1*1D0 !mol/m3 crude estimate 
			!set ghost cell row as fixed concentration:
			!o(nz+3)=100*1024/1e6*1D0 !simulates deep lateral advection/source of O2
		zoo(:)=inmask*1e-5*1D0
		zoo2(:)=inmask*1e-5*1D0
		zoo3(:)=inmask*1e-5*1D0
		!zoo3(:)=0D0
		
        no3(3:nz+2)=0.03*(1-exp(-zm/200))*1D0 ! n increases with depth


		!d(nz+2-50:nz+2)=0D0 !make no D at the start for the last 20*dz meters.
		d(:)=0D0 !start with none
		        !d(3:nz+2)=0.001*(exp(-zm/euz))*1D0 ! start with very little D at depth so it doesn't pile at bottom boundary

end if

!MODIFY
!INITIAL CONDITIONS

!take out p1 and p3:
p1(:)=0D0
!p2(:)=0D0

bo_pa(:)=0D0

!put p1 and p2 back in:
!p1=p3
p2=p3

!11/6/19: ONE DIN POOL: turn off nitrification by zeroing out nitrifiers
!no3(3:nz+2)=0.04*(1-exp(-zm/200))*1D0 ! n increases with depth
nh4 = nh4+no2+no3
nh4 = nh4*inmask

no2(:)=0D0
no3(:)=0D0

bnh4(:)=0D0
bno2(:)=0D0

!zoo:
!zoo(:)=inmask*1D-5
!zoo2(:)=inmask*1D-5
!zoo3(:)=inmask*1e-5
!
zoo(:)=0D0 !just b
zoo2(:)=0D0 !just p
zoo3(:)=0D0 !both p and b

!Apr 7: now in SS
!fill in all the Bs and Ds with IC (later change to real IC):
    !do I=1,ndistr
    !    boall(:,I)=inmask*1e-5!bo(:)*1D0/ndistr
    !    dall(:,I)=0D0!d(:)*1D0/ndistr
    !end do

!zero out the old
!d(:)=0D0
!bo(:)=0D0

!for original start of distr_update:
!bo(:)=inmask*1D-5
!dall(3:nz+2,:)=1D-3
!no3(3:nz+2)=40D-3*(1D0-exp(-zm/200D0))

! do I=1,ndistr
!     boall(:,I)=inmask*1e-6
! !    dall(:,I)=inmask*1e-5
! end do
! !docall = dall*Rb_CN
! 
! !APRIL 2020: just set POM AND DOM:
d = inmask*1e-3
poc = d*Rb_CN
! !
! !dic = inmask
! 
do I=1,ndistr
    docall(:,I) = inmask*1e-3
!     !docall(:,I) = inmask*1e-6
end do
! !then make first few large so total C is about right in whole model
! docall(:,1) = inmask*1e-2
! docall(:,2) = inmask*1e-2
! docall(:,3) = inmask*1e-2
! docall(:,4) = inmask*1e-2
dall = docall/Rb_CN

OPEN(UNIT=5,FILE='time_record.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)

if (recordDaily.eq.1) then

OPEN(UNIT=5,FILE='Iin.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)

OPEN(UNIT=5,FILE='Iz.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)

OPEN(UNIT=5,FILE='time_all.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)

OPEN(UNIT=5,FILE='nh4_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',STATUS='REPLACE')
CLOSE(5)	
OPEN(UNIT=5,FILE='no2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',STATUS='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='no3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',STATUS='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='d_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='dom_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='o_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)	
OPEN(UNIT=5,FILE='bo_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)	
OPEN(UNIT=5,FILE='bo_pa_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)	
OPEN(UNIT=5,FILE='bnh4_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='bno2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='ubo_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='ubnh4_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='ubno2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='up_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='zoo_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='zoo2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='zoo3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='p1_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='xp1_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='p2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='p3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='xp2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='xp3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)
OPEN(UNIT=5,FILE='boall_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)	
OPEN(UNIT=5,FILE='dall_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',status='REPLACE')
CLOSE(5)	


end if 


print *,'Starting time loop:'
do t=1,nt

!LIGHT:
if (dailycycle.eq.1) then !Incoming light daily cycle:
    Iin=Iinmax/2D0*(cos(t*dt*2D0*3.1416D0)+1D0)
else !No daily cycle:
    Iin=Iinmax/2D0
end if

!chl impact on light:
Chlt=(p1*chl2c_p1+p2*chl2c_p2+p3*chl2c_p3)*6.6D0 !molN/m3 *mgChl/mmolC *5molC/molN = gChl

do j=1,nz
    Iz(j+2)=Iin*exp(-zm(j)*(1/euz+sum(Chlt(3:j+2)*a_chlD))) !this one sums the k at each depth
end do

i=dayint2*1000
j=t*dt*1000
    
if (MOD(j,i).eq.0) then
	!trace the integral:
	ind=t*dt/dayint2
	
	!print*,ind


	time(ind)=(t-1)*dt !start at zero- weird just bc i can't just easily write out the vector

	ntot=nh4+no2+no3
	
	sumall(ind,1)=sum(o)*dz !mol/m3 times volume (with dx=1,dy=1)
	sumall(ind,2)=sum(d)*dz
	sumall(ind,3)=sum(dall)*dz
!	sumall(ind,3)=sum(dom)*dz
	sumall(ind,4)=sum(bo)*dz
	sumall(ind,5)=sum(boall)*dz
!	sumall(ind,5)=sum(bo_pa)*dz
	sumall(ind,6)=sum(bnh4)*dz
	sumall(ind,7)=sum(bno2)*dz
	sumall(ind,8)=sum(zoo)*dz
	sumall(ind,9)=sum(ntot)*dz
	sumall(ind,10)=sum(nh4)*dz
	sumall(ind,11)=sum(no2)*dz
	sumall(ind,12)=sum(no3)*dz
	sumall(ind,13)=sum(p1)*dz
	sumall(ind,14)=sum(p2)*dz
	sumall(ind,15)=sum(p3)*dz
	sumall(ind,16)=sum(zoo2)*dz
	sumall(ind,17)=sum(zoo3)*dz
	sumall(ind,18)=sum(docall)*dz
	sumall(ind,19)=sum(dic)*dz
	sumall(ind,20)=sum(poc)*dz

end if

call MYRK(nh4,no3,no2,d,dall,docall,poc,dic,dom,o,zoo,zoo2,zoo3,p1,xp1,p2,xp2,p3,xp3, &
				bo,boall,bo_pa,bnh4,bno2, & 
				knh4A,kno3A,kno2A,kdA,kdallA,kdocallA,kpocA,kdicA,kdomA,koA,kzooA,kzoo2A,kzoo3A, &
                kp1A,kxp1A,kp2A,kxp2A,kp3A,kxp3A, &
				kboA,kboallA,kbo_paA,kbnh4A,kbno2A)

nh4A = nh4 + dt/2D0*knh4A; 
boA = bo + dt/2D0*kboA; 
boallA = boall + dt/2D0*kboallA; 
dA = d + dt/2D0*kdA; 
pocA = poc + dt/2D0*kpocA; 
dallA = dall + dt/2D0*kdallA; 
docallA = docall + dt/2D0*kdocallA; 
dicA = dic + dt/2D0*kdicA; 
oA = o + dt/2D0*koA; 
p2A = p2 + dt/2D0*kp2A;
p3A = p3 + dt/2D0*kp3A;


call MYRK(nh4A,no3A,no2A,dA,dallA,docallA,pocA,dicA,domA,oA,zooA,zoo2A,zoo3A,p1A,xp1A,p2A,xp2A,p3A,xp3A, &
				boA,boallA,bo_paA,bnh4A,bno2A, & 
				knh4B,kno3B,kno2B,kdB,kdallB,kdocallB,kpocB,kdicB,kdomB,koB,kzooB,kzoo2B,kzoo3B, &
                kp1B,kxp1B,kp2B,kxp2B,kp3B,kxp3B, &
				kboB,kboallB,kbo_paB,kbnh4B,kbno2B)
				
nh4B = nh4 + dt/2D0*knh4B; 
boB = bo + dt/2D0*kboB; 
boallB = boall + dt/2D0*kboallB; 
dB = d + dt/2D0*kdB; 
pocB = poc + dt/2D0*kpocB; 
dallB = dall + dt/2D0*kdallB; 
docallB = docall + dt/2D0*kdocallB; 
dicB = dic + dt/2D0*kdicB; 
oB = o + dt/2D0*koB; 
p2B = p2 + dt/2D0*kp2B;
p3B = p3 + dt/2D0*kp3B;

call MYRK(nh4B,no3B,no2B,dB,dallB,docallB,pocB,dicB,domB,oB,zooB,zoo2B,zoo3B,p1B,xp1B,p2B,xp2B,p3B,xp3B, &
				boB,boallB,bo_paB,bnh4B,bno2B, & 
				knh4C,kno3C,kno2C,kdC,kdallC,kdocallC,kpocC,kdicC,kdomC,koC,kzooC,kzoo2C,kzoo3C, &
                kp1C,kxp1C,kp2C,kxp2C,kp3C,kxp3C, &
				kboC,kboallC,kbo_paC,kbnh4C,kbno2C)
				
nh4C = nh4 + dt*knh4C; 
boC = bo + dt*kboC; 
boallC = boall + dt*kboallC; 
dC = d + dt*kdC; 
pocC = poc + dt*kpocC; 
dallC = dall + dt*kdallC; 
docallC = docall + dt*kdocallC; 
dicC = dic + dt*kdicC; 
oC = o + dt*koC; 
p2C = p2 + dt*kp2C;
p3C = p3 + dt*kp3C;

call MYRK(nh4C,no3C,no2C,dC,dallC,docallC,pocC,dicC,domC,oC,zooC,zoo2C,zoo3C,p1C,xp1C,p2C,xp2C,p3C,xp3C, &
				boC,boallC,bo_paC,bnh4C,bno2C, & 
				knh4D,kno3D,kno2D,kdD,kdallD,kdocallD,kpocD,kdicD,kdomD,koD,kzooD,kzoo2D,kzoo3D, &
                kp1D,kxp1D,kp2D,kxp2D,kp3D,kxp3D, &
				kboD,kboallD,kbo_paD,kbnh4D,kbno2D)
	
nh4 = nh4 + dt/6D0*(knh4A + 2D0*knh4B + 2D0*knh4C + knh4D);
!no2 = no2 + dt/6D0*(kno2A + 2D0*kno2B + 2D0*kno2C + kno2D);
!no3 = no3 + dt/6D0*(kno3A + 2D0*kno3B + 2D0*kno3C + kno3D);
bo = bo + dt/6D0*(kboA + 2D0*kboB + 2D0*kboC + kboD);
boall = boall + dt/6D0*(kboallA + 2D0*kboallB + 2D0*kboallC + kboallD);
!bnh4 = bnh4 + dt/6D0*(kbnh4A + 2D0*kbnh4B + 2D0*kbnh4C + kbnh4D);
!bno2 = bno2 + dt/6D0*(kbno2A + 2D0*kbno2B + 2D0*kbno2C + kbno2D);
d = d + dt/6D0*(kdA + 2D0*kdB + 2D0*kdC + kdD);
dall = dall + dt/6D0*(kdallA + 2D0*kdallB + 2D0*kdallC + kdallD);
docall = docall + dt/6D0*(kdocallA + 2D0*kdocallB + 2D0*kdocallC + kdocallD);
dic = dic + dt/6D0*(kdicA + 2D0*kdicB + 2D0*kdicC + kdicD);
poc = poc + dt/6D0*(kpocA + 2D0*kpocB + 2D0*kpocC + kpocD);
o = o + dt/6D0*(koA + 2D0*koB + 2D0*koC + koD);	
!zoo = zoo + dt/6D0*(kzooA + 2D0*kzooB + 2D0*kzooC + kzooD);	
!zoo2 = zoo2 + dt/6D0*(kzoo2A + 2D0*kzoo2B + 2D0*kzoo2C + kzoo2D);	
!zoo3 = zoo3 + dt/6D0*(kzoo3A + 2D0*kzoo3B + 2D0*kzoo3C + kzoo3D);	
p1 = p1 + dt/6D0*(kp1A + 2D0*kp1B + 2D0*kp1C + kp1D);	
p2 = p2 + dt/6D0*(kp2A + 2D0*kp2B + 2D0*kp2C + kp2D);	
p3 = p3 + dt/6D0*(kp3A + 2D0*kp3B + 2D0*kp3C + kp3D);	

if (recordDaily.eq.1) then

!append at every time step:
print*,t*dt		
	OPEN(UNIT=7,FILE='Iin.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
	WRITE(7,*) (Iin)
	CLOSE(7)

OPEN(UNIT=5,FILE='Iz.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (Iz)
CLOSE(5)	

print*,t*dt		
	OPEN(UNIT=7,FILE='time_all.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
	WRITE(7,*) (t*dt)
	CLOSE(7)	

OPEN(UNIT=5,FILE='nh4_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (nh4)
CLOSE(5)	
OPEN(UNIT=5,FILE='no2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (no2)
CLOSE(5)
OPEN(UNIT=5,FILE='no3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (no3)
CLOSE(5)
OPEN(UNIT=5,FILE='d_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (d)
CLOSE(5)
OPEN(UNIT=5,FILE='dom_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (dom)
CLOSE(5)
OPEN(UNIT=5,FILE='o_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (o)
CLOSE(5)	
OPEN(UNIT=5,FILE='bo_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (bo)
CLOSE(5)	
OPEN(UNIT=5,FILE='bo_pa_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (bo_pa)
CLOSE(5)	
OPEN(UNIT=5,FILE='bnh4_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (bnh4)
CLOSE(5)
OPEN(UNIT=5,FILE='bno2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (bno2)
CLOSE(5)
OPEN(UNIT=5,FILE='ubo_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (u_bo)
CLOSE(5)
OPEN(UNIT=5,FILE='ubnh4_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (u_bnh4)
CLOSE(5)
OPEN(UNIT=5,FILE='ubno2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (u_bno2)
CLOSE(5)
OPEN(UNIT=5,FILE='up1_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (u_p1)
CLOSE(5)
OPEN(UNIT=5,FILE='up2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (u_p2)
CLOSE(5)
OPEN(UNIT=5,FILE='up3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (u_p3)
CLOSE(5)
OPEN(UNIT=5,FILE='zoo_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (zoo)
CLOSE(5)
OPEN(UNIT=5,FILE='zoo2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (zoo2)
CLOSE(5)
OPEN(UNIT=5,FILE='zoo3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (zoo3)
CLOSE(5)
OPEN(UNIT=5,FILE='p1_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (p1)
CLOSE(5)	
OPEN(UNIT=5,FILE='xp1_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (xp1)
CLOSE(5)	
OPEN(UNIT=5,FILE='p2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (p2)
CLOSE(5)	
OPEN(UNIT=5,FILE='xp2_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (xp2)
CLOSE(5)	
OPEN(UNIT=5,FILE='p3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (p3)
CLOSE(5)	
OPEN(UNIT=5,FILE='xp3_fa.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
WRITE(5,*) (xp3)
CLOSE(5)	

end if

if (MOD(t*dt,10.00).eq.0) then
	print*,t*dt	
    print*,nh4(3)*1D3    
	OPEN(UNIT=7,FILE='time_record.txt',ACCESS='SEQUENTIAL',BLANK='ZERO',POSITION='APPEND')
	WRITE(7,*) (t*dt)
	CLOSE(7)	
end if

if ((MOD(t*dt,dayint).eq.0).OR.(t*dt.eq.ndays)) then

    
OPEN(UNIT=5,FILE='chl2c_p1.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (chl2c_p1(J),J=1,nz+4)
CLOSE(5)
OPEN(UNIT=5,FILE='chl2c_p2.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (chl2c_p2(J),J=1,nz+4)
CLOSE(5)
OPEN(UNIT=5,FILE='chl2c_p3.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (chl2c_p3(J),J=1,nz+4)
CLOSE(5)
OPEN(UNIT=5,FILE='chl2c.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (chl2c(J),J=1,nz+4)
CLOSE(5)

!diff -- using no3uptakeP as placeholder
    do j=1,nz ; jc=j+2;
        call mydiff(no3,Kz,j,dz,nz,diff)
            no3uptakeP(jc)=diff
    end do
OPEN(UNIT=5,FILE='kno3_diff.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (no3uptakeP(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='wd.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (wd(J),J=1,nz+1)
CLOSE(5)

OPEN(UNIT=5,FILE='kz.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (kz(J),J=1,nz+1)
CLOSE(5)

OPEN(UNIT=5,FILE='Iz.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (Iz(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='nh4_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (nh4(J),J=1,nz+4)
CLOSE(5)
	
OPEN(UNIT=5,FILE='no2_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (no2(J),J=1,nz+4)
CLOSE(5)
	
OPEN(UNIT=5,FILE='no3_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (no3(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='d_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (d(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='dom_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (dom(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='o_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (o(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='bo_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (bo(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='ubo_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (u_bo(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='ubnh4_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (u_bnh4(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='ubno2_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (u_bno2(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='up1_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (u_p1(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='up2_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (u_p2(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='up3_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (u_p3(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='bo_pa_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (bo_pa(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='ubo_pa_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (u_bo_pa(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='bnh4_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (bnh4(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='bno2_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (bno2(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='zoo_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (zoo(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='zoo2_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (zoo2(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='zoo3_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (zoo3(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='p1_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (p1(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='p2_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (p2(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='p3_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (p3(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='time_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO J=1,ind
WRITE(5,*) (time(J))
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='z_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,nz
WRITE(5,*) (zm(I))
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='sumall_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ind
WRITE(5,*) (sumall(I,J),J=1,20)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='boall_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ndistr
WRITE(5,*) (boall(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='uboall_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ndistr
WRITE(5,*) (u_boall(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='dall_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ndistr
WRITE(5,*) (dall(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='docall_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ndistr
WRITE(5,*) (docall(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='dic_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (dic(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='poc_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (poc(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='pDOC_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ndistr
WRITE(5,*) (pDOC(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='pDON_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ndistr
WRITE(5,*) (pDON(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='pNH4_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ndistr
WRITE(5,*) (pNH4(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='eNH4_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ndistr
WRITE(5,*) (eNH4(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='eDIC_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
DO I=1,ndistr
WRITE(5,*) (eDIC(J,I),J=1,nz+4)
END DO
CLOSE(5)

OPEN(UNIT=5,FILE='pPOC_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (pPOC(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='pPON_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (pPON(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='pNH4_bo_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (pNH4_bo(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='eNH4_bo_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (eNH4_bo(J),J=1,nz+4)
CLOSE(5)

OPEN(UNIT=5,FILE='eDIC_bo_f.txt',ACCESS='SEQUENTIAL',BLANK='ZERO')
WRITE(5,*) (eDIC_bo(J),J=1,nz+4)
CLOSE(5)


	end if  !end the mod

end do !!end time loop

print*,'Total CPU time in seconds:'
call CPU_TIME(cputime)
print*,cputime


contains
SUBROUTINE MYQUICK(C,w,j,dz,nz,adv)
implicit none
REAL*8 :: C(nz+4),w(nz+1),adv
REAL*8 :: wp1,wn1,wp,wn,Dy1,Dy2,Dyn1,Fu,Fd
INTEGER :: j,nz,dz
INTEGER :: jc

jc=j+2

        !at top face:
        wp1=(w(j+1)+abs(w(j+1)))/2D0;
        wn1=(w(j+1)-abs(w(j+1)))/2D0;
        !at bottom face:
        wp=(w(j)+abs(w(j)))/2D0;
        wn=(w(j)-abs(w(j)))/2D0;
          
        Dy1=C(jc+2)-2D0*C(jc+1)+C(jc);
        Dy2=C(jc+1)-2D0*C(jc)+C(jc-1);
        Dyn1=C(jc)-2D0*C(jc-1)+C(jc-2);

        Fu=w(j+1)/2D0*(C(jc)+C(jc+1)) - wp1/8D0*Dy2 - wn1/8D0*Dy1;
        Fd=w(j)/2D0*(C(jc-1)+C(jc)) - wp/8D0*Dyn1 - wn/8D0*Dy2;
                
		adv=(Fu-Fd)/dz;  
        
END SUBROUTINE MYQUICK

SUBROUTINE MYDIFF(C,Kz,j,dz,nz,diff)
implicit none
REAL*8 :: C(nz+4),Kz(nz+1),diff
REAL*8 :: Fu,Fd
INTEGER :: j,nz,dz
INTEGER :: jc

jc=j+2
        
        Fu=Kz(j+1)*(C(jc+1)-C(jc))/dz;
        Fd=Kz(j)*(C(jc)-C(jc-1))/dz;
               
        diff=(Fu-Fd)/dz;
        
END SUBROUTINE MYDIFF


SUBROUTINE MYRK(nh4_one,no3_one,no2_one,d_one,dall_one,docall_one,poc_one,dic_one,dom_one,o_one, &
                zoo_one,zoo2_one,zoo3_one, &
				p1_one,xp1_one,p2_one,xp2_one,p3_one,xp3_one, &
				bo_one,boall_one,bo_pa_one,bnh4_one,bno2_one, & 
				knh4_two,kno3_two,kno2_two,kd_two,kdall_two,kdocall_two,kpoc_two,kdic_two,kdom_two,ko_two, &
                kzoo_two,kzoo2_two,kzoo3_two, &
				kp1_two,kxp1_two,kp2_two,kxp2_two,kp3_two,kxp3_two, &
				kbo_two,kboall_two,kbo_pa_two,kbnh4_two,kbno2_two)
				
implicit none
REAL*8, dimension(nz+4), intent(in) :: nh4_one,no3_one,no2_one,d_one,dom_one,o_one, &
				zoo_one,zoo2_one,zoo3_one,p1_one,xp1_one,p2_one,xp2_one,p3_one,xp3_one, &
				bo_one,bo_pa_one,bnh4_one,bno2_one,poc_one,dic_one
REAL*8, dimension(nz+4,ndistr), intent(in) :: boall_one,dall_one,docall_one
REAL*8, dimension(nz+4), intent(out) :: knh4_two,kno3_two,kno2_two,kd_two, &
				kdom_two,ko_two,kzoo_two,kzoo2_two,kzoo3_two, &
                kp1_two,kxp1_two,kp2_two,kxp2_two,kp3_two,kxp3_two, &
				kbo_two,kbo_pa_two,kbnh4_two,kbno2_two, &
                kpoc_two,kdic_two
REAL*8, dimension(nz+4,ndistr), intent(out) :: kboall_two,kdall_two,kdocall_two

    do j=1,nz ; jc=j+2;       
    	call mydiff(nh4_one,Kz,j,dz,nz,diff)
	knh4_two(jc)=diff
	!call mydiff(no2_one,Kz,j,dz,nz,diff)
	!kno2_two(jc)=diff
	!call mydiff(no3_one,Kz,j,dz,nz,diff)
        !kno3_two(jc)=diff
	call mydiff(bo_one,Kz,j,dz,nz,diff)
	!call myquick(bo_one,wd,j,dz,nz,adv) !sinking PA
	kbo_two(jc)=diff
	!call myquick(bo_pa_one,wd,j,dz,nz,adv)
	!call mydiff(bo_pa_one,Kz,j,dz,nz,diff)
        !kbo_pa_two(jc)=-adv+diff
        !call mydiff(bnh4_one,Kz,j,dz,nz,diff)
        !kbnh4_two(jc)=diff
        !call mydiff(bno2_one,Kz,j,dz,nz,diff)
        !kbno2_two(jc)=diff
        call myquick(d_one,wd,j,dz,nz,adv)
        call mydiff(d_one,Kz,j,dz,nz,diff)
        kd_two(jc)=-adv+diff
        call myquick(poc_one,wd,j,dz,nz,adv)
        call mydiff(poc_one,Kz,j,dz,nz,diff)
        kpoc_two(jc)=-adv+diff
        call mydiff(dic_one,Kz,j,dz,nz,diff)
        kdic_two(jc)=diff
        call mydiff(dom_one,Kz,j,dz,nz,diff)
        kdom_two(jc)=diff
        call mydiff(o_one,KzO,j,dz,nz,diff)
        ko_two(jc)=diff
        !call mydiff(zoo_one,Kz,j,dz,nz,diff)
        !kzoo_two(jc)=diff
        call mydiff(p2_one,Kz,j,dz,nz,diff)
        kp2_two(jc)=diff
        call mydiff(p3_one,Kz,j,dz,nz,diff)
        kp3_two(jc)=diff
        !call mydiff(zoo2_one,Kz,j,dz,nz,diff)
        !kzoo2_two(jc)=diff
        !call mydiff(zoo3_one,Kz,j,dz,nz,diff)
        !kzoo3_two(jc)=diff
        do I=1,ndistr
            call mydiff(boall_one(:,I),Kz,j,dz,nz,diff)
            kboall_two(jc,I)=diff
            call mydiff(dall_one(:,I),Kz,j,dz,nz,diff)
            !call myquick(dall_one(:,I),wd,j,dz,nz,adv)
            kdall_two(jc,I)=diff
            call mydiff(docall_one(:,I),Kz,j,dz,nz,diff)
            kdocall_two(jc,I)=diff
        end do
    
end do


bt=bo_one+bnh4_one+bno2_one
pt=p2_one+p3_one
btsq=bo_one*bo_one+bnh4_one*bnh4_one+bno2_one*bno2_one
ptsq=p2_one*p2_one+p3_one*p3_one
    do I=1,ndistr
        bt=bt+boall_one(:,I)
        btsq=btsq+boall_one(:,I)*boall_one(:,I)
    end do

!heterotrophic uptake and growth:
!Note: Put TempFun on growth rates, not uptake rates
po=po_coef*o_one !not really uptake- already normalized into 1/day by diving by Qmin (above)
pNH4max=pn_maxL*(nh4_one/(nh4_one+knL_nox)) 

!POM-consumers
pPON=pd_max*(d_one/(d_one+kd))
pPOC=pd_max*(poc_one/(poc_one+kd*Rb_CN))
Rup = pPOC/pPON !ratio of SPECIFIC uptake

Rup(1:2)=0D0
Rup(nz+3:nz+4)=0D0

!u_bo = min(po*yo_bo,pd*yd_bo)*TempFun
u_bo = min(yd_bo*pPOC, &
           pPON + PNH4max, &
           yo_bo*po) &
           *TempFun

u_Norg = min(yd_bo*pPOC, &
           pPON, &
           yo_bo*po) &
           *TempFun

!actual uptake/growth using NH4
pNH4_bo=MIN(pNH4max,MAX(u_bo - u_Norg, 0D0))

!org uptake may be different if O2 is limiting, for example
!note that both are still SPECIFIC uptake, so mult by B-C or B-N below in eqns
pPOC=1D0/yd_bo*u_bo
pPON=pPOC/Rup !uptake is constrained to this ratio, even if uptake is less due to O2 lim

pPON(1:2)=0D0
pPON(nz+3:nz+4)=0D0

eDIC_bo=pPOC - u_bo
eNH4_bo=pPON - u_Norg

!DOM-consumers: Litchman uptake 


do I=1,ndistr
        
    pDON(:,I)=Vmaxall(I)*(dall(:,I)/(dall(:,I)+kdall(I))) !did 10D3, try 1D3
    pDOC(:,I)=Vmaxall(I)*(docall(:,I)/(docall(:,I)+kdall(I)*Rb_CN)) !did 10D3, try 1D3
    Rup = pDOC(:,I)/pDON(:,I)!ratio of SPECIFIC uptake, so mult each by biomass (C or N) thus Rb_CN

    Rup(1:2)=0D0
    Rup(nz+3:nz+4)=0D0

    u_boall(:,I)=min(yd_ball(I)*pDOC(:,I), & !C limitation
                     pDON(:,I) + pNH4max, & !N limitation
                     yo_ball(I)*po) & !O2 limitation
                     *TempFun

    u_Norg=min(yd_ball(I)*pDOC(:,I), & !C limitation
                     pDON(:,I), & !N limitation
                     yo_ball(I)*po) & !O2 limitation
                     *TempFun

    !actual uptake/growth using NH4
    pNH4(:,I)=MIN(pNH4max,MAX(u_boall(:,I) - u_Norg, 0D0))

    !org uptake may be different if O2 is limiting, for example
    !note that both are still SPECIFIC uptake, so mult by B-C or B-N below in eqns
    pDOC(:,I)=1D0/yd_ball(I)*u_boall(:,I)
    pDON(:,I)=pDOC(:,I)/Rup !uptake is constrained to this ratio, even if uptake is less due to O2 lim
    
    pDON(1:2,:)=0D0
    pDON(nz+3:nz+4,:)=0D0

    eDIC(:,I)=pDOC(:,I) - u_boall(:,I)
    eNH4(:,I)=pDON(:,I) - u_Norg

end do

!nitrifier uptake and growth:
!limnh4b=(nh4_one/(nh4_one+kn))!*TempFun
!limno2b=(no2_one/(no2_one+kn_noo))!*TempFun

pnh4b=pn_max*(nh4_one/(nh4_one+kn)) 
pno2b=pn_max_noo*(no2_one/(no2_one+kn_noo)) 

u_bnh4=min(pnh4b*ynh4_bnh4,po*yo_bnh4)*TempFun
u_bno2=min(pno2b*yno2_bno2,po*yo_bno2)*TempFun

!phytopl uptake and growth rate:
inhibnh4 = exp(-amminhib*nh4_one) !from GUD
!for 2:
limnh4=(nh4_one/(nh4_one+knh4_effp2))!*TempFun
limno2=(no2_one/(no2_one+knox_effp2))*inhibnh4!*TempFun
limno3=(no3_one/(no3_one+knox_effp2))*inhibnh4!*TempFun
!u_p2=umaxp2*TempFun*min(Iz/(Iz+kI),limnh4+limno2+limno3)!*exp(-Iz*.01D0) !multiply by 0.9 for penalty?
PCmax2 = umaxp2*min(1D0,limnh4+limno2+limno3)*TempFun

!again now for 3:
limnh4p3=(nh4_one/(nh4_one+knh4_effp3))!*TempFun
limno2p3=(no2_one/(no2_one+knox_effp3))*inhibnh4!*TempFun
limno3p3=(no3_one/(no3_one+knox_effp3))*inhibnh4!*TempFun
nlimtot=limnh4p3+limno2p3+limno3p3
PCmax3 = umaxp3*min(1D0,nlimtot)*TempFun

!Geider chlorophyll:c based growth rates:
a_Ip2 = phimax*a_chlp2*Iz*convI !mmol C/mol Ein * m2/mg chla * Ein/m2/d = mmol C/mg chla/d
a_Ip3 = phimax*a_chlp3*Iz*convI !mmol C/mol Ein * m2/mg chla * Ein/m2/d = mmol C/mg chla/d
chl2c_p1 = 0D0!max(chl2cmin, min(chl2cmax, chl2cmax/(1D0+chl2cmax*a_I/2D0/PCmax1)))
chl2c_p2 = max(chl2cmin, min(chl2cmax, chl2cmax/(1D0+chl2cmax*a_Ip2/2D0/PCmax2)))*inmask
chl2c_p3 = max(chl2cmin, min(chl2cmax, chl2cmax/(1D0+chl2cmax*a_Ip3/2D0/PCmax3)))*inmask
!chl2c = chl2c_p3
!PC1 = PCmax1*(1D0 - exp(-a_I*chl2c_p1/PCmax1))
PC2 = PCmax2*(1D0 - exp(-a_Ip2*chl2c_p2/PCmax2))
PC3 = PCmax3*(1D0 - exp(-a_Ip3*chl2c_p3/PCmax3))

PC3(nz+2:nz+4)=0D0
PC3(1:2)=0D0
PC2(nz+2:nz+4)=0D0
PC2(1:2)=0D0
!PC1(nz+2:nz+4)=0D0
!PC1(1:2)=0D0

u_p1=0D0!PC1
u_p2=PC2
u_p3=PC3

!grazing: type II
g=gmax*bt/(bt+kg)*TempFun !for zoo
g2=gmax*pt/(pt+kg)*TempFun !for zoo2
g3=0D0!gmax*(bt+pt)/(bt+pt+kg)*TempFun !for zoo3 (goal: ONLY this one)

!lateral N flux:
latnout(:)=0D0
latn=inmask*latnConc
!mlz weighted:
latn(3:nz+2)=latn(3:nz+2)*exp(-zm/mlz)
!out at bottom in dom form:
latnout(nz+2)=sum(latn(3:nz+2))

!mortality:
mortall = mlinB*bo_one + mlinA*bnh4_one + mlinN*bno2_one &
        + mlinp2*p2_one + mlinp3*p3_one + mquad*(btsq + ptsq) &
        + mz*(zoo_one**2 + zoo2_one**2 + zoo3_one**2)
mortallC = Rb_CN*(mlinB*bo_one + mlinA*bnh4_one + mlinN*bno2_one) &
        + Rp_CN*(mlinp2*p2_one + mlinp3*p3_one) + mquad*(btsq*Rb_CN + ptsq*Rp_CN) & !THIS may cause trouble
        + Rz_CN*(mz*(zoo_one**2 + zoo2_one**2 + zoo3_one**2))

do I=1,ndistr
    mortall = mortall + mlinBall(I)*boall_one(:,I) !mquad is in btsq
    mortallC = mortallC + mlinBall(I)*boall_one(:,I)*Rb_CN !mquad is in btsq
end do

!EQUATIONS:

kdic_two = kdic_two & 
            !+ Rb_CN*ec*(u_bo*bo_one + u_bo_pa*bo_pa_one) &
            + Rb_CN*eDIC_bo*bo_one &
            + Rz_CN*(1D0-gam)*g*zoo_one &  !zoo
            + Rz_CN*(1D0-gam)*g2*zoo2_one  & !zoo2
            + Rz_CN*(1D0-gam)*g3*zoo3_one &  !zoo
            - Rb_CN*u_bnh4*bnh4_one & !consumption: NH4 oxidizer 
            - Rb_CN*u_bno2*bno2_one & !consumption: NH4 oxidizer 
            - Rp_CN*u_p1*p1_one & 
            - Rp_CN*u_p2*p2_one &
            - Rp_CN*u_p3*p3_one
            do I=1,ndistr
                kdic_two = kdic_two + Rb_CN*eDIC(:,I)*boall_one(:,I)
            end do

knh4_two = knh4_two &
			+ (eNH4_bo - pNH4_bo)*bo_one &
            + (1-gam)*g*zoo_one &  !zoo
			+ (1-gam)*g2*zoo2_one  & !zoo2
			+ (1-gam)*g3*zoo3_one &  !zoo
			- 1/ynh4_bnh4*u_bnh4*bnh4_one & !consumption: NH4 oxidizer 
            - u_p2*p2_one*limnh4/(limnh4+limno2+limno3+1D-38) & 
            - u_p3*p3_one*limnh4p3/(limnh4p3+limno2p3+limno3p3+1D-38) 
            !+ latn*nh4_one/(nh4_one+no2_one+no3_one+1D-38) !lateral transport  
            do I=1,ndistr
                knh4_two = knh4_two + (eNH4(:,I) - pNH4(:,I))*boall_one(:,I) 
            end do

kno2_two = kno2_two &
			+ eno2_bnh4*u_bnh4*bnh4_one & !source: NH4 oxidizer 
			- 1/yno2_bno2*u_bno2*bno2_one & !sink: aerobic NO2 oxidizer
            - u_p2*p2_one*limno2/(limnh4+limno2+limno3+1D-38) & 
            - u_p3*p3_one*limno2p3/(limnh4p3+limno2p3+limno3p3+1D-38) 
            !+ latn*no2_one/(nh4_one+no2_one+no3_one+1D-38)  !lateral transport  
		
kno3_two = kno3_two &
			+ eno3_bno2*u_bno2*bno2_one & !source: aerobic NO2 oxidizer
            - u_p2*p2_one*limno3/(limnh4+limno2+limno3+1D-38) & !no3uptakeP
            - u_p3*p3_one*limno3p3/(limnh4p3+limno2p3+limno3p3+1D-38) 
            !+ latn*no3_one/(nh4_one+no2_one+no3_one+1D-38)  !lateral transport  
	   						

kbo_two= kbo_two +  bo_one*(u_bo - mlinB*TempFun - mquad*bo_one*TempFun - g*zoo_one/(bt+1D-38) &
         - g3*zoo3_one/(bt+pt+1D-38)) !

do I=1,ndistr         
    kboall_two(:,I) = kboall_two(:,I) + boall_one(:,I)*(u_boall(:,I) - mlinBall(I)*TempFun &
            - mquad*boall_one(:,I)*TempFun &
            - g*zoo_one/(bt+1D-38) - g3*zoo3_one/(bt+pt+1D-38))
end do


kbnh4_two= kbnh4_two +  bnh4_one*(u_bnh4 - mlinA*TempFun - mquad*bnh4_one*TempFun - g*zoo_one/(bt+1D-38) &
           - g3*zoo3_one/(bt+pt+1D-38))
kbno2_two= kbno2_two +  bno2_one*(u_bno2 - mlinN*TempFun - mquad*bno2_one*TempFun - g*zoo_one/(bt+1D-38) &
           - g3*zoo3_one/(bt+pt+1D-38))

kp2_two = kp2_two + p2_one*(u_p2 - mlinp2*TempFun - mquad*p2_one*TempFun - g2*zoo2_one/(pt+1D-38) - g3*zoo3_one/(bt+pt+1D-38))
kp3_two = kp3_two + p3_one*(u_p3 - mlinp3*TempFun - mquad*p3_one*TempFun - g2*zoo2_one/(pt+1D-38) - g3*zoo3_one/(bt+pt+1D-38))

kzoo_two=kzoo_two + gam*g*zoo_one - mz*zoo_one*zoo_one*TempFun
kzoo2_two=kzoo2_two + gam*g2*zoo2_one - mz*zoo2_one*zoo2_one*TempFun
kzoo3_two=kzoo3_two + gam*g3*zoo3_one - mz*zoo3_one*zoo3_one*TempFun

!now d is POM and dall is DOM
kd_two= kd_two &
		+ (1D0-mortf)*mortall*TempFun &
		- alphaI*pPON*bo_one  !sink: 1 heterotroph

kpoc_two= kpoc_two &
		+ (1D0-mortf)*mortallC*TempFun &
		- alphaI*Rb_CN*pPOC*bo_one  !sink: 1 heterotroph


!DOM:
do I=1,ndistr
    !DON:
    kdall_two(:,I) = kdall_two(:,I)  &
		+ (Vmaxall_w(I)/sum(Vmaxall_w))*( & !this one is for the weighted distribution
		+ (alphaI-1D0)*pPON*bo_one & !sink: 1 heterotroph
        + mortf*mortall*TempFun) &
        - pDON(:,I)*boall_one(:,I)
    
    !DOC:
    kdocall_two(:,I) = kdocall_two(:,I) & 
		+ (Vmaxall_w(I)/sum(Vmaxall_w))*( & !this one is for the weighted distribution
		+ (alphaI-1D0)*Rb_CN*pPOC*bo_one &  !sink: 1 heterotroph
        + mortf*mortallC*TempFun) &
        - Rb_CN*pDOC(:,I)*boall_one(:,I)
end do

ko_two = ko_two &
		+ RredO*(u_p2*p2_one + u_p3*p3_one) & !pp production
        - RredO*(1D0-gam)*g*zoo_one & !zoo use
		- RredO*(1D0-gam)*g2*zoo2_one & !zoo2 use
		- RredO*(1D0-gam)*g3*zoo3_one & !zoo2 use
		- 1D0/yo_bo*u_bo*bo_one & !het B use
		- 1D0/yo_bnh4*u_bnh4*bnh4_one & !NH4 oxid use
		- 1D0/yo_bno2*u_bno2*bno2_one & !NO2 oxid use
		+ koverh*(o2sat-o_one)*eqmask & !air-sea
		+ t_o2relax*(o2satdeep-o_one)*inmask  !relaxation at depth (lateral flux)
            do I=1,ndistr
                ko_two = ko_two - 1D0/yo_ball(I)*u_boall(:,I)*boall_one(:,I)
            end do


END SUBROUTINE MYRK

END PROGRAM EZM


