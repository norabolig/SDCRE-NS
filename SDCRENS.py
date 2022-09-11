# Space Debris Collisional Rate Equations with N Shells (SDCRENS)
# 31 AUG 2022
# by norabolig (A. Boley)
#
#
import numpy as np
import sys
import copy

# preliminary setup
# read density model
atmfile=open("cira-2012.dat","r")
header=atmfile.readline()
zmodel=[]
denmodelL=[]
denmodelM=[]
denmodelHL=[]
denmodelHS=[]
for line in atmfile:
  alt,low,med,highL,highS=line.split()
  zmodel.append(float(alt))
  denmodelL.append(float(low))
  denmodelM.append(float(med))
  denmodelHL.append(float(highL))
  denmodelHS.append(float(highS))

atmfile.close()

zmodel=np.array(zmodel)*1000  # puts in metres
denmodelL=np.array(denmodelL)
denmodelM=np.array(denmodelM)
denmodelHL=np.array(denmodelHL)
denmodelHS=np.array(denmodelHS)

logdenL = np.log10(denmodelL)
logdenM = np.log10(denmodelM)
logdenHL = np.log10(denmodelHL)
logdenHS = np.log10(denmodelHS)
logz = np.log10(zmodel)

#
# read solar cycle template (using F10.7 as the solar activity index)
#
f107file = open("solar_cycle_table36_cira2012.dat","r")
header=f107file.readline()
f107_mo=[]
for line in f107file:
   mo,tmp,f107,tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp=line.split()
   f107_mo.append(float(f107))
f107file.close()

f107_mo=np.array(f107_mo) 

# GM of Earth
GM = 3.9860044418e14

# Radius of Earth at Equator
RE = 6378.137e3

# yr conversion
yr2sec = 365.25*24*3600

#
# This is the basic shell class
#
class NType:

    CD = 2.2
    AoM = None # m^2/kg
    tau = None # deorbit in seconds
    N = 0 # number of particles
    Tlife = None # lifetime in seconds
    size = None # diameter in m
    area = None
    denatm=None
    z=None
    v=None
    dz=None

    def setZ(self,z,dz):
        self.z = z
        self.dz = dz

    def setV(self):
        self.v = np.sqrt(GM/(RE+self.z))

    def setSize(self,s):
        self.size=s
        self.area=np.pi*s*s

# This is the lifetime of the satellites
    def setTlife(self,Tlife):
        self.Tlife=Tlife

# This is the drag time
    def setTau(self,denatm):
        adot = -self.CD*self.AoM*denatm*(self.z+RE)*self.v
        try: self.tau = np.abs(self.dz/adot)
        except:
           print("Error in setting Tlife. Likely adot is zero")
           sys.exit(-1)

#
# get the atmospheric density
#
def setDenatm(alt,t,mo0,setF107=None):

    i=int((alt-100e3)/20e3)
    if i > len(zmodel)-2: i=len(zmodel)-2
    if i < 0: i=0

    try:
       logalt = np.log10(alt)
    except:
       logalt = 0.
    mo_frac = (t/yr2sec)*12 + mo0

    mo = mo_frac % 144

    #print("Cycled mo = {} ".format(mo))

    moID = int(mo)

    if setF107==None: 
       moID1 = moID+1
       if moID1>143:moID1=0
       F107 = f107_mo[moID] + (f107_mo[moID1]-f107_mo[moID])*(mo-moID)
    else: F107 = setF107

    #print("F107 = {}".format(F107))

    if F107 <= 65:
       d = 10.**(  logdenL[i]+(logdenL[i+1]-logdenL[i])/(logz[i+1]-logz[i])*(logalt-logz[i]) )
    elif F107 <= 140:
       d0 = 10.**(  logdenL[i]+(logdenL[i+1]-logdenL[i])/(logz[i+1]-logz[i])*(logalt-logz[i]) )
       d1 = 10.**(  logdenM[i]+(logdenM[i+1]-logdenM[i])/(logz[i+1]-logz[i])*(logalt-logz[i]) )
       d = d0 + (d1-d0)*(F107-65.)/75.
    elif F107 <= 250:
       d0 = 10.**(  logdenM[i]+(logdenM[i+1]-logdenM[i])/(logz[i+1]-logz[i])*(logalt-logz[i]) )
       d1 = 10.**(  logdenHL[i]+(logdenHL[i+1]-logdenHL[i])/(logz[i+1]-logz[i])*(logalt-logz[i]) )
       d = d0 + (d1-d0)*(F107-140.)/110.
    else:
       d = 10.**(  logdenHL[i]+(logdenHL[i+1]-logdenHL[i])/(logz[i+1]-logz[i])*(logalt-logz[i]) )
    return d

#
# Initialize the grid using the input file
#
def initSystem(Nshell,LowShell,HiShell,dz,z,finit):

    if finit==None: finit="orbital.distro"

    Sats = []
    DSats= []
    DebL = []
    DebS = []

    Lam = np.zeros(Nshell)

    for i in range(Nshell): 
         Sats.append(NType())
         DSats.append(NType())
         DebL.append(NType())
         DebS.append(NType())

    fh = open(finit,'r')
    SatSize,LargeSize,SmallSize,SatAM,LargeAM,SmallAM=fh.readline().split(",")

    SatSize = float(SatSize)
    LargeSize = float(LargeSize)
    SmallSize = float(SmallSize)
    SatAM = float(SatAM)
    LargeAM = float(LargeAM)
    SmallAM = float(SmallAM)

    for i in range(Nshell):
         Sats[i].setZ(z[i],dz)
         Sats[i].setV()
         Sats[i].AoM = SatAM
         Sats[i].setSize(SatSize)

         DSats[i].setZ(z[i],dz)
         DSats[i].setV()
         DSats[i].AoM = SatAM
         DSats[i].setSize(SatSize)

         DebL[i].setZ(z[i],dz)
         DebL[i].setV()
         DebL[i].setSize(LargeSize)
         DebL[i].AoM = LargeAM

         DebS[i].setZ(z[i],dz)
         DebS[i].setV()
         DebS[i].setSize(SmallSize)
         DebS[i].AoM = SmallAM

    for line in fh:
        alt,S,D,Large,Small,SatLaunch=line.split(",")
        alt = float(alt)
        S   = float(S)
        D   = float(D)
        Large=float(Large)
        Small=float(Small)
        SatLaunch=float(SatLaunch)/yr2sec

        id = int((alt-LowShell)/dz)
        if id < 0 or id > Nshell-1:
            print("initialization altitude out of bounds")
            sys.exit(-1)

        Sats[id].N = S
        DSats[id].N = D
        DebL[id].N = Large
        DebS[id].N = Small
        Lam[id] = SatLaunch

    fh.close()
    return Sats,DSats,DebL,DebS,Lam
        

#
# Driving script for running evolution
#

def simSystem(Nshell=120,LowShell=300e3,HiShell=1500e3,dtyr=0.01,TENDyr=10,alpha1=0.0,alpha2=0.0,alpha3=0.2,beta=1,P=0.95,NewDebL=300,NewDebS=16500,SprayDebS=10,SatLifeTimeYr=5,mo0=0,idSnap=25,setF107=None,INITFILE=None,OUTFILE="multi.out"):
    dz = (HiShell-LowShell)/Nshell
    z  = np.linspace(LowShell,HiShell,Nshell,endpoint=False)+dz*0.5

    Sats,DSats,DebL,DebS,Lam = initSystem(Nshell,LowShell,HiShell,dz,z,INITFILE)

    Sats_p = copy.deepcopy(Sats)
    DSats_p = copy.deepcopy(DSats)
    DebL_p = copy.deepcopy(DebL)
    DebS_p = copy.deepcopy(DebS)

    dt = dtyr * yr2sec
    t = 0
    trecord=0
    TEND = TENDyr * yr2sec

    trecordDelta = TEND/1000

    denatm = np.zeros(Nshell)
    fluxDSats = np.zeros(Nshell+1)
    fluxDebL  = np.zeros(Nshell+1)
    fluxDebS  = np.zeros(Nshell+1)

    tsnap=[]
    SatsSnap=[]
    DSatsSnap=[]
    DebLSnap=[]
    DebSSnap=[]
    DSatsTauSnap=[]

    for i,N in enumerate(DebS): 
       denatm[i] = setDenatm(z[i],t,mo0=0,setF107=65)
       print(N.z,denatm[i],Sats[i].N,DSats[i].N,DebL[i].N,DebS[i].N)

    volume = np.zeros(Nshell)
    for i in range(Nshell): volume[i] = 4*np.pi * ( (RE+z[i]+dz*0.5)**3 - (RE+z[i]-dz*0.5)**3)/3.

    vcolfac = 4./3.
    while t < TEND:

       #STEP 1: PREDICTOR
       for i in range(Nshell):
           denatm[i] = setDenatm(z[i],t,mo0=mo0,setF107=setF107)
           Sats[i].setTlife(SatLifeTimeYr*yr2sec)
           DSats[i].setTau(denatm[i])
           DebL[i].setTau(denatm[i])
           DebS[i].setTau(denatm[i])

           fluxDSats[i] = -DSats[i].N/DSats[i].tau
           fluxDebL[i]  = -DebL[i].N/DebL[i].tau
           fluxDebS[i]  = -DebS[i].N/DebS[i].tau 

       for i in range(Nshell):
           
           vcoll = Sats[i].v * vcolfac/volume[i]

           SatSat  = alpha1 * Sats[i].N * Sats[i].N  * 4 * Sats[i].area* vcoll
           SatDSat = alpha2 * Sats[i].N * DSats[i].N * 4 * Sats[i].area* vcoll
           DSatDSat=         DSats[i].N * DSats[i].N * 4 * Sats[i].area* vcoll
           SatDebL = alpha3 *     Sats[i].N * DebL[i].N  * Sats[i].area* vcoll 
           SatDebS =              Sats[i].N * DebS[i].N  * Sats[i].area* vcoll
           DSatDebL=             DSats[i].N * DebL[i].N  * Sats[i].area* vcoll
           DSatDebS=             DSats[i].N * DebS[i].N  * Sats[i].area* vcoll


           dotDebL = NewDebL*(2*(SatSat+SatDSat+DSatDSat) + SatDebL + DSatDebL) + fluxDebL[i] - fluxDebL[i+1]
           dotDebS = NewDebS*(2*(SatSat+SatDSat+DSatDSat) + SatDebL + DSatDebL) + (DSatDebS+SatDebS)*SprayDebS + fluxDebS[i] - fluxDebS[i+1]

           dotSats = Lam[i] - beta*SatDebS -2*SatSat - SatDSat - SatDebL - Sats[i].N/Sats[i].Tlife
           dotDSats=        + beta*SatDebS -2*DSatDSat - SatDSat - DSatDebL + fluxDSats[i] - fluxDSats[i+1] + (1-P)*Sats[i].N/Sats[i].Tlife
           
           Sats_p[i].N  = Sats[i].N  + dotSats * dt
           DSats_p[i].N = DSats[i].N + dotDSats* dt
           DebL_p[i].N  = DebL[i].N  + dotDebL * dt
           DebS_p[i].N  = DebS[i].N  + dotDebS * dt

       #STEP 2: Corrector
       for i in range(Nshell):
           denatm[i] = setDenatm(z[i],t+dt,mo0=mo0,setF107=setF107)
           Sats_p[i].setTlife(SatLifeTimeYr*yr2sec)
           DSats_p[i].setTau(denatm[i])
           DebL_p[i].setTau(denatm[i])
           DebS_p[i].setTau(denatm[i])

           fluxDSats[i] = -DSats_p[i].N/DSats_p[i].tau
           fluxDebL[i]  = -DebL_p[i].N/DebL_p[i].tau
           fluxDebS[i]  = -DebS_p[i].N/DebS_p[i].tau 

       for i in range(Nshell):
           
           vcoll = Sats[i].v * vcolfac/volume[i]

           SatSat  = alpha1 * Sats_p[i].N * Sats_p[i].N  * 4 * Sats[i].area* vcoll
           SatDSat = alpha2 * Sats_p[i].N * DSats_p[i].N * 4 * Sats[i].area* vcoll
           DSatDSat=         DSats_p[i].N * DSats_p[i].N * 4 * Sats[i].area* vcoll
           SatDebL = alpha3 *     Sats_p[i].N * DebL_p[i].N  * Sats[i].area* vcoll 
           SatDebS =              Sats_p[i].N * DebS_p[i].N  * Sats[i].area* vcoll
           DSatDebL=             DSats_p[i].N * DebL_p[i].N  * Sats[i].area* vcoll
           DSatDebS=             DSats_p[i].N * DebS_p[i].N  * Sats[i].area* vcoll


           dotDebL = NewDebL*(2*(SatSat+SatDSat+DSatDSat) + SatDebL + DSatDebL) + fluxDebL[i] - fluxDebL[i+1]
           dotDebS = NewDebS*(2*(SatSat+SatDSat+DSatDSat) + SatDebL + DSatDebL) + (DSatDebS+SatDebS)*SprayDebS + fluxDebS[i] - fluxDebS[i+1]

           dotSats = Lam[i] - beta*SatDebS -2*SatSat - SatDSat - SatDebL - Sats[i].N/Sats[i].Tlife
           dotDSats=        + beta*SatDebS -2*DSatDSat - SatDSat - DSatDebL + fluxDSats[i] - fluxDSats[i+1] + (1-P)*Sats[i].N/Sats[i].Tlife
           
           Sats[i].N  = (Sats_p[i].N  + Sats[i].N  + dotSats * dt)*0.5
           DSats[i].N = (DSats_p[i].N + DSats[i].N + dotDSats* dt)*0.5
           DebL[i].N  = (DebL_p[i].N  + DebL[i].N  + dotDebL * dt)*0.5
           DebS[i].N  = (DebS_p[i].N  + DebS[i].N  + dotDebS * dt)*0.5

            
       t+=dt 
       if t > trecord:

           tsnap.append(t/yr2sec)
           SatsSnap.append(Sats[idSnap].N)
           DSatsSnap.append(DSats[idSnap].N)
           DebLSnap.append(DebL[idSnap].N)
           DebSSnap.append(DebS[idSnap].N)
           DSatsTauSnap.append(DSats[idSnap].tau/yr2sec) 
           trecord += trecordDelta

    sumSat = 0
    sumDSat = 0
    sumDebL = 0
    sumDebS = 0
    for i,N in enumerate(DebS): 
       sumSat+=Sats[i].N
       sumDSat+=DSats[i].N
       sumDebL+=DebL[i].N
       sumDebS+=DebS[i].N

       denatm[i] = setDenatm(z[i],t,mo0=0,setF107=65)

       print(N.z,denatm[i],Sats[i].N,DSats[i].N,DebL[i].N,DebS[i].N)

    print("Sums: {} Sats, {} DSats, {} DebL, {} DebS".format(sumSat,sumDSat,sumDebL,sumDebS))

    import matplotlib.pylab as plt

    plt.figure()
    plt.plot(tsnap,SatsSnap,label="Sats")
    plt.plot(tsnap,DSatsSnap,label="Derelict")
    plt.plot(tsnap,DebLSnap,label="DebrisL")
    plt.legend()
    plt.yscale('log')
    plt.ylim(1,1e5)
    plt.xlim(0,TENDyr)
    plt.show()

    dfout = open(OUTFILE,"w")
    dfout.write("# time(yr) Sats DSats DebL DebS DSats.tau(yr)\n")
    for i in range(len(tsnap)):
       dfout.write("{} {} {} {} {} {}\n".format(tsnap[i],SatsSnap[i],DSatsSnap[i],DebLSnap[i],DebSSnap[i],DSatsTauSnap[i]))

    dfout.close()
    
