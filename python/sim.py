#!/usr/bin/python3
from prettytable import PrettyTable
import math
import matplotlib.pyplot as plt



copper_res = 1.72 * math.pow(10,-8) #ohm/m

wire_dict = {"18AWG" : 1.02362,
             "20AWG" : 0.8128,
             "22AWG" : 0.64262,
             "24AWG" : 0.51054,
             "26AWG" : 0.40386,
             "28AWG" : 0.32004,
             "30AWG" : 0.254}

class Projectile():
    def __init__(self,mass):
        self.init = True
        self.locx   = 0
        self.v      = 0
        self.mass   = mass
        self.ur     = 4000
        
    def update(self,acceleration,time):
        return 0

class Coil():
    def __init__(self,coilConfig,cid=False,length=False,core_dia=False,out_dia=False,turns=False,gauge=False,locx=False):
        
        #if config  dictionary is not present 
        if(coilConfig == False):
            self.init       = True
            self.cid        = cid;
            self.len        = length
            self.core_dia   = core_dia
            self.out_dia    = out_dia
            self.turns      = turns 
            self.locx       = locx
            self.gauge      = gauge
        else:
            self.parse_config(coilConfig)
        
        #constants or derived values 

        self.u0         = 4*math.pi*pow(10,-7)
        self.area       = math.pi*pow(self.out_dia/2,2)
        self.wireLen    = (2*math.pi*(self.out_dia/2)) * self.turns 
        self.wireDia    = (wire_dict[self.gauge]/1000)  
        self.wireArea   = math.pi*pow(self.wireDia/2,2)
        self.seriesR    = copper_res * self.wireLen/self.wireArea
        self.airGap     = self.out_dia - self.core_dia

        self.dumpStats()
    
    def parse_config(self,config):
            self.init       = True
            self.cid        = config["cid"];
            self.len        = config["length"]
            self.core_dia   = config["core_dia"]
            self.out_dia    = config["out_dia"]
            self.turns      = config["turns"]
            self.locx       = config["locx"]
            self.gauge      = config["gauge"]

    def calcSeriesR(self):
        return self.seriesR
        
    def dumpStats(self):
        x = PrettyTable()
        x.field_names = ["Parameter", "Value", "Description"]
        x.add_row(["Coil ID",self.cid,"Coil Identifier"])
        x.add_row(["Coil Length",self.len,"Coil Length (m)"])
        x.add_row(["Coil Inner Diameter",self.core_dia,"Coil Inner Diameter (Projectile) (m)"])
        x.add_row(["Coil Outer Diameter",self.out_dia,"Coil Outer Diameter (Barrel) (m)"])
        x.add_row(["Coil Cross-sectional Area",self.area,"Coil Cross-sectional Area (Barrel) (m)"])
        x.add_row(["Coil Air Gap",self.airGap,"Coil Air Gap"])
        x.add_row(["Coil turns",self.turns,"Number of Wire Turns in Coil"])
        x.add_row(["Coil wire length",self.wireLen,"Length of wire in coil (m)"])
        x.add_row(["Coil wire gauge",self.gauge,"wire gauge"])
        x.add_row(["Coil wire diameter",self.wireDia,"wire diameter (m)"])
        x.add_row(["Coil wire area",self.wireArea,"wire area (m^2)"])
        x.add_row(["Coil Series Resistance",self.seriesR,"Series Resistance of Coil (Ohm)"])
        
        x.add_row(["Coil x-axis location",self.locx,"start of coil location in sim"])
        
        x.add_row(["Vaccum Permeability Constant",self.u0,"H/m"])
         
        print("\n")
        print(x)
        print("\n")

    def calcInductance(self,core_locx,ur):

        #projectile has not entered core
        if(self.locx >= core_locx):
            normx = self.len;

        #projectile inside core 
        elif(core_locx - self.locx <= self.len):
            normx = self.len - (core_locx - self.locx)
        else:
            normx = self.len; 

        #print(normx)        
        L0      = (ur*self.u0*pow(self.turns,2)*self.area)/self.len
        L1      = L0 / ur
        alpha   = math.log(ur)
        
        if(normx == 0):
            return L0
        else: 
            Lx  = L0 * math.exp(-(alpha/self.len)*(normx))
            return(Lx) 

        return 0


class CoilCircuit():
    def __init__(self,stage,V0,Cap,switchR,proj,coilConfig):
        self.init = True 
        self.id   = stage
        self.on   = False

        self.V0      = V0
        self.C       = Cap
        self.SwitchR = switchR
        
        #derived values / initial conditions
        self.Coil    = Coil(coilConfig) 
        self.L       = self.Coil.calcInductance(0,proj.ur)
        self.R       = self.SwitchR + self.Coil.calcSeriesR()
        print(self.L) 
        self.cur_v   = self.V0
        self.cur_i   = 0
        self.time    = 0
        

    def getCurrent(self,time):
 
        self.Alpha   = (self.R/(2*self.L))
        self.W       = 1/(math.sqrt(self.L * self.C))

        #over-damped condition
        if(self.Alpha > self.W):
            self.S1      = -self.Alpha + math.sqrt(pow(self.Alpha,2) - pow(self.W,2))
            self.S2      = -self.Alpha - math.sqrt(pow(self.Alpha,2) - pow(self.W,2))

            it = ((self.V0 * math.exp(-self.Alpha * time))/(2*self.L*math.sqrt(pow(self.Alpha,2) - pow(self.W,2))))*(math.exp(-math.sqrt(pow(self.Alpha,2) - pow(self.W,2))*time) - math.exp(math.sqrt(pow(self.Alpha,2) - pow(self.W,2))*time ))

        #under-damped condition
        elif(self.Alpha < self.W):
            self.Wd = math.sqrt(pow(self.W,2) - pow(self.Alpha,2))
            it = ((-self.V0 * math.exp(-self.Alpha * time))/(self.L * self.Wd))*math.sin(self.Wd * time)

        #critically damped condition
        else:
            it  = (-self.V0*time*math.exp(-self.Alpha * time))/self.L

        return it

    def getForce(self,i,core_locx):
        #if projectile is not in coil, assume no force applied
        if(core_locx < self.Coil.locx):
            return 0 

        F = (pow(self.Coil.turns*i,2)*self.Coil.u0*self.Coil.area)/(2*pow(self.Coil.airGap,2))
        return F

    def updateCoilCircuit(self,time,proj):
        #update inductance with position of projectile 
        #print(proj.locx)
        self.L       = self.Coil.calcInductance(proj.locx,proj.ur)

        #get current 
        i = self.getCurrent(time)

        #get force 
        f = self.getForce(i,proj.locx)
    
        return(i,f,self.L)

                   
    def StepCircuit(self,proj):
        return 0;
   

class Simulation():
    def __init__(self,V0,C,cfg,stages, stageSpacing):
        self.init = True
        self.V0   = V0
        self.C    = C
        self.numStages = stages
        self.spacing   = stageSpacing 
        self.proj = Projectile(0.00569)
   
        self.coils = []
        #create n stages of coils 
        for i in range(self.numStages):
            tmp = cfg
            tmp["cid"] = i
            tmp["locx"] = i*(self.spacing + tmp["length"])
            self.coils.append(CoilCircuit(i,200,0.006,0.00014,self.proj,tmp))

    def Exec(self,dur): 

        time=0
        cntr=0
        accel   = [(0,0)]
        veloc   = [(0,0)]
        posit   = [(0,0)]
        stageResults = []
        coilTime = 0
        prefire = 0.02
        
        #iterate across all stages, calculating the changes in coil and projectile properties over time
        for i in range(self.numStages):
            samples = []

            #stop coil run when projectile reaches midway point in col 
            while(self.proj.locx < (self.coils[i].Coil.locx + self.coils[i].Coil.len/2)):
                
                #if were not in the coil yet, don't fire
                if(self.proj.locx < (self.coils[i].Coil.locx - prefire)):
                    coilTime = 0
                
                #1) get coil updates
                current,force,inductance = self.coils[i].updateCoilCircuit(coilTime,self.proj) 
        
                #2) Update Projectile properties
                acceleration = force/self.proj.mass
                accel.append((acceleration,time))
        
                if(time == 0):
                    velocity = 0
                    position = 0
                else:
                    velocity = ((accel[cntr-1][0] + accel[cntr][0])/2) * (accel[cntr][1] - accel[cntr-1][1]) + veloc[cntr-1][0]
                    position = ((veloc[cntr-1][0] + veloc[cntr][0])/2) * (veloc[cntr][1] - veloc[cntr-1][1]) + posit[cntr-1][0]
            

                veloc.append((velocity,time))
                posit.append((position,time))
                self.proj.locx = position

                #log data 
                tmp = [current,force,acceleration,velocity,inductance,coilTime,self.proj.locx,time]
        
                samples.append(tmp)
                time = time + 0.000001
                coilTime = coilTime + 0.000001
                cntr = cntr + 1
            

            ret = sim.processStageData(samples)
            stageResults.append(ret)
            
            print(self.displayData(ret))

        x = PrettyTable()
        x.field_names = ["Stage","On Time","Peak Current","Exit Velocity"] 
        st = 0
        for s in stageResults:
            
            peakC = 0
            for data in s:
                if (abs(data[2]) > peakC):
                    peakC = abs(data[2])

            x.add_row([st,s[len(s)-1][0],peakC,s[len(s)-1][6]])
            #print(s[len(s)-1])

            st = st + 1
        print(x)

        return 0

    def processStageData(self,data):
        #tmp = [current,force,acceleration,velocity,inductance,time,proj.locx]
 
        results = []
        isum    = []
        for i in range(len(data)):
            if (i == 0):
                isum.append(0)
                continue 
            
            tmp = (((data[i-1][0] + data[i][0])/2) * (data[i-1][5] - data[i][5])) + isum[i-1]
            #print(tmp)
            isum.append(tmp)

        for i in range(len(data)):
            if(i == 0):
                tmp = (data[0][5],self.V0,data[0][0],data[0][4],data[0][1],data[0][2],data[0][3],data[0][6],data[0][7])
            
            else: 
                voltage = self.V0 - (1/self.C)*isum[i]
                tmp = (data[i][5],voltage,data[i][0],data[i][4],data[i][1],data[i][2],data[i][3],data[i][6],data[i][7])
            
            results.append(tmp)
        
        return(results)

    def displayData(self,results):
        x = PrettyTable()
        #x.field_names = ["Cap Voltage", "Circuit Current", "Force", "TIme"]
        x.field_names = ["Absolute Time", "Coil Time","Cap Voltage", "Circuit Current", "Coil Inductance","Force","Projectile Acceleration", "Projectile Velocity", "Projectile Position"]
        for r in results: 
            #x.add_row([r[0],r[1],r[2],r[3]])
            x.add_row([r[8],r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7]])

        return(x)

if __name__ == "__main__":
   
    cfg = {}
    cfg["cid"]      = 0
    cfg["length"]   = 0.0254
    cfg["core_dia"] = 0.006
    cfg["out_dia"]  = 0.009525
    cfg["turns"]    =  80
    cfg["gauge"]    = "26AWG"
    cfg["locx"]     = 0

    sim = Simulation(200,0.006,cfg,16,0.0508)

    sim.Exec(0.001) 

    
