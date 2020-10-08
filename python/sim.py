#!/usr/bin/python3
from prettytable import PrettyTable
import math




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

    def updatePosition(self,locx):
        self.locx = locx
        return 0 

    def updateSpeed(self,v):
        self.v = v  
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
            normx = 1;

        #projectile inside core 
        elif(core_locx - self.locx <= self.len):
            normx = self.len - (core_locx - self.locx)
        else:
            normx = 1; 

        
        L0      = (ur*self.u0*pow(self.turns,2)*self.area)/self.len
        L1      = L0 / ur
        alpha   = math.log(ur)
        
        if(normx == 1):
            return L1
        else: 
            Lx  = L0 * math.exp(-(alpha/self.len)*(core_locx))
            return(Lx) 

        return 0


class CoilCircuit():
    def __init__(self,V0,Cap,switchR,proj,coilConfig):
        self.init = True 
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
            it = ((-self.V0 * math.exp(-self.Aplha * time))/(self.L * self.Wd))*math.sin(self.Wd * time)

        #critically damped condition
        else:
            it  = (-self.V0*time*math.exp(-self.Alpha * time))/self.L

        return it

    def getForce(self,i):
        F = (pow(self.Coil.turns*i,2)*self.Coil.u0*self.Coil.area)/(2*pow(self.Coil.airGap,2))
        return F

    def updateCoilCircuit(self,time,proj)
        #update inductance with position of projectile 
        self.L       = self.Coil.calcInductance(proj.locx,proj.ur)

        #get current 
        i = self.getCurrent(time)

        #get force 
        f = self.getForce(i)

    def processData(self,data):

        results = []
        isum    = []
        for i in range(len(data)):
            if (i == 0):
                isum.append(0)
                continue 
            
            tmp = (((data[i-1][0] + data[i][0])/2) * (data[i-1][2] - data[i][2])) + isum[i-1]
            #print(tmp)
            isum.append(tmp)

        for i in range(len(data)):
            if(i == 0):
                tmp = (self.V0,data[0][0],data[0][1],data[0][2])
            
            else: 
                voltage = self.V0 - (1/self.C)*isum[i]
                tmp = (voltage,data[i][0],data[i][1],data[i][2])
            
            results.append(tmp)
        
        x = PrettyTable()
        x.field_names = ["Cap Voltage", "Circuit Current", "Force", "TIme"]
        for r in results: 
            x.add_row([r[0],r[1],r[2],r[3]])

        return(x)
                


    
    def StepCircuit(self,proj):
        return 0;



class Coilgun():
    def __init__(self):
        self.init = True


class Simulation():
    def __init__(self):
        self.init = True

    def Exec(self): 
        print("hello")
        return 0

if __name__ == "__main__":
    sim = Simulation()
    sim.Exec()

    proj = Projectile(0.007)
    
    cfg = {}
    cfg["cid"]      = 0
    cfg["length"]   = 0.0254
    cfg["core_dia"] = 0.006
    cfg["out_dia"]  = 0.009525
    cfg["turns"]    = 200
    cfg["gauge"]    = "26AWG"
    cfg["locx"]     = 0.0254

    #c = Coil(False,0,0.0254,0.006,0.009525,200,"20AWG",0.0254)
    c = CoilCircuit(200,0.006,0.00014,proj,cfg)
   
    i=0
    samples = []
    while (i < 0.00050):
        tmp = c.getForce(i,proj)
        samples.append(tmp)
        #print(tmp)
        i = i + 0.000001
    
    
    ret = c.processData(samples)

    print(ret)
