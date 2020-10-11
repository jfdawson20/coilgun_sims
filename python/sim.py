#!/usr/bin/python3
from prettytable import PrettyTable
import math
import matplotlib.pyplot as plt

coil_profiles = []

#copper wire table for calculating resistance based on wire size
copper_res = 1.72 * math.pow(10,-8) #ohm/m
wire_dict = {"18AWG" : 1.02362,
             "20AWG" : 0.8128,
             "22AWG" : 0.64262,
             "24AWG" : 0.51054,
             "26AWG" : 0.40386,
             "28AWG" : 0.32004,
             "30AWG" : 0.254}

#small class/container to track position data / properties
class Projectile():
    def __init__(self,mass):
        self.init           = True
        self.locx           = 0
        self.velocity       = 0
        self.acceleration   = 0
        self.mass           = mass
        self.ur             = 100 #carbon steel
        
        
    def update(self,acceleration,time):
        return 0

#coil class
#represents the physical properties of a single coil 
#user can either pass in the coil configuration via a configuration dictionary (preffered) or for manual testing
#the coil properties can be passed in via separate arguments 
#coils have a specified location on the x-axis as provided by the "locx" config parameter. The "locx" 

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
        #else, parse configuration dictionary 
        else:
            self.parse_config(coilConfig)
        
        #calculate constants or derived values 
        #based on user provided data 

        self.u0             = 4*math.pi*pow(10,-7)
        self.area           = math.pi*pow(self.out_dia/2,2)
        self.wireLen        = (2*math.pi*(self.out_dia/2)) * self.turns 
        self.wireDia        = (wire_dict[self.gauge]/1000)  
        self.wireArea       = math.pi*pow(self.wireDia/2,2)
        self.wireAreaMils   =  (self.wireArea * 1973525241.77) 
        self.seriesR        = copper_res * self.wireLen/self.wireArea
        self.airGap         = self.out_dia - self.core_dia
        self.maxSLT         = math.floor(self.len / self.wireDia)
        self.layers         = math.ceil(self.turns / self.maxSLT)
    
        self.maxOnTime = 0.0297 * math.log10((self.maxt + 234)/(self.ambt + 234)) * (pow((self.wireAreaMils) ,2)/pow(self.maxc,2))

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
            self.ambt       = config["ambt"]
            self.maxt       = config["maxt"]
            self.maxc       = config["maxc"]

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
        x.add_row(["Maximum Turns in a single layer",self.maxSLT,"Turns"])
        x.add_row(["Number of coil layers",self.layers,"layers"])
        x.add_row(["Ambient Coil Temperature",self.ambt,"C"])
        x.add_row(["Max Coil Temperature",self.maxt,"C"])
        x.add_row(["Max Coil Current",self.maxc,"A"])
        x.add_row(["Max On Time",self.maxOnTime,"s"])

         
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


#coil circuit, a combination of a coil, capacitor, switch series resistance
#contains functiions to calculate current, magnetic force, and coil inductance as a function of 
#a provided projectile position. a coil circuit has a location on the x-axis defined by the 
#locx parameter in the associated coilConfig

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
         
        self.cur_v   = self.V0
        self.cur_i   = 0
        self.time    = 0
        
        self.maxCurrent, self.maxTime = self.timeToMaxI()
        print("Max Current %f, Max Time %f\n" % (self.maxCurrent,self.maxTime))

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

    def timeToMaxI(self):
        i = 0
        maxTime = 0
        max_current = 0
        while(i < 0.01):
            current = abs(self.getCurrent(i))
            if(current > max_current):
                max_current = current
                maxTime = i
            i = i +0.000001

        return (max_current,maxTime)

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

    
    #execute simulation with specified parameters
    def Exec(self,dur,prefire_en=False): 

        time=0
        cntr=0
        #arrays of velocity,acceleration,position data (data,timestamp)
        #these are global and are tracked for the full sim 
        accel   = [(0,0)]
        veloc   = [(0,0)]
        posit   = [(0,0)]
        stageResults = []

        #iterate across all stages, calculating the changes in coil and projectile properties over time
        for i in range(self.numStages):
            samples = []
            coil_ovr    = 0
            prefire     = 0
            coilTime    = 0
            #if we are past the first stage , assume projectile is moving and calculate optimum turn on point
            if(prefire_en == True):
                if(i != 0 and (self.proj.locx < self.coils[i].Coil.locx)):
                    distanceToCoil = self.coils[i].Coil.locx - self.proj.locx
                    timeToCoil = distanceToCoil / self.proj.velocity    
                
                    #if the projectile will hit the coil before optimal current is reached, just enable the coil.
                    if(timeToCoil < self.coils[i].maxTime):
                        coil_ovr = 1
                
                    else:
                        prefire = -(self.timeToCoil * self.proj.velocity)
    
            else:
                prefire     = 0
                coil_ovr    = 0


            #stop coil run when projectile reaches midway point in col 
            while(self.proj.locx < (self.coils[i].Coil.locx + self.coils[i].Coil.len/2)):
                tmp = {}
                                    
                #if were not in the coil yet, don't fire   
                if(coil_ovr == 0):
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
                self.proj.locx          = position
                self.proj.velocity      = velocity
                self.proj.acceleration  = acceleration

                #log data 
                #tmp = [current,force,acceleration,velocity,inductance,coilTime,self.proj.locx,time]
                tmp["CURRENT"]          = current
                tmp["FORCE"]            = force
                tmp["ACCELERATION"]     = acceleration
                tmp["VELOCITY"]         = velocity
                tmp["INDUCTANCE"]       = inductance 
                tmp["COIL_TIME"]         = coilTime
                tmp["PROJ_LOCX"]        = self.proj.locx
                tmp["ABS_TIME"]          = time

                samples.append(tmp)
                time = time + 0.000001
                coilTime = coilTime + 0.000001
                cntr = cntr + 1
            

            ret = sim.processStageData(samples,i)
            stageResults.append(ret)
            
            print(self.displayData(ret))
            print("\n\n")
        x = PrettyTable()
        x.field_names = ["Stage","On Time","Max Coil Temp","Peak Current","Exit Velocity"] 
        st = 0
        for s in stageResults:
            
            peakC = 0
            for data in s:
                if (abs(data["CURRENT"]) > peakC):
                    peakC = abs(data["CURRENT"])

            x.add_row([st,s[len(s)-1]["COIL_TIME"],s[len(s)-1]["COIL_TEMP"],peakC,s[len(s)-1]["VELOCITY"]])
            #print(s[len(s)-1])

            st = st + 1
        print(x)
        
        self.plotData(stageResults)

        return 0

    def processStageData(self,data,stage): 
        results = []

        #calculate current integral to back calculate capacitor voltage
        #calculate temperature rise of coil 
        isum    = []
        coil_temp = []
        for i in range(len(data)):
            if (i == 0):
                isum.append(0)
                coil_temp.append(self.coils[stage].Coil.ambt)
                continue 
            
            #calculate current integral
            cur_int = (((data[i-1]["CURRENT"] + data[i]["CURRENT"])/2) * (data[i-1]["COIL_TIME"] - data[i]["COIL_TIME"])) + isum[i-1]
            isum.append(cur_int)
            
            #calculate coil temperature
            temper  = (coil_temp[i-1] + 234) * pow(10,(pow(data[i]["CURRENT"],2)/pow(self.coils[stage].Coil.wireAreaMils,2))*(data[i]["COIL_TIME"] - data[i-1]["COIL_TIME"])/0.0297) - 234
            coil_temp.append(temper)

        #build new results dictionary to return
        for i in range(len(data)):
            tmp = {}
            if(i == 0):
                tmp["CAP_VOLTAGE"]  = self.V0
            else: 
                voltage = self.V0 - (1/self.C)*isum[i]
                tmp["CAP_VOLTAGE"]  = voltage

                
            tmp["COIL_TIME"]    = data[i]["COIL_TIME"]
            tmp["CURRENT"]      = data[i]["CURRENT"]
            tmp["ABS_TIME"]     = data[i]["ABS_TIME"]
            tmp["INDUCTANCE"]   = data[i]["INDUCTANCE"]
            tmp["FORCE"]        = data[i]["FORCE"]
            tmp["ACCELERATION"] = data[i]["ACCELERATION"]
            tmp["VELOCITY"]     = data[i]["VELOCITY"]
            tmp["PROJ_LOCX"]    = data[i]["PROJ_LOCX"]
            tmp["COIL_TEMP"]    = coil_temp[i]
           
            results.append(tmp)
        
        return(results)

    def displayData(self,results):
        x = PrettyTable()
        #x.field_names = ["Cap Voltage", "Circuit Current", "Force", "TIme"]
        x.field_names = ["Absolute Time (s)", "Coil Time (s)","Coil Temp (C)","Cap Voltage (V)", "Circuit Current (A)", "Coil Inductance (L)","Force (N)","Projectile Acceleration (m/s^2)", "Projectile Velocity (m/s)", "Projectile Position (m)"]
        for r in results: 
            x.add_row([r["ABS_TIME"],r["COIL_TIME"],r["COIL_TEMP"],r["CAP_VOLTAGE"],r["CURRENT"],r["INDUCTANCE"],r["FORCE"],r["ACCELERATION"],r["VELOCITY"],r["PROJ_LOCX"]])

        return(x)

    def plotData(self,stageRes):
        
        abstime = []
        current = []
        veloc   = []
        pos     = []
        acc     = []

        fig , axs = plt.subplots(2, 2)
       
        #first graph coil current vs time
        axs[0, 0].set_xlabel('Absolute Time (s)')
        axs[0, 0].set_ylabel('Coil(s) Current (A)')

        for s in stageRes:
            marker = 0 
            for data in s:
                if(data["FORCE"] != 0 and marker == 0):
                    axs[0, 0].axvline(x=data["ABS_TIME"])
                    marker = 1

                abstime.append(data["ABS_TIME"])
                current.append(data["CURRENT"])
                veloc.append(data["VELOCITY"])
                pos.append(data["PROJ_LOCX"])
                acc.append(data["ACCELERATION"])
        
            axs[0, 0].plot(abstime,current)
            axs[1, 0].plot(abstime,veloc)
            axs[0, 1].plot(abstime,pos)
            axs[1, 1].plot(abstime,acc)

    
        axs[1, 0].set_xlabel('Absolute Time (s)')
        axs[1, 0].set_ylabel('Projectile Velocity (m/s)')

        axs[0, 1].set_xlabel('Absolute Time (s)')
        axs[0, 1].set_ylabel('Projectile Position (m)')

        axs[1, 1].set_xlabel('Absolute Time (s)')
        axs[1, 1].set_ylabel('Projectile Acceleration (m/s^2)')


        plt.show()

if __name__ == "__main__":
   

    for i in range(10):
        cfg = {}
        cfg["cid"]      = 0
        cfg["length"]   = 0.0381
        cfg["core_dia"] = 0.006
        cfg["out_dia"]  = 0.009525
        cfg["turns"]    =  100 + (i*50)
        cfg["gauge"]    = "26AWG"
        cfg["locx"]     = 0
        cfg["ambt"]     = 25
        cfg["maxt"]     = 60
        cfg["maxc"]     = 500

        coil_profiles.append(cfg)

    sim = Simulation(200,0.006,coil_profiles[0],8,0.0127)

    sim.Exec(0.001,True)

    
