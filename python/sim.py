#!/usr/bin/python3
from prettytable import PrettyTable
import math
import matplotlib.pyplot as plt
import time as Time
import argparse 

coil_profiles = []



#copper wire table for calculating resistance based on wire size
copper_res = 1.72 * math.pow(10,-8) #ohm/m
wire_dict = {"18AWG" : 1.02362,
             "20AWG" : 0.8128,
             #"20AWG" : 0.85,
             "22AWG" : 0.64262,
             "24AWG" : 0.51054,
             "26AWG" : 0.40386,
             #"26AWG" : 0.43,
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
      
        self.wireDia        = (wire_dict[self.gauge]/1000)  
        self.wireArea       = math.pi*pow(self.wireDia/2,2)
        self.wireAreaMils   =  (self.wireArea * 1973525241.77) 
        self.maxSLT         = math.floor(self.len / self.wireDia) 
        self.layers         = math.ceil(self.turns / self.maxSLT)
 
        self.wireLen        = 0
        layer               = 0
        for i in range(self.turns+1):
            if(i == 0):
                layer = 0
            else:
                layer = math.ceil(i/self.maxSLT) - 1 

            self.wireLen = self.wireLen + (2*math.pi*((self.core_dia + layer*self.wireDia*2)/2))
        
        self.seriesR        = copper_res * self.wireLen/self.wireArea
        
        self.u0             = 4*math.pi*pow(10,-7)
        self.wireDepth      = self.layers * self.wireDia 
        
        #mean outer diameter
        self.out_dia        = self.core_dia + self.wireDepth
        
        self.inner_r        = (self.core_dia)/2
        self.outer_r        = (self.core_dia + self.layers*self.wireDia*2)/2


        #mean coil area used for inductance calcs
        self.area           = math.pi*pow(self.core_dia/2,2)

        #projectile area used for force calcs 
        self.forceArea      = math.pi*pow(self.proj_dia/2,2)

       
        #self.airGap         = self.out_dia - self.core_dia
        self.airGap         = self.core_dia - self.proj_dia

   
        self.maxOnTime = 0.0297 * math.log10((self.maxt + 234)/(self.ambt + 234)) * (pow((self.wireAreaMils) ,2)/pow(self.maxc,2))

        self.dumpStats()
    
    def parse_config(self,config):
            self.init       = True
            self.cid        = config["cid"];
            self.len        = config["length"]
            self.core_dia   = config["core_dia"]
            #self.out_dia    = config["out_dia"]
            self.turns      = config["turns"]
            self.locx       = config["locx"]
            self.gauge      = config["gauge"]
            self.ambt       = config["ambt"]
            self.maxt       = config["maxt"]
            self.maxc       = config["maxc"]
            self.profile    = config["profile"]
            self.proj_dia   = config["proj_dia"]

    def calcSeriesR(self):
        return self.seriesR
        
    def dumpStats(self):
        ind = self.calcInductance(0,100)
        print(ind)

        x = PrettyTable()
        x.field_names = ["Parameter", "Value", "Description"]
        x.add_row(["Coil ID",self.cid,"Coil Identifier"])
        x.add_row(["Air Core Inductance",ind,"Coil Inductance (air-core) (H)"])
        x.add_row(["Coil Length",self.len,"Coil Length (m)"])
        x.add_row(["Coil Inner Diameter",self.core_dia,"Barrel outer Diameter (m)"])
        x.add_row(["Projectile Diameter",self.proj_dia,"Projectile Diameter (m)"])
        x.add_row(["Coil Outer Diameter",self.out_dia,"Mean diameter of multilayer coil (m)"])
        x.add_row(["Coil Inner Radius",self.inner_r,"inner radius of coil (m)"])
        x.add_row(["Coil Outer Radius",self.outer_r,"outer radius of coil (m)"])
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
        #old single layer calc
        #L0      = (ur*self.u0*pow(self.turns,2)*self.area)/self.len
        
        #l(uH)  = 31.6*N
         
        L1      = (31.6 * pow(self.turns,2) * pow(self.out_dia/2,2))/(6*(self.out_dia/2) + 9*self.len + 10*(self.wireDepth)) 
        L1      = L1 / pow(10,6)

        L0      = L1 * ur
        alpha   = math.log(ur)
        
        Lx  = L0 *  math.exp(-(alpha/self.len)*(normx))
        
        return(Lx) 

        #return 0


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
         
        #self.cur_v   = self.V0
        #self.cur_i   = 0
        #self.time    = 0
        
        self.maxCurrent, self.maxTime = self.timeToMaxI()
        #print("Max Current %f, Max Time %f\n" % (self.maxCurrent,self.maxTime))
    
    def dumpCircuitConfig(self):
        print("Voltage: %f\n" % self.V0)
        print("Capacitance: %f\n" % self.C)
        self.Coil.dumpStats()

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

    def timeToTargetI(self,target):
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
        F = (pow(self.Coil.turns*i,2)*self.Coil.u0*self.Coil.forceArea)/(2*pow(self.Coil.airGap,2))

        if(core_locx < self.Coil.locx or ((core_locx - self.Coil.locx) > self.Coil.len/2)):
            return 0#1/pow((self.Coil.locx - core_locx),2)

            
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
        self.switchR = 0.00014
        self.numStages = stages
        self.spacing   = stageSpacing 
        self.proj = Projectile(0.00569)
        
        self.coils = []
        #create n stages of coils 
        for i in range(self.numStages):
            tmp = cfg
            tmp["cid"] = i
            tmp["locx"] = i*(self.spacing + tmp["length"])
            tmp["turns"] = cfg["turns"]-20

            self.coils.append(CoilCircuit(i,self.V0,self.C,self.switchR,self.proj,tmp))

    
    def updateCoil(self,coilID, config):
        tmp = config
        tmp["cid"] = coilID
        tmp["locx"] = coilID*(self.spacing + tmp["length"])

        self.coils[coilID] = CoilCircuit(coilID,self.V0,self.C,self.switchR,self.proj,tmp)

        return 0 
    
    def getCoilConfig(self):
        cfgs = []
        for coil in self.coils:
            cfgs.append((coil.Coil.cid,coil.Coil.profile))

        return cfgs

    def dumpConfig(self):
        for coils in self.coils:
            coils.dumpCircuitConfig()

    #execute simulation with specified parameters
    def Exec(self,prefire_en=False): 

        time=0
        cntr=0
        #arrays of velocity,acceleration,position data (data,timestamp)
        #these are global and are tracked for the full sim 
        accel   = [(0,0)]
        veloc   = [(0,0)]
        posit   = [(0,0)]
        stageResults = []
        self.proj = Projectile(0.00569)
        print(self.proj.locx)

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
                        prefire = (timeToCoil * self.proj.velocity)
    
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
            
        return stageResults


    #execute simulation with specified parameters
    def Exec2(self,prefire_en=False): 

        time=0
        cntr=0
        #arrays of velocity,acceleration,position data (data,timestamp)
        #these are global and are tracked for the full sim 
        accel   = [(0,0)]
        veloc   = [(0,0)]
        posit   = [(0,0)]
        stageResults = []
        coilTimes    = []
        distanceToCoils = []
        timeToCoils     = []
        coil_ovrs       = []
        currents        = []
        forces          = []
        inductances     = []
        samples         = []
        simSamples      = []
        stageActive     = []

        for i in range(self.numStages):
            distanceToCoils.append(0)
            timeToCoils.append(0)
            coil_ovrs.append(0)
            coilTimes.append(0)
            samples.append([])
            currents.append(0)
            forces.append(0)
            inductances.append(0)
            stageActive.append(0)
        

        self.proj = Projectile(0.00569)
        #print(self.proj.locx)

        
        while(self.proj.locx < 0.5):
            activeForce = 0
            tmp = {}
            #get the status of all coils in the system
            for i in range(self.numStages):

                if(self.proj.velocity == 0): 
                    distanceToCoils[i] = self.coils[i].Coil.locx 
                    timeToCoils[i]     = 0
                else:
                    distanceToCoils[i]   = self.coils[i].Coil.locx - self.proj.locx
                    timeToCoils[i]      = distanceToCoils[i]/self.proj.velocity
                
              
                #1) get coil updates
                currents[i],forces[i],inductances[i] = self.coils[i].updateCoilCircuit(coilTimes[i],self.proj) 
                #print(i,currents[i],forces[i],inductances[i])

                #only one coil will be applying force at a time
                activeForce = activeForce + forces[i]
                #print(forces[i]) 
            
            #2) Update Projectile properties
            acceleration = activeForce/self.proj.mass
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
            samples = []
            for i in range(self.numStages):
                tmp = {}
                #print(currents[i])
                if(stageActive[i] == 1):
                    tmp["CURRENT"]          = currents[i]
                    tmp["FORCE"]            = forces[i]
                    tmp["ACCELERATION"]     = acceleration
                    tmp["VELOCITY"]         = velocity
                    tmp["INDUCTANCE"]       = inductances[i] 
                    tmp["COIL_TIME"]        = coilTimes[i]
                    tmp["PROJ_LOCX"]        = self.proj.locx
                    tmp["ABS_TIME"]         = time
                else:
                    tmp["CURRENT"]          = currents[i]
                    tmp["FORCE"]            = 0
                    tmp["ACCELERATION"]     = 0
                    tmp["VELOCITY"]         = 0
                    tmp["INDUCTANCE"]       = inductances[i]
                    tmp["COIL_TIME"]        = coilTimes[i]
                    tmp["PROJ_LOCX"]        = self.proj.locx
                    tmp["ABS_TIME"]         = time
 

                samples.append(tmp)

            #update counters 
            for i in range(self.numStages):
                #if we are within prefire range (e.g. projectile will reach coil at or before optimal current peak) and the projectile has also not left the coil
                #enable the coil
                
                if(prefire_en == True):
                    if(timeToCoils[i] <= self.coils[i].maxTime and (self.proj.locx - self.coils[i].Coil.locx <= self.coils[i].Coil.len/2)):
                        coilTimes[i] = coilTimes[i] + 0.000001
                
                    #else set the time to zero which effectivly shuts off the coil
                    else:
                        coilTimes[i] = 0
                else:
                    if(self.proj.locx >= self.coils[i].Coil.locx and  self.proj.locx - self.coils[i].Coil.locx <= self.coils[i].Coil.len/2):
                        coilTimes[i] = coilTimes[i] + 0.000001
                    else:
                        coilTimes[i] = 0
 
                if(self.proj.locx >= self.coils[i].Coil.locx and  self.proj.locx - self.coils[i].Coil.locx <= self.coils[i].Coil.len/2):
                    stageActive[i] = 1
                else:
                    stageActive[i] = 0
            
            #for x in samples:
            #    print(x)
            #    print("\n")

            simSamples.append(samples)

            #Time.sleep(0.01) 
            time = time + 0.000001
            cntr = cntr + 1
            
        
        for i in range(self.numStages):
            sam = []
            for s in simSamples:
                sam.append(s[i])
            
            ret = sim.processStageData(sam,i)
            stageResults.append(ret)

        #for i in range(self.numStages):
        #    ret = sim.processStageData(samples[i],i)
        #    stageResults.append(ret)
            
        return stageResults


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
        
        abstimes = []
        currents = []
        forces   = []
        velocs   = []
        poss     = []
        accs     = []
       
        fig , axs = plt.subplots(3, 2)
         
        #first graph coil current vs time
        axs[0, 0].set_xlabel('Absolute Time (s)')
        axs[0, 0].set_ylabel('Coil(s) Current (A)')
        i = 0
        #for each stage of data
        for s in stageRes:
            abstime = []
            current = []
            force   = []
            veloc   = []
            pos     = []
            acc     = []
 
            #iterate through all datapoints
            marker = 0 
            for data in s:
                #mark entry of projectile into coil 
                if(data["FORCE"] != 0 and marker == 0):
                    axs[0, 0].axvline(x=data["ABS_TIME"])
                    axs[2, 0].axvline(x=data["ABS_TIME"])
                    marker = 1

                current.append(data["CURRENT"])
                veloc.append(data["VELOCITY"])
                acc.append(data["ACCELERATION"])
                force.append(data["FORCE"])
                
                pos.append(data["PROJ_LOCX"])
                abstime.append(data["ABS_TIME"])
 
            abstimes.append(abstime)
            currents.append(current)
            forces.append(force)
            poss.append(pos)
            accs.append(acc)
            velocs.append(veloc)
        
        #print(abstimes[0])

        for i in range(self.numStages):
            axs[0, 0].plot(abstimes[0],currents[i])
            axs[1, 0].plot(abstimes[0],velocs[i])
            axs[1, 1].plot(abstimes[0],accs[i])
            axs[2, 0].plot(abstimes[0],forces[i])

        axs[0, 1].plot(abstimes[0],poss[0])

        axs[1, 0].set_xlabel('Absolute Time (s)')
        axs[1, 0].set_ylabel('Projectile Velocity (m/s)')

        axs[0, 1].set_xlabel('Absolute Time (s)')
        axs[0, 1].set_ylabel('Projectile Position (m)')

        axs[1, 1].set_xlabel('Absolute Time (s)')
        axs[1, 1].set_ylabel('Projectile Acceleration (m/s^2)')

        axs[2, 0].set_xlabel('Absolute Time (s)')
        axs[2, 0].set_ylabel('Coil(s) Force (N)')

           
        plt.show()

    def showResults(self,res,verbose=False,plots=False):
        
        if(verbose == True):
            for ret in res:
                print(self.displayData(ret))
                print("\n\n")
 
        x = PrettyTable()
        x.field_names = ["Stage","On Time","Max Coil Temp","Exit Cap Voltage","Peak Current","Min Inductance","Exit Velocity"] 
        st = 0
        for s in res:
            
            peakC = 0
            peakV = 0 
            peakT = 0
            minV  = self.V0
            for data in s:
                if (abs(data["CURRENT"]) > peakC):
                    peakC = abs(data["CURRENT"])

                if(data["VELOCITY"] > peakV):
                    peakV = data["VELOCITY"]

                if(data["COIL_TEMP"] > peakT):
                    peakT = data["COIL_TEMP"]

                if(data["CAP_VOLTAGE"] < minV):
                    minV = data["CAP_VOLTAGE"]


            x.add_row([st,s[len(s)-1]["COIL_TIME"],peakT,minV,peakC,s[0]["INDUCTANCE"],peakV])
            #print(s[len(s)-1])

            st = st + 1
        print(x)
       
        if(plots == True):
            self.plotData(res)


    def Optimize(self): 
         
        trial_results = []
        run = self.Exec(True)
        trial_results.append(run)
        
        max_veloc = run[self.numStages-1][len(run[self.numStages-1])-1]["VELOCITY"]
        max_cfg = self.getCoilConfig()
        max_stage_cfg = coil_profiles[0]
    
        for stages in range(2):
            for coils in coil_profiles:
                print("Optimizing Stage %i, Coil Profile %f\n" % (stages,coils["profile"]))
                for cfgs in self.getCoilConfig():
                    print(cfgs)

                #self.dumpConfig()

                #update coil configuration
                self.updateCoil(stages,coils)

                #run the sim 
                run = self.Exec(True)
                trial_results.append(run)
                
                #check results. If final stage velocity is higher than last recorded max. update 
                if(run[self.numStages-1][len(run[self.numStages-1])-1]["VELOCITY"] > max_veloc):
                    max_veloc       = run[len(run)-1]["VELOCITY"]
                    max_cfg         = getCoilConfig()
                    max_stage_cfg   = coils

                self.showResults(run,verbose=False)

            self.updateCoil(stages,max_stage_cfg)

        print("RESULTS\n")
        print("Max Velocity: %f\n" % (max_veloc))
        for c in max_cfg: 
            print(c)

        return 0

    def coilCalculator(self,targetR, corelen, core_dia):
        for gauges in wire_dict:
            wireDia        = (wire_dict[gauges]/1000)  
            wireArea       = math.pi*pow(wireDia/2,2)
            wireAreaMils   =  (wireArea * 1973525241.77) 
            maxSLT         = math.floor(corelen / wireDia) 
        
            wireLen = (targetR * wireArea)/copper_res
            
            turns = wireLen / (2*math.pi*(core_dia/2))
            print("%s : %f" % (gauges,turns))
 
    #dynamically build a list of coil configs that fall within the range
    def buildCoilConfigs(self,min_i, max_i, length, shaft_dia, max_coil_dia, proj_dia, num):


        return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Multi-Stage Coilgun simulator/optimizer")
    parser.add_argument('-o','--op',dest='simOp',help="Simulation Operation to Perform")
    parser.add_argument('-p','--pre',dest='preShoot',help="When simulating, pre-activate each coil before projectile arrival: AUTO=auto calculate optimal firing time, DISABLE = do not use preshoot activation",default="AUTO")



    args = parser.parse_args()

    if(args.simOp == "exec"):
        preshoot_en = False
        if(args.preShoot == "AUTO"):
            preshoot_en = True


        for i in range(10):
            cfg = {}
            cfg["profile"]  = i
            cfg["cid"]      = 0
            cfg["length"]   = 0.035
            cfg["proj_dia"] = 0.006
            cfg["core_dia"] = 0.0098
            cfg["turns"]    = 350#80 + (i*10)
            cfg["gauge"]    = "22AWG"
            cfg["locx"]     = 0
            cfg["ambt"]     = 25
            cfg["maxt"]     = 60
            cfg["maxc"]     = 500

            coil_profiles.append(cfg)

        sim = Simulation(200,0.006,coil_profiles[0],8,0.0127)

        #ret = sim.Optimize() 

        ret = sim.Exec2(preshoot_en)
        sim.showResults(ret,plots=True,verbose=False)

    elif(args.simOp == "calc"):
        sim.coilCalculator(0.5,0.035,0.0098)
