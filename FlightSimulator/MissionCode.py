from fsHeader import *

def MissionCode(files):

    main_simdf = pd.DataFrame(columns=["Time (sec)","Dist (nM)","Wt (lbm)","MACH","Alt (ft)","Q (psf)","KTAS","KCASideal (conflate with KIAS)","AOA","CL","CD","L/D","M(L/D)","Drag (lbf)","Thrust (lbf)","Lin.Acc (ft/s^2)","ROC (ft/min)","percent induced","FF (lb/hr)","gamma (deg)","SR (nM/lbm)","KEAS (nm/hr)","HEADING (deg)","WINDS (nm/hr)","GROUND SPEED (nm/hr)","ESAD (nM)","MODE","PLA","TSFC","STATIC TEMP (R)","STAG TEMP (R)","RAD EQUI TEMP (R)","REF. HEAT. RATE (BTU/FT^2/SEC)","UNIT RE # (/FT)","Nz"])
    m = MissionClass(files)
    
    # variable inits
    t = 0

    Alt = 0
    Mach = 0
    minPLA = 0
    maxPLA = 0
    dist = 0


    flag = 0
    flag2 = 0
    flag3 = 0
    creditflag = 0
    starttime = 0
    wt = 100
    WF = 50
    itr = True # Iteration boolean

    # Begin Primary Loop
    with open(files.mission_file,'r') as f:
        while line:=f.readline().upper():
            if not itr:
                break
            if wt< WF:
                itr = False

            # Search for comments (leading *)
            if line.startswith("*"):
                continue

            if line.startswith("STOP"): # Not sure if continue is the correct line here
                continue
            if "SET" in line:
                flag = 0
            if "LOG" in line:
                print(t)
                print(wt)
                print(Alt)
                print((Mach*1000)/1000)
                print(line)
                continue
            
            ## Begin of main instruction set
            if "CLIMB" in line:
                flag = 1
                PLA = maxPLA
                creditflag = 1
                starttime = t
                startdist = dist
            if "CLIMB_CONST_MACH" in line:
                flag = 1
                creditflag = 1
                starttime = t
                startdist = dist
            if "CLIMB_CONST_KIAS" in line:
                flag = 2
                creditflag = 1
                starttime = t
                startdist = dist
            if "DESCEND" in line:
                flag = 1
                PLA = minPLA
                creditflag = 1
                starttime = t
                startdist = dist
            if "DESCEND_CONST_MACH" in line:
                flag = 1
                PLA = minPLA
                creditflag = 1
                starttime = t
                startdist = dist
            if "DESCEND_CONST_KIAS" in line:
                flag = 2
                PLA = minPLA
                creditflag = 1
                starttime = t
                startdist = dist
            if "ACCEL" in line:
                flag = 3
                PLA = maxPLA
                creditflag = 1
                starttime = t
                startdist = dist
            if "DECEL" in line:
                flag = 3
                PLA = minPLA
                creditflag = 1
                starttime = t
                startdist = dist
            if "LEVEL" in line:
                flag = 4
                leftright = 0
                creditflag = 1
                starttime = t
                startdist = dist
            if "CONST_CL" in line:
                flag = 5
                leftright = 0
                creditflag = 1
                starttime = t
                startdist = dist
                CLtarget = 0
            if "CONST_AOA" in line:
                flag = 6
                leftright = 0
                creditflag = 1
                starttime = t
                startdist = dist
                AoAtarget = 0
            if "GROUND_RUNUP" in line:
                flag = 7
                PLA = maxPLA
                creditflag = 1
                starttime = t
                startdist = dist
            if flag == 0:
                if "NENG" in line:
                    neng = float(line[len("NENG"):])
                    if neng < 0:
                        neng = 1
                if "W_START" in line:
                    Wi = float(line[len("W_START"):])
                    wt = Wi
                if "W_END" in line:
                    Wf = float(line[len("W_END"):])
                if "WEIGHT" in line:
                    wt = float(line[len("WEIGHT"):])
                if "MACH" in line:
                    Mach = float(line[len("MACH"):])
                    machold = Mach
                if "APU_BURN" in line:
                    APU_burn = float(line[len("APU_BURN"):])
                if "ALTITUDE" in line:
                    Alt = float(line[len("ALTITUDE"):])
                if "DROP" in line:
                    wt -= float(line[len("DROP"):])
                if "PLA" in line:
                    PLA = float(line[len("PLA"):])
                if "LOAD_FACTOR" in line:
                    Nz = float(line[len("LOAD_FACTOR"):])
                if "WINDS" in line:
                    winds = float(line[len("WINDS"):]) * 1.689
                if "w" in line:
                    gross2net = float(line[len("w"):]) / 100
                if "DT" in line:
                    dt = float(line[len("DT"):])
                    if dt < 0.001:
                        dt = 0.001
                    elif dt > 60:
                        dt = 60
                if "DELTA_CD" in line:
                    cdincrement = float(line[len("DELTA_CD"):])
                    optflag = 0
                    # Debug Print
                    print("CD_INCREMENT = {}".format(cdincrement))

                # Input of wind file
                if "WIND_FILE" in line:
                    wind_file = line[len("WIND_FILE"):].strip()
                    files.wind_file = wind_file
                    m.read_wind_data(wind_file)
                if "ZERO_WIND_DATA" in line:
                    wind_file = ""
                    winds = 0
                    heading = 0
                    m.zero_wind_data()
                if "AERO_FILE" in line:
                    aero_file = line[len("AERO_FILE"):].strip()
                    files.aero_file = aero_file
                    m.read_aero_data(aero_file)
                
                if "PROP_FILE" in line:
                    prop_file = line[len("PROP_FILE"):].strip()
                    files.prop_file = prop_file
                    m.read_propulsion_data(prop_file)
                    
                

    return

def main():
    st = time()
    class files:
        mission_file = "D:/CollegeAssignments/2021Fall/AEE468-Aircraft_Systems_Design/X15 Data/mission.inp"
        aero_file = "D:/CollegeAssignments/2021Fall/AEE468-Aircraft_Systems_Design/X15 Data/x15-aerodata.out"
        prop_file = "D:/CollegeAssignments/2021Fall/AEE468-Aircraft_Systems_Design/X15 Data/ISP235sec-rocket-engine.txt"
    MissionCode(files)
    et = time()
    print(f"File took {round(et-st,5)} seconds to run")
if __name__ == '__main__':
    main()
