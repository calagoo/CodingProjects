import os
from fsHeader import *
from time import sleep, time
import matplotlib.pyplot as plt

# Todo: 
# Remove match case so < python3.10 is compatible


def MissionCode(files):

    # Turns on debug prints
    debug = True

    main_simdf = pd.DataFrame(columns=["Time (sec)", "Dist (nM)", "Wt (lbm)", "MACH", "alt (ft)", "Q (psf)", "KTAS", "KCASideal (conflate with KIAS)", "AOA", "CL", "CD", "L/D", "M(L/D)", "Drag (lbf)", "Thrust (lbf)", "Lin.Acc (ft/s**2)", "ROC (ft/min)", "percent induced",
                              "FF (lb/hr)", "gamma (deg)", "SR (nM/lbm)", "KEAS (nm/hr)", "HEADING (deg)", "WINDS (nm/hr)", "GROUND SPEED (nm/hr)", "ESAD (nM)", "MODE", "PLA", "TSFC", "STATIC TEMP (R)", "STAG TEMP (R)", "RAD EQUI TEMP (R)", "REF. HEAT. RATE (BTU/FT**2/SEC)", "UNIT RE # (/FT)", "Nz"])
    m = MissionClass(files)

    # lists
    class FlightData:
        time = []
        alt = []
        weight = []
        mach = []
        dist = []
        t_over_50mi = 0
        t_over_100km = 0
        t_over_350F = 0
        t_over_1000F = 0
        peak_stag_T = 0
        peak_eq_wall_T = 0


    # variable inits
    t = 0
    Wi = 0

    ix = 2  # Possibly fix this, setting to 2 only to comply with VBA script
    iy = 50
    iy1 = 25

    minPLA = 0
    maxPLA = 0

    t = 0
    wt = Wi
    mach = 0
    dmach = 0
    alt = 0
    AoA = 0
    dh = 0

    roc = 0
    dist = 0
    creditdist = 0

    bank = 0
    bank_rad = 0

    neng = 0
    fuelburned = 0
    winds = 0
    heading = 0
    creditfuelburned = 0

    m.PLA = maxPLA
    THRUST = 0
    TSFC = 999
    drag = 0
    CL = 0
    CD = 0.0001
    cdincrement = 0
    Nz = 1
    flag = 0
    Flag1 = 0

    creditflag = 1
    linaccel = 0
    machold = mach

    qmax = 9999999
    qmin = 0

    aero_flag = False
    prop_flag = False
    first_flag = True

    fpa = 0
    ESAD = 0

    ke0 = 1
    constq = -1
    constKEAS = -1

    altitude_tolerance = 100
    mach_tolerance = 0.01

    gross2net = 0
    vfreefall = 0
    dalt = 0
    velV = 0
    Fvert = 0
    vfts = 0
    accelV = 0
    APU_burn = 0

    # Stop Conditions
    stop_dist = 9999
    stop_time = 9999
    stop_alt = 99999
    stop_mach = 9999
    stop_wt = 0
    stop_ROC = 999999
    stop_kias = 999999
    stop_gamma = 999

    GAMMA = 0

    loitertime = 0
    alt_flag = 1
    dt = 1  # Set Default Time Step

    tover50mi = 0
    tover100km = 0
    tover350F = 0
    tover1000F = 0
    peakStagT = 0
    peakEqWT = 0

    flownMaxWt = 0
    flownMaxLift = 0
    flownNzmax = 0

    flag = 0
    flag2 = 0
    flag3 = 0
    creditflag = 0
    starttime = 0
    Wf = 0
    itr = True  # Iteration boolean
    exitcode = False
    linenum = 0

    # Begin Primary Loop
    with open(files.mission_file, 'r') as f:
        while wt >= Wf and not exitcode:
            itr = True
            while itr and (line := f.readline()):
                if debug:
                    print(line)
                linenum += 1
                # if itr is False:
                #     break

                # Same format as VBA script
                # Space before and after flight commands
                # This is so similar commands arent mixed ie CONST_KIAS and CLIMB_CONST_KIAS
                line = ' ' + line.upper().strip() + ' '
                # Search for comments (leading *)
                if line.startswith("*"):
                    continue
                if " END " in line:
                    itr = False
                    exitcode = True
                    print("NORMAL TERMINATION")
                    break

                if " STOP " in line:  # Not sure if continue is the correct line here
                    itr = False
                
                if " SET " in line:
                    flag = 0
                
                if " LOG " in line:
                    print(t)
                    print(wt)
                    print(alt)
                    print(mach)
                    print(line)
                    continue

                # Begin of main instruction set
                if " CLIMB " in line:
                    flag = 1
                    m.PLA = maxPLA
                    creditflag = 1
                    starttime = t
                    startdist = dist
                if " CLIMB_CONST_MACH " in line:
                    flag = 1
                    creditflag = 1
                    starttime = t
                    startdist = dist
                if " CLIMB_CONST_KIAS " in line:
                    flag = 2
                    creditflag = 1
                    starttime = t
                    startdist = dist
                if " DESCEND " in line:
                    flag = 1
                    m.PLA = minPLA
                    creditflag = 1
                    starttime = t
                    startdist = dist
                if " DESCEND_CONST_MACH " in line:
                    flag = 1
                    m.PLA = minPLA
                    creditflag = 1
                    starttime = t
                    startdist = dist
                if " DESCEND_CONST_KIAS " in line:
                    flag = 2
                    m.PLA = minPLA
                    creditflag = 1
                    starttime = t
                    startdist = dist
                if " ACCEL " in line:
                    flag = 3
                    m.PLA = maxPLA
                    creditflag = 1
                    starttime = t
                    startdist = dist
                if " DECEL " in line:
                    flag = 3
                    m.PLA = minPLA
                    creditflag = 1
                    starttime = t
                    startdist = dist
                if " LEVEL " in line:
                    flag = 4
                    leftright = 0
                    creditflag = 1
                    starttime = t
                    startdist = dist
                if " CONST_CL " in line:
                    flag = 5
                    leftright = 0
                    creditflag = 1
                    starttime = t
                    startdist = dist
                    CLtarget = 0
                if " CONST_AOA " in line:
                    flag = 6
                    leftright = 0
                    creditflag = 1
                    starttime = t
                    startdist = dist
                    AoAtarget = 0
                if " GROUND_RUNUP " in line:
                    flag = 7
                    m.PLA = maxPLA
                    creditflag = 1
                    starttime = t
                    startdist = dist

                # flag 0 is the initialized variables (although they can be changed throughout the program)
                if flag == 0:
                    if " NENG " in line:
                        neng = float(line[len(" NENG "):])
                        print(neng)
                        if neng < 0:
                            neng = 1
                    if " W_START " in line:
                        Wi = float(line[len(" W_START "):])
                        wt = Wi
                    if " W_END " in line:
                        Wf = float(line[len(" W_END "):])
                    if " WEIGHT " in line:
                        wt = float(line[len(" WEIGHT "):])
                    if " MACH " in line:
                        mach = float(line[len(" MACH "):])
                        machold = mach
                    if " APU_BURN " in line:
                        APU_burn = float(line[len(" APU_BURN "):])
                    if " ALTITUDE " in line:
                        alt = float(line[len(" ALTITUDE "):])
                    if " DROP " in line:
                        wt -= float(line[len(" DROP "):])
                    if " PLA " in line:
                        m.PLA = float(line[len(" PLA "):])
                    if " LOAD_FACTOR " in line:
                        Nz = float(line[len(" LOAD_FACTOR "):])
                    if " WINDS " in line:
                        winds = float(line[len(" WINDS "):]) * 1.689
                    if " w " in line:
                        gross2net = float(line[len(" w "):]) / 100
                    if " DT " in line:
                        dt = float(line[len(" DT "):])
                        if dt < 0.001:
                            dt = 0.001
                        elif dt > 60:
                            dt = 60
                    if " DELTA_CD " in line:
                        cdincrement = float(line[len(" DELTA_CD "):])
                        optflag = 0
                        if debug:
                            print("CD_INCREMENT = {}".format(cdincrement))

                    # Input of wind file
                    if " WIND_FILE " in line:
                        wind_file = line[len(" WIND_FILE "):].strip()
                        files.wind_file = wind_file
                        m.read_wind_data(wind_file)
                        if debug:
                            print("Wind File Changed")

                    if " ZERO_WIND_DATA " in line:
                        wind_file = ""
                        winds = 0
                        heading = 0
                        m.zero_wind_data()
                    if " HEADING " in line:
                        heading = float(line[len(" HEADING "):])

                    if " AERO_FILE " in line:
                        aero_file = line[len(" AERO_FILE "):].strip()
                        files.aero_file = aero_file
                        m.read_aero_data(aero_file)
                        if debug:
                            print("Aero File Changed")

                    if " PROP_FILE " in line:
                        prop_file = line[len(" PROP_FILE "):].strip()
                        files.prop_file = prop_file
                        m.read_propulsion_data(prop_file)
                        if debug:
                            print("Prop File Changed")
                # END FLAG 0

                else:
                    if " DT " in line:
                        dt = float(line[len(" DT "):])
                        if dt < 0.001:
                            dt = 0.001
                        elif dt > 10:
                            dt = 10
                    if " PLA " in line:
                        m.PLA = float(line[len(" PLA "):])
                    if " CL_TARGET " in line:
                        CLtarget = float(line[len(" CL_TARGET "):])
                    if " AOA_TARGET " in line:
                        AoAtarget = float(line[len(" AOA_TARGET "):])
                    if " BANK " in line:
                        bank = float(line[len(" BANK "):])
                        bank_rad = bank * np.pi / 180
                    if " NO DISTANCE CREDIT " in line:
                        creditflag = 0

                    # "STOP" Checking and Parsing
                    if " STOP " in line:
                        stop_line = line.strip().split(" ")
                        # Checking Format
                        if len(stop_line) != 3:
                            print("Line Error: ln {} ({})".format(
                                linenum, line.strip()))
                        else:
                            if stop_line[0] != "STOP":
                                print("Line Error: ln {} ({})".format(
                                    linenum, line.strip()))
                            if not ("<" in stop_line[1] or ">" in stop_line[1]):
                                print("Line Error: ln {} ({})".format(
                                    linenum, line.strip()))
                            if not stop_line[2].isnumeric():
                                print("Line Error: ln {} ({})".format(
                                    linenum, line.strip()))
                        # Done Checking Format
                        arg = stop_line[1]
                        val = float(stop_line[2])

                        # Reseting Stop Values
                        stop_time,stop_alt,stop_mach,stop_dist,stop_wt,stop_ROC = 99999,999999,99999,99999,0,999999
                        stop_kias,stop_batterycharge,stop_gamma,stop_heading_max,stop_heading_min = 999999,0,999,9999,-9999

                        # Parsing
                        match arg:
                            case "TIME>":
                                stop_time = val
                            case "ALT>":
                                stop_alt = val
                            case "ALT<":
                                stop_alt = -val
                            case "GAMMA>":
                                stop_gamma = val
                            case "DIST>":
                                stop_dist = val
                            case "MACH>":
                                stop_mach = val
                            case "MACH<":
                                stop_mach = -val
                            case "WT<":
                                stop_wt = val
                            case "RELATIVE_TIME>":
                                stop_time = starttime + val
                            case "RELATIVE_DIST>":
                                stop_dist = val + startdist
                            case "ROC>":
                                stop_ROC = val
                            case "ROC<":
                                stop_ROC = -val
                            case "KIAS>":
                                stop_kias = val
                            case "KIAS<":
                                stop_kias = -val
                        #
                    # End of "STOP" Checking and Parsing
                    if " CONST_KIAS " in line:
                        constKIAS = float(line[len(" CONST_KIAS "):])
            if flag >= 1:
                # I leave my steps on the top of my loops. The -1 allows it to loop the same.
                ix -= 1
                itr = True
                while itr and not exitcode:

                    ## Adding values to lists
                    FlightData.time.append(t)
                    FlightData.alt.append(alt)
                    FlightData.weight.append(wt)
                    FlightData.mach.append(mach)
                    FlightData.dist.append(dist)
                    ##
                    if THRUST < 5:
                        pass
                    if ix % 50 and debug:  # For Debugging
                        # try:
                        #     print(f"{THRUST    =    }")
                        # except:
                        #     pass
                        # sleep(.05)
                        pass
                    ix += 1  # Timestep essentially. Used in VBA to set the cell row

                    # Used in place of a GoTo function
                    bypass = False

                    # call Bill Mason's 1976 standard atmosphere
                    StaticTmp, p, R, a, MU, TS, RR, PP, RM, qm = m.Atmos76(alt)

                    # calculate q
                    q = qm * mach ** 2
                    if q <= 0:
                        q = 0.001       # THIS MAY BE A PROBLEM IN FUTURE - JACK
                        if debug:
                            print("Warning: Dynamic Pressure = {}".format(q))

                    # velocity in kts and mach
                    V = (a * mach) * 3600 / 6080  # knots
                    vfts = a * mach
                    velV = vfts * np.sin(GAMMA)
                    velH = vfts * np.cos(GAMMA)

                    KEAS = 661 * np.sqrt(q / 1481)
                    KCAS = m.calcKCAS(mach, alt)

                    # Radius of the earth, cylindrical earth model
                    # Todo update to wgs ellipsoid model
                    radius_earth = 6371.009 * 1000 * 39.4 / 12

                    # net lift required is that of the weight of the vehicle minus the centrifugal acceleration due to flight speed
                    # add to include load factor
                    netlift = wt * \
                        ((1 - vfts ** 2 / (alt + radius_earth) / 32.174)) + \
                        wt * (Nz - 1)
                    # AoA not declared yet !!!
                    CL = (netlift * np.cos(GAMMA) - THRUST *
                          np.sin(np.radians(AoA))) / q / m.Sref
                    if CL < 0:
                        CL = 0

                    CLmaxx = m.interpolateCLmax(mach)
                    if CL > CLmaxx:
                        CL = CLmaxx

                    if flag == 5:
                        CL = CLtarget
                    if flag == 6:
                        CL = m.interpolateCL(mach, AoAtarget)
                    if flag == 7:
                        CL = 0
                    CD = m.interpolateCD(
                        mach, CL) + m.interpolateReynoldsNumberEffect(mach, alt) + cdincrement
                    CD0 = m.interpolateCD(
                        mach, 0) + m.interpolateReynoldsNumberEffect(mach, alt) + cdincrement
                    pctind = (CD - CD0) / CD

                    vfts = (a * mach)
                    V = vfts * 3600 / 6080
                    ktas = V  # in ktas
                    KEAS = np.sqrt(q / 1481) * 1116 * 3600 / 6080  # in KEAS
                    drag = CD * m.Sref * q
                    AoA = m.interpolateAoA(mach, CL)
                    if flag == 6:
                        AoA = AoAtarget

                    flag2 = 0  # = 0 level; =1 climb; = 2 accel/decel; =3 combined accel-climb

                    match flag:
                        case 1 | 2:
                            Mode = "CLIMB_CONST_MACH" if flag == 1 else "CLIMB_CONST_KEAS"
                            dt = 1
                            flag2 = 1
                        case 3:                                             # ACCEL/DECEL LOGIC
                            Mode = "ACCEL/DECEL"
                            dt = 1
                            flag2 = 2
                        case 4:                                             # LEVEL FLIGHT
                            Mode = "LEVEL"
                            dt = 1  # default for typical aircraft is 1 minute, now 1 -sec
                            flag2 = 0
                        case 5:                                             # CONSTANT CL "GLIDE" - bank angle is enabled!
                            Mode = "CONST_CL (BANK=" + str(bank) + ")"
                            dt = 1  # default for typical aircraft is 1 minute, now 1 -sec
                            flag2 = 3
                        case 6:                                             # CONSTANT AoA "GLIDE" - bank angle is enabled!
                            Mode = "CONST_AOA (BANK=" + str(bank) + ")"
                            dt = 1  # default for typical aircraft is 1 minute, now 1 -sec
                            flag2 = 3
                        case 7:
                            Mode = "GROUND_RUNUP"
                            dt = 1
                            flag2 = 7

                    match flag2:
                        case 0:  # = 0 level
                            THRUST = drag

                            dmach = 0
                            dalt = 0

                            T1 = m.interpolateThrust(mach, alt, 1)*neng
                            if THRUST > T1:
                                dt = 1
                                THRUST = T1
                                # put all energy into climb/descend
                                dalt = 1 * \
                                    (a * mach) * ((THRUST * np.cos(np.radians(AoA)
                                                                   ) - drag - netlift * np.sin(GAMMA)) / wt) * dt

                            TSFC = m.interpolateTSFC(mach, alt, THRUST/neng)

                        case 1:  # = 1 climb
                            # set fixed power per PLA variable
                            THRUST = m.interpolateThrust(mach, alt, m.PLA)*neng
                            TSFC = m.interpolateTSFC(mach, alt, THRUST / neng)
                            dmach = 0                                               # hold constant speed

                            k0 = 1
                            if alt < 36089 and flag == 1:
                                k0 = 1 / (1 - 0.133184 * mach ** 2)
                            if alt < 36089 and flag == 2:
                                k0 = 1 / (1 + 0.566816 * mach ** 2)
                            if alt > 36089 and flag == 2:
                                k0 = 1 / (1 + 0.7 * mach ** 2)

                            if flag == 1:
                                dmach = 0
                            else:
                                StaticTmp, p, R, a, MU, TS, RR, PP, RM, qm = m.Atmos76(
                                    alt + dh)
                                dmach = np.sqrt(q / qm) - mach

                            SEThrust = ((THRUST - drag) / wt)
                            if SEThrust > 0:
                                # gross to net conversion - 2/9/2017
                                SEThrust -= gross2net
                            # put all energy into climb/descend
                            dalt = k0 * (a * mach) * SEThrust * dt

                        case 2:  # = 2 accel/decel
                            # set fixed power per PLA variable
                            THRUST = m.interpolateThrust(
                                mach, alt, m.PLA) * neng
                            TSFC = m.interpolateTSFC(mach, alt, THRUST / neng)
                            dmach = (THRUST - drag) / wt * 32.2 * dt / \
                                a  # put all energy into accel/decel
                            dalt = 0

                        case 3:  # = 3 combined accel-climb
                            # set fixed power per PLA variable
                            THRUST = m.interpolateThrust(
                                mach, alt, m.PLA) * neng
                            TSFC = m.interpolateTSFC(mach, alt, THRUST / neng)
                            lift = CL * q * m.Sref * np.cos(bank_rad)
                            drag = CD * q * m.Sref
                            Fvert = lift * \
                                np.cos(GAMMA) + THRUST * np.sin(GAMMA +
                                                                np.radians(AoA)) - drag * np.sin(GAMMA) - wt
                            Fhorz = THRUST * \
                                np.cos(GAMMA + np.radians(AoA)) - drag * \
                                np.cos(GAMMA) - lift * np.sin(GAMMA)
                            accelV = Fvert / wt  # in gee's
                            accelH = Fhorz / wt  # in gee's
                            velV = velV + 32.2 * accelV * dt  # in ft/sec
                            velH = velH + 32.2 * accelH * dt  # in ft/sec
                            vfts = np.sqrt(
                                velV * velV + velH * velH)  # in ft/sec
                            dmach = vfts / a - mach  # change in mach
                            dh = velV * dt  # in ft/sec
                            bypass = True

                        case 7:
                            # set fixed power per PLA variable
                            THRUST = m.interpolateThrust(
                                mach, alt, m.PLA) * neng
                            TSFC = m.interpolateTSFC(mach, alt, THRUST / neng)
                            dalt = 0                                                # hold alt
                            dmach = (THRUST - drag) / wt * 32.2 * \
                                dt / a            # add speed
                    if not bypass:
                        # Recalc lift and climbing values
                        if CL == CLmaxx:
                            lift = CL * q * m.Sref
                            vfreefall = vfreefall + \
                                (lift - wt) / wt * 32.2 * dt
                            dh = dalt + vfreefall * dt
                        else:
                            dh = dalt
                            vfreefall = 0
                        if flag2 == 7:
                            dh = 0

                    roc = (dh / dt) * 60

                    # compute delta weight for time step + add fuel burn for APU
                    dw = TSFC * THRUST * dt / 3600 + APU_burn * dt / 3600

                    if dw < 0.001:
                        dw = 0.000000001

                    mach = mach + dmach

                    if TSFC > 999:
                        print("FATAL ERROR: Engine Table Interpolation Error")

                    # increment state variables
                    alt = alt + dh
                    t = t + dt
                    wt = wt - dw

                    # winds aloft file
                    if m.wind_file != "":
                        # winds externally are in nm/hr, internally are in ft/sec
                        winds = m.interpolateWind(alt, heading) * 1.689

                    # in nM <== revised 3/18/2017
                    dist = dist + (a * mach * dt *
                                   np.cos(GAMMA) - winds * dt) / 6080
                    ESAD = ESAD + (a * mach * dt * np.cos(GAMMA)
                                   ) / 6080   # in nM
                    fuelburned = fuelburned + dw
                    FF = dw / dt * 3600

                    # calculate new gamma in radians - limit to 1.4 radians (~80 degrees)
                    GAMMA = np.arcsin(dh / (a * mach * dt))  # rev 3/18/2017
                    if GAMMA > 1.4:
                        GAMMA = 1.4
                    if GAMMA < -1.4:
                        GAMMA = -1.4
                    if creditflag == 1:
                        creditdist = creditdist + \
                            (a * mach * dt - winds * dt) / 6080  # in nM
                    if creditflag == 1:
                        creditfuelburned = creditfuelburned + dw  # in lb

                    linaccel = ((a * mach) - (a * machold)) / dt

                    if linaccel > 400:
                        # 10/23/2021 - error trap if haywire - T
                        print("Linear Acceleration Over Limit: Stopping...")
                        sys.exit()
                    machold = mach

                    # Stop conditions
                    if dist >= stop_dist:
                        itr = False
                        print("DISTANCE STOP")
                    if t >= stop_time:
                        itr = False
                        print("TIME STOP")
                    if stop_alt >= 0 and alt >= stop_alt:
                        itr = False
                        print("ALTITUDE STOP")
                    if stop_alt < 0 and alt <= -stop_alt:
                        itr = False
                        print("ALTITUDE BELOW ZERO:")
                        if debug:
                            print(f"{alt} < {-stop_alt}")
                    if stop_mach > 0 and mach >= stop_mach:
                        itr = False
                        print("MACH STOP")
                    if stop_mach < 0 and mach <= -stop_mach:
                        itr = False
                        print("MACH BELOW ZERO")
                    if (GAMMA / (np.pi/180)) > stop_gamma:
                        itr = False
                        print("GAMMA STOP")
                    if wt <= stop_wt:
                        itr = False
                        print("WEIGHT STOP")
                    if stop_ROC > 0 and roc >= stop_ROC:
                        itr = False
                        print("CLIMB RATE STOP")
                    if stop_ROC < 0 and roc <= -stop_ROC:
                        itr = False
                        print("CLIMB RATE STOP")
                    if stop_kias > 0 and KCAS >= stop_kias:
                        itr = False
                        print("KEAS STOP")
                    if stop_kias < 0 and KCAS <= -stop_kias:
                        itr = False
                        print("KEAS STOP")

                    # check if nonsense altitudes or speeds (revised 8/30/2021)
                    if alt < 0:
                        alt = 0
                        itr = False
                        exitcode = True
                        print("ALTITUDE BELOW SEA-LEVEL!")
                    if alt > 400000:
                        itr = False
                        exitcode = True
                        print("ALTITUDE ABOVE 400K FT!")
                    if mach > 4:
                        itr = False
                        exitcode = True
                        print("MACH ABOVE 4!")
                    if mach < 0:
                        itr = False
                        exitcode = True
                        print("MACH BELOW ZERO")

                    # check if weight goes below zero fuel weight
                    if wt < Wf:
                        itr = False
                        print("WT BELOW MIN")

                    ########################################################################
                    # work aerothermal cases for sphere nose                               #
                    # do basic aerothermal stuff here                                      #
                    # Thermodynamics and Heat Transfer - Takahashi & Rodi & Griffin - 2021 #
                    ########################################################################

                    # Thermally Perfect Gas Model
                    cp, cv, gam = m.PerfectGas(StaticTmp)
                    StagTmp = StaticTmp * (1 + (gam - 1) / 2 * mach ** 2)
                    RefHeatRate = m.calcSpherT(
                        1, RR, mach, a)  # 1ft referencesphere
                    emissivity = 1  # CHANGE ME .. No!
                    EquiWall = m.calcRadEqu1(
                        RefHeatRate, mach, a, StaticTmp, emissivity)
                    Re_ft = RM * mach

                    if alt > 264000:
                        tover50mi += dt
                    if alt > 328300:
                        tover100km += dt
                    if EquiWall > (350 + 459):
                        tover350F += dt
                    if EquiWall > (1000 + 459):
                        tover1000F += dt
                    if peakStagT <= StagTmp:
                        peakStagT = StagTmp
                    if peakEqWT <= EquiWall:
                        peakEqWT = EquiWall



                    # In VBA we begin setting the cells to values here

                    if flownMaxWt <= wt:
                        flownMaxWt = wt
                    if flownMaxLift <= lift:
                        flownMaxLift = lift
                    if flownNzmax <= lift / wt:
                        flownNzmax = lift / wt


                    FlightData.t_over_50mi = tover50mi
                    FlightData.t_over_100km = tover100km
                    FlightData.t_over_350F = tover350F
                    FlightData.t_over_1000F = tover1000F
                    FlightData.peak_stag_T = peakStagT
                    FlightData.peak_eq_wall_T = peakEqWT


                


    return FlightData


def main():
    os.system('cls')
    st = time()

    class files:
        mission_file = "D:/CollegeAssignments/2021Fall/AEE468-Aircraft_Systems_Design/X15 Data/mission.inp"
        aero_file = "D:/CollegeAssignments/2021Fall/AEE468-Aircraft_Systems_Design/X15 Data/x15-aerodata.out"
        prop_file = "D:/CollegeAssignments/2021Fall/AEE468-Aircraft_Systems_Design/X15 Data/ISP235sec-rocket-engine.txt"
    MissionData = MissionCode(files)
    et = time()
    print('\n-----------RESULTS-----------')
    print(f"File took {round(et-st,5)} seconds to run")

    print()
    print("MISSION DATA START:")
    print(f"t_over_50mi         =   {MissionData.t_over_50mi}")
    print(f"t_over_100km        =   {MissionData.t_over_100km}")
    print(f"t_over_350F         =   {MissionData.t_over_350F}")
    print(f"t_over_1000F        =   {MissionData.t_over_1000F}")
    print(f"peak_stag_T         =   {MissionData.peak_stag_T}")
    print(f"peak_eq_wall_T      =   {MissionData.peak_eq_wall_T}")


    plt.figure()
    plt.plot(MissionData.alt)
    plt.title('Altitude')
    plt.grid()

    plt.figure()
    plt.plot(MissionData.time)
    plt.title('time')
    plt.grid()
    
    plt.figure()
    plt.plot(MissionData.mach)
    plt.title('mach')
    plt.grid()
    
    plt.figure()
    plt.plot(MissionData.dist)
    plt.title('dist')
    plt.grid()
    
    plt.figure()
    plt.plot(MissionData.weight)
    plt.title('weight')
    plt.grid()

    plt.show()



if __name__ == '__main__':
    main()
