'''
Mission Code - for 2021 - rev Oct 13, 2021


Keyword driven mission code for Aircraft Performance Estimation

includes winds and gross-to-net - for obstacle clearance'   minor bug fixes 3/18/2017
idealized calibrated airspeed mod 6/22/2018
winds aloft mode 3/12/2019
fix supersonic KCAS convergence 11/21/2019
aerothermal mods summer 2021 (w/ Jack Griffin)
out-of-bounds checking 8/30/2021
time over 50-miles, time over 100 km, 350F, 1000F
revised error trapping 11/1/2021


T. Takahashi


Rewritten in Python by Robert Rodriguez
'''
# Pub Variables -- Change Later
import numpy as np
import sys
import os
from time import time
os.system('cls')

global t, Pstatic, R, a, MU, TS, RR, PP, RM, qm, cp_hot, cv_hot, gamma_hot, cp_inf, cv_inf, gamma_inf, e, qw_cold, h_tot, maxPLA, minPLA
global machprop,altprop,plaprop,nmachE,naltE,nplaE,thrustprop,FFprop,machtrim,nmach,CLtrim,nAoA,AoATrim


## Interpolation Section






def interpolateAoA(Mach, CL):

    if machtrim(1) >= Mach:
        i0 = 1
    else:
        if machtrim(nmach) <= Mach:
            i0 = nmach - 1
        else:
            for i in np.arange(1, nmach - 1):
                if Mach >= machtrim(i) and Mach < machtrim(i + 1): i0 = i

    x3 = (Mach - machtrim(i0)) / (machtrim(i0 + 1) - machtrim(i0))

    if CL < CLtrim(i0, 1):
        j0 = 1
    else:
        j0 = nAoA - 1
        for j in np.arange(1,nAoA - 1):
            if CL >= CLtrim(i0, j) and CL < CLtrim(i0, j + 1):
                j0 = j

    x1 = (CL - CLtrim(i0, j0)) / (CLtrim(i0, j0 + 1) - CLtrim(i0, j0))

    if CL < CLtrim(i0 + 1, 1):
        j1 = 1
    else:
        j1 = nAoA - 1
        for j in np.arange(1, nAoA - 1):
            if CL >= CLtrim(i0 + 1, j) and CL < CLtrim(i0 + 1, j + 1):
                j1 = j

    x2 = (CL - CLtrim(i0 + 1, j1)) / (CLtrim(i0 + 1, j1 + 1) - CLtrim(i0 + 1, j1))

    i1 = x1 * (AoATrim(i0, j0 + 1) - AoATrim(i0, j0)) + AoATrim(i0, j0)
    i2 = x2 * (AoATrim(i0 + 1, j1 + 1) - AoATrim(i0 + 1, j1)) + AoATrim(i0 + 1, j1)
    
    return i1 + x3 * (i2 - i1)
    


        
## Reading Data Files
def readuntil(file,input_string):
    '''
    Reads lines from file until it finds the input_string, then breaks the while loop
    '''
    while line:=file.readline():    # Uses walrus operator to stop reading lines at the end of the file
        if input_string in line:
            return
    print("Error reading file:",aero_file)
    sys.exit(1)

def read_aero_data(aero_file):
    global AeroData

    class AeroData:

        ## Initialize Lists
        # Note that the _res lists are used in order to organize the main lists in a form similar to a matrix (it is a 2D List/Array)


        # First Loop (CHANGE IN DRAG COEFFICIENT FROM CRUISE ALTITUDE):
        AltRe = np.zeros((50))
        machRe = np.zeros((50))
        DeltaCDRe = np.zeros((50,50))
        #

        # Second Loop (DRAG POLARS):
        machtrim = np.zeros((50))
        AoATrim = np.zeros((50,50))
        CLtrim = np.zeros((50,50))
        CDtrim = np.zeros((50,50))
        #

        # Third Loop (BUFFET DATA):
        MACHbuffet = np.zeros((50))
        CLbuffet = np.zeros((50))
        #

        try:
            with open(aero_file,'r') as f:
                readuntil(f,"* PARSED INPUTS")
                readuntil(f,"* EDET RESULTS")
                readuntil(f,"REFERENCE AREA")
                Sref = float(f.readline()[32:].split()[0])

                readuntil(f,"NUMCL")
                nAoA = int(f.readline())
                readuntil(f,"NUMMACH")
                nmach = int(f.readline())
                readuntil(f,"NALT")
                nalt = int(f.readline())

                readuntil(f,"CHANGE IN DRAG")
                f.readline()    # Skip a line


                machmax = -999
                machmin = 999
                altmax = -1

                for j in np.arange(0, nalt):
                    for i in np.arange(0, nmach):
                        # From inside out, this next line reads the next line, splits it by spaces, uses a map object to convert all entities of the list into floats, then converts it back to a list
                        Alt, Mach, deltacd = list(map(lambda x: float(x),f.readline().split()))

                        AltRe[j] = Alt
                        machRe[i] = Mach
                        DeltaCDRe[i][j] = deltacd

                        if Mach < machmin:
                            machmin = Mach
                        if Mach > machmax:
                            machmax = Mach
                        if Alt > altmax:
                            altmax = Alt

                readuntil(f,"DRAG POLARS")
                readuntil(f,"* MACH")
                CLmax = 0
                for i in np.arange(0, nmach):
                    for j in np.arange(0, nAoA):
                        Mach, AoA, CL, CD = list(map(lambda x: float(x),f.readline().split()))
                        
                        machtrim[i] = Mach
                        AoATrim[i][j] = AoA
                        CLtrim[i][j] = CL
                        CDtrim[i][j] = CD

                        if CL > CLmax:
                            CLmax = CL
                
                readuntil(f,"BUFFET")

                nmachBuffet = nmach
                for i in np.arange(1, nmachBuffet):
                    line = f.readline()
                    MACHbuffet[i] = float(line[:10])
                    if line[35:] == '':
                        CLbuffet[i] = 1
                    else:
                        CLbuffet[i] = float(line[35:])

            # Convert lists to numpy array
            AltRe = np.array(AltRe)
            machRe = np.array(machRe)
            DeltaCDRe = np.array(DeltaCDRe)
            machtrim = np.array(machtrim)
            AoATrim = np.array(AoATrim)
            CLtrim = np.array(CLtrim)
            CDtrim = np.array(CDtrim)
            
            ## Debug Prints
            # print(f"{AltRe}")
            # print(f"{machRe}")
            # print(f"{DeltaCDRe}")
            # print(f"{machtrim}")
            # print(f"{AoATrim}")
            # print(f"{CLtrim}")
            # print(f"{CDtrim}")
            ## Debug Prints

        except Exception as e:
            print("Error reading file:",aero_file)
            print(e)
    
    return AeroData

def read_propulsion_data(prop_file):
    global PropData
    class PropData:

        ## Initialize Lists Using np.zeros to make an array of size (x,y,z)
        machprop = np.zeros((50))
        altprop = np.zeros((50))
        plaprop = np.zeros((50))
        thrustprop = np.zeros((50,50,50))
        FFprop = np.zeros((50,50,50))

        try:
            with open(prop_file,'r') as f:
                readuntil(f,"NPLA")
                nplaE = int(f.readline())
                readuntil(f,"NMACH")
                nmachE = int(f.readline())
                readuntil(f,"NALT")
                naltE = int(f.readline())
                nAoAE = 1

                readuntil(f,"DATA")

                machmax1 = -999
                maxPLA = -999
                minPLA = 999
                altmax = -1
                for i in np.arange(0, nmachE):
                    for j in np.arange(0, naltE):
                        for K in np.arange(0, nplaE):
                            # read five column data
                            Mach, Alt, PLA, THRUST, TSFC = list(map(lambda x: float(x),f.readline().split()))
                            machprop[i] = Mach
                            altprop[j] = Alt
                            plaprop[K] = PLA
                            thrustprop[i][j][K] = THRUST
                            FFprop[i, j, K] = THRUST * TSFC
                            
                            if Mach > machmax1:
                                machmax1 = Mach
                            if Alt > altmax:
                                altmax = Alt
                            
                            if PLA > maxPLA:
                                maxPLA = PLA
                            if PLA < minPLA:
                                minPLA = PLA
            if machmax1 < AeroData.machmax:
                AeroData.machmax = machmax1
            
            mincruisemach = 0.3
            maxcruisemach = machmax1 - 0.01
            mincruisealt = 10000
            maxcruisealt = altmax - 1000
            minloiteralt = 5000
            maxloiteralt = altmax - 1000
            minloitermach = 0.3
            maxloitermach = machmax1 - 0.01

        except Exception as e:
            print("Error reading file:",prop_file)
            print(e)
    return PropData


