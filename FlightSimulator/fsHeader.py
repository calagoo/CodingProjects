import numpy as np
import pandas as pd
import sys
import os
from scipy import interpolate
from time import time

os.system("cls")


## Using a MissionClass because this is how global variables are (generally) set in python. Using "Global" is not the most ideal solution.
## Also its good practice with OOP.

class MissionClass:
    def __init__(self,files) -> None:
    
        # files
        self.aero_file = files.aero_file
        self.prop_file = files.prop_file


        # aero data
        self.machtrim = np.full((50),np.nan)
        self.AoATrim = np.full((50,50),np.nan)
        self.CLtrim = np.full((50,50),np.nan)
        self.CDtrim = np.full((50,50),np.nan)
        self.AltRe = np.full((50),np.nan)
        self.machRe = np.full((50),np.nan)
        self.DeltaCDRe = np.full((50,50),np.nan)
        self.MACHbuffet = np.full((50),np.nan)
        self.CLbuffet = np.full((50),np.nan)
        self.CPoint1 = np.full((50),np.nan)
        self.CPoint2 = np.full((50),np.nan)
        self.CPoint3 = np.full((50),np.nan)

        #propulsiondata
        self.machprop = np.full((50),np.nan)
        self.altprop = np.full((50),np.nan)
        self.plaprop = np.full((50),np.nan)
        self.alfprop = np.full((50),np.nan)
        self.thrustprop = np.full((50,50,50),np.nan)
        self.FFprop = np.full((50,50,50),np.nan)


        #windsdatafor
        self.windAlt = np.full((400),np.nan)
        self.windNSpeed = np.full((400),np.nan)
        self.windESpeed = np.full((400),np.nan)

        #otherpublicvariables
        self.GAMMA = 0
        self.radiusearth = 0
        self.netlift = 0

        self.maxPLA = 0
        self.minPLA = 0
        self.nmachE = 0
        self.naltE = 0
        self.nplaE = 0
        self.nalfE = 0
        self.Sref = 0
        self.Wi = 0
        self.Wf = 0
        self.CLmax = 0
        self.nmach = 0
        self.nAoA = 0
        self.nalt = 0
        self.neng = 0
        self.nmachBuffet = 0
        self.xcolumn = 0
        self.ycolumn = 0
        self.zcolumn = 0
        self.mode2d3d = 0
        self.machmax = 0
        self.machmax1 = 0
        self.altmax = 0
        self.machmin = 0

        self.naltW = 0

        self.cdincrement = 0

        self.winds = 0
        self.iFF = 0
        self.loitertime = 0

        self.Pi = 3.14159265359
        self.DtoR = 3.14159265359/180
        
        ## Init Functions
        self.read_aero_data()
        self.read_propulsion_data()

        ## Check if wind file attached, if so run the "read_wind_data()" function
        if hasattr(files,"wind_file"):
            self.wind_file = files.wind_file
            self.read_wind_data()
        else:
            self.wind_file = False
            self.zero_wind_data()
            print('No Wind File Present')

    ## Reading Data Files
    def readuntil(self,file,input_string):
        '''
        Reads lines from file until it finds the input_string, then breaks the while loop
        '''
        while line:=file.readline():    # Uses walrus operator to stop reading lines at the end of the file
            if input_string in line:
                return
        print("Error reading file:",self.aero_file)
        sys.exit(1)

    def read_aero_data(self):
        try:
            with open(self.aero_file,'r') as f:
                self.readuntil(f,"* PARSED INPUTS")
                self.readuntil(f,"* EDET RESULTS")
                self.readuntil(f,"REFERENCE AREA")
                self.Sref = float(f.readline()[32:].split()[0])

                self.readuntil(f,"NUMCL")
                self.nAoA = int(f.readline())
                self.readuntil(f,"NUMMACH")
                self.nmach = int(f.readline())
                self.readuntil(f,"NALT")
                self.nalt = int(f.readline())

                self.readuntil(f,"CHANGE IN DRAG")
                f.readline()    # Skip a line


                self.machmax = -999
                self.machmin = 999
                self.altmax = -1

                for j in np.arange(0, self.nalt):
                    for i in np.arange(0, self.nmach):
                        # From inside out, this next line reads the next line, splits it by spaces, uses a map object to convert all entities of the list into floats, then converts it back to a list
                        Alt, Mach, deltacd = list(map(lambda x: float(x),f.readline().split()))

                        self.AltRe[j] = Alt
                        self.machRe[i] = Mach
                        self.DeltaCDRe[i][j] = deltacd

                        if Mach < self.machmin:
                            self.machmin = Mach
                        if Mach > self.machmax:
                            self.machmax = Mach
                        if Alt > self.altmax:
                            self.altmax = Alt

                self.readuntil(f,"DRAG POLARS")
                self.readuntil(f,"* MACH")
                self.CLmax = 0
                for i in np.arange(0, self.nmach):
                    for j in np.arange(0, self.nAoA):
                        Mach, AoA, CL, CD = list(map(lambda x: float(x),f.readline().split()))
                        
                        self.machtrim[i] = Mach
                        self.AoATrim[i][j] = AoA
                        self.CLtrim[i][j] = CL
                        self.CDtrim[i][j] = CD

                        if CL > self.CLmax:
                            self.CLmax = CL
                
                self.readuntil(f,"BUFFET")

                self.nmachBuffet = self.nmach
                for i in np.arange(0, self.nmachBuffet):
                    line = f.readline()
                    self.MACHbuffet[i] = float(line[:10])
                    if line[35:] == '':
                        self.CLbuffet[i] = 1
                    else:
                        self.CLbuffet[i] = float(line[35:])

        except Exception as e:
            print("Error reading file:",self.aero_file)
            print(e)
        
    def read_propulsion_data(self):

        try:
            with open(self.prop_file,'r') as f:
                self.readuntil(f,"NPLA")
                self.nplaE = int(f.readline())
                self.readuntil(f,"NMACH")
                self.nmachE = int(f.readline())
                self.readuntil(f,"NALT")
                self.naltE = int(f.readline())
                nAoAE = 1

                self.readuntil(f,"DATA")

                self.machmax1 = -999
                self.maxPLA = -999
                self.minPLA = 999
                self.altmax = -1
                for i in np.arange(0, self.nmachE):
                    for j in np.arange(0, self.naltE):
                        for K in np.arange(0, self.nplaE):
                            # read five column data
                            Mach, Alt, PLA, THRUST, TSFC = list(map(lambda x: float(x),f.readline().split()))
                            self.machprop[i] = Mach
                            self.altprop[j] = Alt
                            self.plaprop[K] = PLA
                            self.thrustprop[i][j][K] = THRUST
                            self.FFprop[i, j, K] = THRUST * TSFC
                            
                            if Mach > self.machmax1:
                                self.machmax1 = Mach
                            if Alt > self.altmax:
                                self.altmax = Alt
                            
                            if PLA > self.maxPLA:
                                self.maxPLA = PLA
                            if PLA < self.minPLA:
                                self.minPLA = PLA
            if self.machmax1 < self.machmax:
                self.machmax = self.machmax1
            
            self.mincruisemach = 0.3
            self.maxcruisemach = self.machmax1 - 0.01
            self.mincruisealt = 10000
            self.maxcruisealt = self.altmax - 1000
            self.minloiteralt = 5000
            self.maxloiteralt = self.altmax - 1000
            self.minloitermach = 0.3
            self.maxloitermach = self.machmax1 - 0.01

        except Exception as e:
            print("Error reading file:",self.prop_file)
            print(e)

    def read_wind_data(self):
        try:
            with open(self.wind_file,'r') as f:
                self.readuntil(f,"N_TOKENS")
                self.naltW = int(f.readline())
                f.readline()    # Skips two lines
                f.readline()

                for i in np.arange(0,self.naltW):
                    # read five column data
                    Alt, WIND, heading, nspeed, espeed = list(map(lambda x: float(x),f.readline().split()))
                    self.windAlt[i] = Alt
                    self.windNSpeed[i] = nspeed
                    self.windESpeed[i] = espeed
        except Exception as e:
            print("Error reading file:",self.wind_file)
            print(e)

    def Atmos76(self,Z):
        '''
        *********** 1976 STANDARD ATMOSPHERE SUBROUTINE **********

            Mason's BASIC program

            W.H. Mason
            Department of Aerospace and Ocean Engineering
            Virginia Tech, Blacksburg, VA 24061
            email: mason@aoe.vt.edu


        extended to higher altitudes by T Takahashi, July 2021

            z  - input altitude, in feet

            output:
                            units: English
            t  - temp.               deg R
            p  - pressure            lb/ft**2
            r  - density (rho)       slug/ft**3
            a  - speed of sound      ft/sec
            mu - viscosity           slug/<ft sec)

            ts - t/t at sea level
            rr - rho/rho at sea level
            pp - p/p at sea level

            rm - Reynolds number per Mach per unit of length
            qm - dynamic pressure/Mach**2


        Altitude in an aircraft is generally measured by the hydrostatic equation:
        p=rho*g*h, where p is the pressure at the point of measurement, rho is the
        density at the point of measurement, g is the acceleration due to gravity at the
        point and h is the height from a reference to that point (the reference is generally
        taken as sea level).

        Aircraft use the hydrostatic equation to determine the height/ altitude because
        pressure can be easily measured with a pitot tube that planes have. So using a
        pitot tube the airplanes measure the pressure and with that they can put it into
        the equation and solve for the height.

        However, gravity is not the same at dself.iFFerent altitudes and changes with respect
        to the altitude. It is very dself.iFFicult for an airplane to measure gravity in the air.
        Therefore airplanes generally measure geopotential altitude. The geopotential
        altitude uses gravity at sea level and takes it to be constant. Whereas geometric
        altitude uses gravity at the point of measurement.
        Therefore p = rho * g0 * h(geopotential)
        where g0 is the gravity at sea-level and h(geopotential) is the geopotential
        altitude
        and
        p = rho * g * h(geometric)
        where g is the gravity at the point of measurement and h(geometric) is the
        geometric altitude or the actual height above sea-level

        at low altitudes, we often conflate geometric (Z) and geopotential (h) altitudes seeing they are so similar

        we may need to rethink "g" ratios
        '''

        KK = 0
        K = 34.163195
        C1 = 0.0003048
        TL = 518.67
        PL = 2116.22
        RL = 0.0023769
        Al = 1116.45
        ML = 0.00000037373
        BT = 0.000000030450963

        h = C1 * Z / (1 + C1 * Z / 6356.766)

        # standard atmosphere through 32km
        if h <= 11:
            t = 288.15 - 6.5 * h
            PP = (288.15 / t) ** (-K / 6.5)
        elif h <= 20:
            t = 216.65
            PP = 0.22336 * np.exp(-K * (h - 11) / 216.65)
        elif h < 32:
            t = 216.65 + (h - 20)
            PP = 0.0540329 * (216.65 / t) ** K

        # new code TTT 7/2021
        elif h < 47:
            t = 228.65 + 2.8 * (h - 32)
            PP = 0.0085666784 * (228.65 / t) ** (K / 2.8)
        elif h < 51:
            t = 270.65
            PP = 0.0010945601 * np.exp(-K * (h - 47) / 270.65)
        elif h < 71:
            t = 270.65 - 2.8 * (h - 51)
            PP = 0.00066063531 * (270.65 / t) ** (-K / 2.8)
        elif h < 84:
            t = 214.65 - 2 * (h - 71)
            PP = 0.000039046834 * (214.65 / t) ** (-K / 2)
        else:
            t = 186.946
            PP = 0.00000368501 * np.exp(-K * (h - 84) / 186.946)

        RR = PP / (t / 288.1)
        MU = BT * t ** 1.5 / (t + 110.4)
        TS = t / 288.15
        a = Al * np.sqrt(TS)
        t = TL * TS
        R = RL * RR
        p = PL * PP
        RM = R * a / MU
        qm = 0.7 * p
        return t, p, R, a, MU, TS, RR, PP, RM, qm
        

    def calcHot1(self,Tmp_hot, e, qw_cold, h_tot):
        sigma = 4.76114E-13  # boltzman constant BTU/(sec ft**2 R**4)
        cp_hot, cv_hot, gamma_hot = self.PerfectGas(Tmp_hot)
        qw_hot_left = qw_cold - cp_hot * Tmp_hot / \
            h_tot * qw_cold  # BTU/(sec ft**2)
        qw_hot_right = e * sigma * Tmp_hot ** 4  # BTU/(sec ft**2)
        return qw_hot_right - qw_hot_left  # BTU/(sec ft**2)
    
    def PerfectGas(self,Tinf): # thermally perfect gas
        AirMW = 28.97 # Air Molecular Weight lbm/lbmol
        Ru = 1.9858 # Universal Gas Constant BTU*lb/mol/R
        R = Ru / AirMW # Gas Const for Air in Btu/lbm-R

        # Values from Chemical and Process Thermodynamics 3/E by Kyle, B. G
        a = 6.713
        B = 0.0002609
        c = 0.000000354
        d = -0.00000000008052

        cp = (a + B * Tinf + c * Tinf ** 2 + d * Tinf ** 3) / AirMW # (BTU/lbm R)
        cv = cp - R # BTU/(lbm R)
        gam = cp / cv # unitless
        return cp, cv, gam

    def calcKCAS(self,Mach, Alt):
        # iterative KCAS solver - revised Nov 21, 2019
        t, Pstatic, R, a, MU, TS, RR, PP, RM, qm = self.Atmos76(Alt)

        if Mach <= 1:
            Qc = Pstatic * ((1 + (1.4 - 1) / 2 * Mach ** 2) ** (1.4 / (1.4 - 1)) - 1)
            KCAS = 661.45 * np.sqrt(2 / (1.4 - 1) *
                                    ((Qc / 2116.22 + 1) ** ((1.4 - 1) / 1.4) - 1))
        else:
            Qc = Pstatic * (((1.4 + 1) / 2 * Mach ** 2) ** (1.4 / (1.4 - 1)) * ((1.4 + 1) / (
                1 - 1.4 + 2 * 1.4 * Mach ** 2)) ** (1 / (1.4 - 1)) - 1)  # herrington eqn 5.31 p. 5.22
            KCAS = 0
            for VC1 in np.arange(2, 661):
                # herrington eqn 5.33, p. 5.23
                Qctarget = 2116.22 * ((1 + 0.2 * (VC1 / 661.4) ** 2) ** 3.5 - 1)
                if Qctarget >= Qc:
                    KCAS = VC1
                    break

            if KCAS == 0:
                for VC2 in np.arange(665, 2200):
                    # Herrington eqn. 5.34, p. 5.23
                    Qctarget = 2116.22 * \
                        (166.921 * (VC2 / 661) ** 7 / (7 * (VC2 / 661) ** 2 - 1) ** 2.5 - 1)
                    if Qctarget >= Qc:
                        KCAS = VC2
                        break
        return KCAS

    def calcKTAS(self,Mach, Alt):
        t, Pstatic, R, a, MU, TS, RR, PP, RM, qm = self.Atmos76(Alt)
        return a * Mach * 3600 / 6076

    def calcRadEqu1(self,q_cold, Mach, a, StaticTmp, e):
        '''
        Detra stagnation point heating
        '''

        cp_inf, cv_inf, gamma_inf = self.PerfectGas(StaticTmp)

        # a is speed of sound ft/s
        # e is the Emissivity
        # StaticTmp in R
        # q_cold in

        unit_con = 25037  # converts ft^2/sec^2 to BTU/lbm
        h_tot = cp_inf * StaticTmp + (Mach * a) ^ (2) / 2 / unit_con  # BTU/lbm

        Tmp_min = 200  # Rankine
        Tmp_max = 3500  # Using a VERY wide range for the search, feel free to tighten up for efficiency

        return self.goldSection(Tmp_min, Tmp_max, e, q_cold, h_tot)
    
    def goldSection(self,xLow, xUpp, e, qw_cold, h_tot):  #Finding the radiative equilibrium temperature

        fLow = self.calcHot1(xLow, e, qw_cold, h_tot)
        fUpp = self.calcHot1(xUpp, e, qw_cold, h_tot)
        if fLow * fUpp > 0:
            print ("Bad initial guess for golden section")
            
        x1g = (1 - 0.38197) * xLow + 0.38197 * xUpp
        f1g = np.abs(self.calcHot1(x1g, e, qw_cold, h_tot))
        x2g = 0.38197 * xLow + (1 - 0.38197) * xUpp
        f2g = np.abs(self.calcHot1(x2g, e, qw_cold, h_tot))
        deltaX = 0.1 # i want the change in T to be less then 0.1 R
            
        Kiter = 3
        while (Kiter < 50 and np.abs(xLow - xUpp) > deltaX):
            Kiter = Kiter + 1
            if f1g > f2g:
                xLow = x1g
                fLow = f1g
                x1g = x2g
                f1g = f2g
                x2g = 0.38197 * xLow + (1 - 0.38197) * xUpp
                f2g = np.abs(self.calcHot1(x2g, e, qw_cold, h_tot))
            else:
                xUpp = x2g
                fUpp = f2g
                x2g = x1g
                f2g = f1g
                x1g = (1 - 0.38197) * xLow + 0.38197 * xUpp
                f1g = np.abs(self.calcHot1(x1g, e, qw_cold, h_tot))
            
        if Kiter > 49:
            print ("Gold Section did not converge")
        else:
            return (x2g + x1g) / 2

    def calcSpherT(self,radius, RR, Mach, a):
        #Detra stagnation point heating

        Rn = radius # ft
        Uco = 26082 # ft/sec   this is the circular orbit velocity
        qo = 17600 / (Rn) ^ 0.5 * (RR) ^ 0.5 * (Mach * a / Uco) ^ 3.15 #BTU/(sec ft^2) # the cold surface heating
        return qo #BTU/(sec ft^2)
        #if you wanted a swept cylinder AKA a leading edge
        #qo = qo * Cos(sweep * DtoR) * Cos(AoA * DtoR) #BTU/(ft^2 sec)# AoA and sweep need to be 0 to be a sphere, when they are non-zero it is a swept cylinder

    def zero_wind_data(self):
        self.naltW = 4
        for i in np.arange(0,self.naltW):
            self.windAlt[i] = (i)*20000
            self.windNSpeed[i] = 0
            self.windESpeed[i] = 0

    ################# Interpolations #################

    def interp2d(self,x,y,z,a,b):
        '''
        A 2D interpolation function that takes in:
        x : 1D array
        y : 1D or 2D array
        z : 2D array

        a : interpolation value 1
        b : interpolation value 2

        Currently the method of interpolation is only Bilinear
        '''
        # Check x,y,z
        if len(list(np.shape(x))) > 1:
            print("ERROR: x must be a 1D array")
            return
        if len(list(np.shape(y))) > 2:
            print("ERROR: y must be a 1D or 2D array")
            return
        if len(list(np.shape(z))) != 2:
            print("ERROR: z must be a 2D array")
            return

        if x[0] >= a:
            i0 = -1
        elif x[-1] <= a:
            i0 = len(x) - 1
        else:
            # This line filters the x list by values less than a.
            # it then finds the max value of the new list and finds where it is in the original array
            i0 = np.where(x==(max(list(filter(lambda i: i <= a, x)))))[0][0]

        # Check y
        if len(list(np.shape(y))) > 1:
            if i0 == -1:
                y = y[0]
            else:
                y = y[i0]
        
        y = y[~np.isnan(y)]

        if y[0] >= b:
            j0 = -1
        elif y[-1] <= b:
            j0 = len(y) - 1
        else:
            j0 = np.where(y==(max(list(filter(lambda i: i <= b, y)))))[0][0]

        new_x = x
        if i0 == -1: 
            new_x = [x[0],x[1]]
        elif i0 == len(x) - 1:
            new_x = [x[-2],x[-1]]
        else:
            new_x = [x[i0],x[i0+1]]

        new_y = y
        if j0 == -1: 
            new_y = [y[0],y[1]]
        elif j0 == len(y) - 1:
            new_y = [y[-2],y[-1]]
        else:
            new_y = [y[j0],y[j0+1]]

        new_z = z[0]
        if i0 == -1:
            ix = 0
        if j0 == -1:
            jx = 0
        if i0 == len(x) - 1:
            ix = len(x) - 2
        if j0 == len(y) - 1:
            jx = len(y) - 2
        if -1 < i0 < len(x) - 1:
            ix = i0
        if -1 < j0 < len(y) - 1:
            jx = j0

        new_z = [[z[ix][jx],z[ix][jx+1]],[z[ix+1][jx],z[ix+1][jx+1]]]

        # Setting up in bilinear interpolation format from: https://x-engineer.org/bilinear-interpolation/
        Q11 = new_z[0][0]
        Q12 = new_z[0][1]
        Q21 = new_z[1][0]
        Q22 = new_z[1][1]

        x1 = new_x[0]
        x2 = new_x[1]
        y1 = new_y[0]
        y2 = new_y[1]

        r1 = lambda x : (Q11 * (x2 - x) / (x2 - x1)) + (Q21 * (x - x1) / (x2 - x1))
        r2 = lambda x : (Q12 * (x2 - x) / (x2 - x1)) + (Q22 * (x - x1) / (x2 - x1))
        p = lambda x,y: (r1(x)*(y2 - y) / (y2 - y1)) + (r2(x)*(y - y1) / (y2 - y1))

        return p(a,b)

    def interpolateCL(self,Mach,AoA):

        # Fix data to put into interp2d function
        xar = self.machtrim[:self.nmach]
        yar = self.AoATrim[:self.nAoA]
        # List comprehension of the z array to take the first n values from the each row
        zar = np.array([i[:self.nAoA] for i in self.CLtrim[:self.nmach]])

        # Remove nans from z array then reshape it so it is a 2d array again
        zar = (zar[~np.isnan(zar)]).reshape(self.nmach,self.nAoA)

        iCL = self.interp2d(xar,yar,zar,Mach,AoA)
        return iCL

    def interpolateCD(self,Mach,CL):
        if Mach > 4:
            return 99
        # Fix data to put into interp2d function
        xar = self.machtrim[:self.nmach]
        yar = self.CLtrim[:self.nAoA]
        # List comprehension of the z array to take the first n values from the each row
        zar = np.array([i[:self.nAoA] for i in self.CDtrim[:self.nmach]])

        # Remove nans from z array then reshape it so it is a 2d array again
        zar = (zar[~np.isnan(zar)]).reshape(self.nmach,self.nAoA)

        iCD = self.interp2d(xar,yar,zar,Mach,CL)
        return iCD

    def interpolateAoA(self,Mach, CL):
        # Fix data to put into interp2d function
        xar = self.machtrim[:self.nmach]
        yar = self.CLtrim[:self.nAoA]
        # List comprehension of the z array to take the first n values from the each row
        zar = np.array([i[:self.nAoA] for i in self.AoATrim[:self.nmach]])

        # Remove nans from z array then reshape it so it is a 2d array again
        zar = (zar[~np.isnan(zar)]).reshape(self.nmach,self.nAoA)

        iAoA = self.interp2d(xar,yar,zar,Mach,CL)
        return iAoA

    def interpolateTSFC(self,Mach, Alt, THRUST):

        if THRUST == 0: 
            return 999
        
        for PLA0 in  np.arange(self.maxPLA, self.minPLA,-0.01):
            PLA0 = int(PLA0 * 100 + 0.5) / 100
            t0 = self.interpolateThrust(Mach, Alt, PLA0)
            if t0 <= THRUST:
                break

        self.PLA = PLA0

        if self.iFF / t0 > 999:
            return 999

    def interpolateThrust(self,Mach, Alt, PLA):
        # This function is directly out of the VBA file. Do not want to make a 3 way interpolator just yet...

        # outofbounds self.flag
        # = 0 : data o.k.
        # = 1 : data marginal, overspeed
        # = 2 : data marginal, overalt
        # = 3 : data corrupted, out-of-bounds

        self.flag = 0


        if self.machprop[0] > Mach or self.machprop[-1] < Mach: 
            self.iFF = 99999
            return 0.001
        if self.altprop[0] > Alt or self.altprop[-1] < Alt: 
            self.iFF = 99999
            return 0.001
        if self.plaprop[0] > PLA or self.plaprop[-1] < PLA: 
            self.iFF = 99999
            return 0.001
        
        i0,j0,k0 = 1,1,1
        
        if self.machprop[-1] <= Mach:
            i0 = self.nmachE - 2             ## CHECK HERE
            self.flag = 1
        else:
            for i in np.arange(0, self.nmachE - 2):
                if self.machprop[i] <= Mach and self.machprop[i + 1] > Mach:
                    i0 = i

        if self.altprop[self.naltE] <= Alt:
            j0 = self.naltE - 2
            self.flag = 2
        else:
            for j in np.arange(0,self.naltE - 2):
                if self.altprop[j] <= Alt and self.altprop[j + 1] > Alt:
                    j0 = j
        
        if PLA == self.plaprop[1]:
            k0 = 2 # GoTo 10
        elif PLA == self.plaprop[self.nplaE]:
            k0 = self.nplaE # GoTo 10
        else:
            k0 = 2
            for K in np.arange(2, self.nplaE):
                if self.plaprop[K - 1] < PLA and self.plaprop[K] >= PLA:
                    k0 = K
        
        # if the lower hinge point in the interp is bad.... then the whole point is bad
        # 10:
        if self.thrustprop[i0, j0, k0] < -990 or self.thrustprop[i0 + 1, j0, k0] < -990 or self.thrustprop[i0, j0 + 1, k0] < -990:
            self.iFF = 999999
            return 0.0001

        deltamach = (Mach - self.machprop[i0])
        if self.thrustprop[i0 + 1, j0, k0] == -999:
            dthrustdMach = 0
            dffdmach = 0
            self.flag = 1
            dtsfcdmach = 0
        deltaalt = (Alt - self.altprop[j0])
        deltapla = (PLA - self.plaprop[k0 - 1])

        if k0 == 0:
            self.iFF = 999999
            return 0.0001

        # Weighted Interpolation Based Upon Mach
        fx0 = deltamach / (self.machprop[i0 + 1] - self.machprop[i0])
        fx1 = deltaalt / (self.altprop[j0 + 1] - self.altprop[j0])
        fx2 = deltapla / (self.plaprop[k0] - self.plaprop[k0 - 1])
                
        k000 = (1 - fx2) * (1 - fx1) * (1 - fx0)
        k001 = fx2 * (1 - fx1) * (1 - fx0)
        k010 = (1 - fx2) * fx1 * (1 - fx0)
        k011 = fx2 * fx1 * (1 - fx0)
        k100 = (1 - fx2) * (1 - fx1) * fx0
        k101 = fx2 * (1 - fx1) * fx0
        k110 = (1 - fx2) * fx1 * fx0
        k111 = fx2 * fx1 * fx0

        total_thrust = self.thrustprop[i0, j0, k0 - 1] * k000 + self.thrustprop[i0, j0, k0] * k001 + \
                        self.thrustprop[i0, j0 + 1, k0 - 1] * k010 + self.thrustprop[i0, j0 + 1, k0] * k011 + \
                        self.thrustprop[i0 + 1, j0, k0 - 1] * k100 + self.thrustprop[i0 + 1, j0, k0] * k101 + \
                        self.thrustprop[i0 + 1, j0 + 1, k0 - 1] * k110 + self.thrustprop[i0 + 1, j0 + 1, k0] * k111
                        
        total_ff = self.FFprop[i0, j0, k0 - 1] * k000 + self.FFprop[i0, j0, k0] * k001 + \
                        self.FFprop[i0, j0 + 1, k0 - 1] * k010 + self.FFprop[i0, j0 + 1, k0] * k011 + \
                        self.FFprop[i0 + 1, j0, k0 - 1] * k100 + self.FFprop[i0 + 1, j0, k0] * k101 + \
                        self.FFprop[i0 + 1, j0 + 1, k0 - 1] * k110 + self.FFprop[i0 + 1, j0 + 1, k0] * k111
        


        self.iFF = total_ff
        
        return total_thrust

    def interpolateReynoldsNumberEffect(self,Mach,Alt):
        # Fix data to put into interp2d function
        xar = self.machRe[:self.nmach]
        yar = self.AltRe[:self.nalt]
        # List comprehension of the z array to take the first n values from the each row
        zar = np.array([i[:self.nalt] for i in self.DeltaCDRe[:self.nmach]])
        
        # Remove nans from z array then reshape it so it is a 2d array again
        zar = (zar[~np.isnan(zar)]).reshape(self.nmach,self.nalt)

        iRe = self.interp2d(xar,yar,zar,Mach,Alt)
        return iRe

    def interpolateCLmax(self,Mach):
        xar = self.MACHbuffet[:self.nmachBuffet]
        yar = self.CLbuffet[:self.nmachBuffet]
        # A 1D interpolation using extrapolation setting
        iCLmax = interpolate.interp1d(xar,yar,fill_value="extrapolate")(Mach)
        return iCLmax

    #### Clean Up Functions here, then paste them in
    def interpolateWind(self,Alt,heading):
        xar = self.windAlt[:self.naltW]
        
        yar = self.windESpeed[:self.naltW]
        iWindE = interpolate.interp1d(xar,yar,fill_value="extrapolate")(Alt)
        
        yar = self.windNSpeed[:self.naltW]
        iWindN = interpolate.interp1d(xar,yar,fill_value="extrapolate")(Alt)
        
        # np.radians(heading) also available
        heading_rad = heading*(np.pi/180)
        iWind = np.cos(heading_rad) * iWindN + np.sin(heading_rad) * iWindE
        return
    ####