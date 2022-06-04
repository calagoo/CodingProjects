import numpy as np
import matplotlib.pyplot as plt
import math
import os

# true/falses
loopdebug = False


os.system('cls')

mass_planet = 10
mass_sat = 1

satx = 0
satvx = 2
satax = 0
saty = 3
satvy = 1
satay = 0

satx_arr = []
saty_arr = []

planetx = 0
planety = 0

time = 0
dtime = .1
maxtime = 50
loopnum = 0

while True:
    loopnum += 1
    time += dtime
    distancex = (planetx - satx)
    distancey = (planety - saty)
    distance = (distancex**2 + distancey**2)**.5


    force = (mass_planet*mass_sat)/(distance**2)
    theta = np.arctan2(saty-planety,satx-planetx)
    beta = (math.pi)/2-theta

    satax = -math.sin(beta)*force/mass_sat
    satvx = satvx + satax * dtime
    satx = satx + satvx * dtime
    satay = -math.cos(beta)*force/mass_sat
    satvy = satvy + satay * dtime
    saty = saty + satvy * dtime

    satx_arr.append(satx)
    saty_arr.append(saty)

    if loopdebug:
        print('Loop',loopnum)
        print('__xaxis__')
        print('distx=',round(distancex,4))
        print('accelx=',round(satax,4))
        print('velx=',round(satvx,4))
        print('satX=',round(satx,4))
        print('__yaxis__')
        print('disty=',round(distancey,4))
        print('accely=',round(satay,4))
        print('vely=',round(satvy,4))
        print('saty=',round(saty,4))
        print('\n\n')

    plt.cla()
    plt.xlim([-10,10])
    plt.ylim([-10,10])
    plt.text(-7,7,'X_Velocity: ' + str(round(satvx,5)) + '\nY_Velocity: ' + str(round(satvy,5)))
    plt.text(-7,9,'Time: ' + str(round(time,5)))
    plt.arrow(satx,saty,satax,satay,color='r')
    plt.arrow(satx,saty,satax,0,color='b')
    plt.arrow(satx,saty,0,satay,color='g')

    plt.plot(satx,saty,'bo')
    plt.plot(planetx,planety,'ro')
    plt.plot(satx_arr,saty_arr,'y--')
    plt.draw()
    plt.pause(.03)

    if time >= maxtime:
        break

