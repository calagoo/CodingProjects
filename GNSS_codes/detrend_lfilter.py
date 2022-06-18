from numpy import nan
import numpy as np
from scipy import signal
import time
'''
Function for detrending using a Buttersworth filter method. 
'''
def detrend_lfilter(x,y,y_int):
    b, a = signal.butter(3, 0.05)
    zi = signal.lfilter_zi(b, a)
    z, _ = signal.lfilter(b, a, y_int, zi=zi*y_int[0])
    z2, _ = signal.lfilter(b, a, z, zi=zi*z[0])
    yz = signal.filtfilt(b, a, y_int)
    ydet = yz-y_int
    return ydet,yz

if __name__ == '__main__':
    # Testing
    y_int = [
1,2,3,4,5,6,7,8,9,10,11,12,13,14
]
    y=0
    x=0
    detrend_lfilter(x,y,y_int)