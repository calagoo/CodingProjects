import os
import time
import numpy as np
import pandas as pd
import matplotlib.dates as plt
import matplotlib.pyplot as plt
from detrend_lfilter import detrend_lfilter

def android_GNSS_log_reader(filename):
    class messageID:
        class raw:
            list = []
            head = ''
            df = 0
        class accel:
            list = []
            head = ''
            df = 0
        class gyro:
            list = []
            head = ''
            df = 0
        class mag:
            list = []
            head = ''
            df = 0
        class fix:
            list = []
            head = ''
            df = 0
        class nav:
            list = []
            head = ''
            df = 0
        class status:
            list = []
            head = ''
            df = 0
    ids = [messageID.raw,messageID.accel,messageID.gyro,messageID.mag,messageID.fix,messageID.nav,messageID.status]

    with open(filename,"r") as f:
        for line in f:
            # Checks if header (helps runtime)
            if line.startswith('#'):
                if line.startswith("# Version:"):
                    ## This contains info about the platform taking data
                    res = line[2:].split()
                    info = {res[0]:res[1],res[2]:res[3],res[4]:res[5],res[6]:' '.join(res[7:])}
                ## These lines contains the header info
                if line.startswith("# Raw"):messageID.raw.head = line[2:].strip().split(",")[1:]
                if line.startswith("# Accel"):messageID.accel.head = line[2:].strip().split(",")[1:]
                if line.startswith("# Gyro"):messageID.gyro.head = line[2:].strip().split(",")[1:]
                if line.startswith("# Mag"):messageID.mag.head = line[2:].strip().split(",")[1:]
                if line.startswith("# Fix"):messageID.fix.head = line[2:].strip().split(",")[1:]
                if line.startswith("# Nav"):messageID.nav.head = line[2:].strip().split(",")[1:]
                if line.startswith("# Status"):messageID.status.head = line[2:].strip().split(",")[1:]
                ##
            ## These lines collects the raw data
            if line.startswith("Raw"):messageID.raw.list.append(line.strip().split(',')[1:])
            if line.startswith("Accel"):messageID.accel.list.append(line.strip().split(',')[1:])
            if line.startswith("Gyro"):messageID.gyro.list.append(line.strip().split(',')[1:])
            if line.startswith("Mag"):messageID.mag.list.append(line.strip().split(',')[1:])
            if line.startswith("Fix"):messageID.fix.list.append(line.strip().split(',')[1:])
            if line.startswith("Nav"):messageID.nav.list.append(line.strip().split(',')[1:])
            if line.startswith("Status"):messageID.status.list.append(line.strip().split(',')[1:])
            ##
    f.close()
    for id in ids:
        # For all of the classes, replace empty with NaN's, and convert to float (other than strings)
        id.df = pd.DataFrame(id.list,columns=id.head).replace('',np.NaN).astype(float,errors='ignore')
    # Section to convert the ms to a pandas datetime format
    messageID.raw.df["datetime"] = pd.to_datetime(messageID.raw.df["utcTimeMillis"],unit="ms")
    messageID.accel.df["datetime"] = pd.to_datetime(messageID.accel.df["utcTimeMillis"],unit="ms")
    messageID.gyro.df["datetime"] = pd.to_datetime(messageID.gyro.df["utcTimeMillis"],unit="ms")
    messageID.mag.df["datetime"] = pd.to_datetime(messageID.mag.df["utcTimeMillis"],unit="ms")
    messageID.status.df["datetime"] = pd.to_datetime(messageID.status.df["UnixTimeMillis"],unit="ms")
    #
    class sat_const:
        gps_id = messageID.raw.df.loc[(messageID.raw.df["ConstellationType"]==1)]["Svid"].unique()
    return messageID, sat_const

def main():
    filename = "E:/Python/GNSS_KAGGLE/training/2020-05-14-US-MTV-1/Pixel4_GnssLog.txt"
    messageID, sat_const = android_GNSS_log_reader(filename)
    


    for id in sat_const.gps_id:
        res = messageID.raw.df.loc[(messageID.raw.df["Svid"]==id) & (messageID.raw.df["ConstellationType"]==1)]["AccumulatedDeltaRangeMeters"].replace(0,np.NaN)
        try:
            ydet,yz = detrend_lfilter(0,0,res.reset_index(drop=True).interpolate())
        except:
            print(f"Sat G{int(id):02d} could not be parser")
            continue
        xrng = range(len(res))
        xrng_z = range(len(yz))
        xrng_det = range(len(ydet))

        fig,axes = plt.subplots(2)
        fig.suptitle("Detrending Delta Range from {} {}".format("G",id))
        axes[0].plot(xrng[100:-100],res[100:-100],c="b",label="Raw Data")
        axes[0].plot(xrng_z[100:-100],yz[100:-100],c="k",label="Fitted Data")
        axes[1].plot(xrng_det[100:-100],ydet[100:-100],color="b",alpha=.5,linestyle='-',marker="o",markersize=3,label="Detrended Data")
        
        for ax in axes:
            ax.grid()
            ax.legend()
        plt.tight_layout()
    return

if __name__ == "__main__":
    os.system('cls')
    start_time = time.time()
    main()
    print("Ran in {} seconds".format(round(time.time()-start_time,3)))
    plt.show()