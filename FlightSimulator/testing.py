# from fsHeader import *
from fsHeader_class import *
class files:
    aero_file = "D:/CollegeAssignments/2021Fall/AEE468-Aircraft_Systems_Design/X15 Data/x15-aerodata.out"
    prop_file = "D:/CollegeAssignments/2021Fall/AEE468-Aircraft_Systems_Design/X15 Data/ISP235sec-rocket-engine.txt"
    # wind_file = "test"
st = time()

mission = MissionClass(files)

tx,ty = 1,.4

print(mission.interpolateAoA(tx,ty))

et = time()

# MACH     ALFA      CL      CD
# 0.1		1	0.05965	0.06399
# 0.1		2.5	0.1489	0.06945
print(f"{round(et-st,5)} seconds to run")