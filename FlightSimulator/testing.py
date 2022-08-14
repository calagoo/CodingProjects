# from fsHeader import *
from fsHeader_class import *
class files:
    aero_file = "D:/CollegeAssignments/2021Fall/AEE468-Aircraft_Systems_Design/X15 Data/x15-aerodata.out"
    prop_file = "D:/CollegeAssignments/2021Fall/AEE468-Aircraft_Systems_Design/X15 Data/ISP235sec-rocket-engine.txt"
    # wind_file = "test"
st = time()

mission = MissionClass(files)
print(mission.interpolateCL(.5999,5))

et = time()



print(f"{round(et-st,5)} seconds to run")