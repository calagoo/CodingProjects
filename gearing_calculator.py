gears = 4
gear_ratio_L1 = 3/1 # 3:1
gear_ratio_L2 = 3/1 # 3:1
gear_ratio_H1 = 2/1 # 5:2
gear_ratio_H2 = 9/5 # 2:1

gear_1 = gear_ratio_L1*gear_ratio_L2
gear_2 = gear_ratio_H1*gear_ratio_L2
gear_3 = gear_ratio_L1*gear_ratio_H2
gear_4 = gear_ratio_H1*gear_ratio_H2

print(f"{gear_1 =   }")
print(f"{gear_2 =   }")
print(f"{gear_3 =   }")
print(f"{gear_4 =   }")