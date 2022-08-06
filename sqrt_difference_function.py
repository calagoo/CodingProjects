import time
import numpy as np
import math
import os

os.system('cls')
def sqrt(x):
    if x == 0:
        return 0
    diff = 1
    eps = 1e-10
    b = x/2
    a=0
    while (diff > eps):
        a = (x/b + b)/2
        # diff = abs((b*b)/x-1)
        diff = abs(b-a)
        b = a
    # return int(b)   # Returns truncated value
    return b        # Returns the actual value

x=1028780
# def timing(f):
#     st = time.time()
#     for i in range(100000):
#         f
#     print(time.time()-st)
# timing(np.sqrt(x))
# timing(math.sqrt(x))
# timing(sqrt(x))

print(np.sqrt(x))
print(math.sqrt(x))

print(sqrt(x))