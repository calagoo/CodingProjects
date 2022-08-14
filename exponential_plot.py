import matplotlib.pyplot as plt
import numpy as np

w = 96
h = 32

xc = w/2
yc = h/2

x0 = int(0*w/3)
xf = int(1*w/3)
y0 = int(1)
yf = int(h-1)

y_l = []
y_l2 = []
steps = .02

for x in np.arange(x0,xf,steps):
    res = h*((x-x0)**2/(xf-x0)**2)
    res2 = h*((x-x0)**2/(xf-x0)**2)+10
    y_l.append(res)
    y_l2.append(res2)
    # print(x,res)
    # print(x,res2)
x = np.arange(x0,xf,steps)

plt.scatter(w-x,h-np.array(y_l),s=10)
plt.scatter((w-x),h-np.array(y_l2),s=10)

plt.show()