## Fastest Path
import os
import time
import random
import pandas as pd
from copy import copy
import matplotlib.pyplot as plt

'''
This function creates a random array of points that will be used for the pathing function.
It also selects the first and last set of xy coordinates.
The commented section allows selecting your own points.
Currently, since this is random.seeded, it will always create the same points (depending on how many points selected)
'''
def point_locations(points):
    def rand_xy():
        x = random.randint(0,10)
        y = random.randint(0,10)
        return x,y

    xs = []
    ys = []
    xy = []
    idx = 0
    while idx < points:
        x,y = rand_xy()
        for i in range(len(xs)):
            if x == xs[i] and y == ys[i]:
                x,y = rand_xy()
                continue
        xy.append([x,y])
        idx += 1
    # # Uncomment this section to make your own points
    # xy = [
    #     [0,0],[4,0],[4,4],[0,4],[1,1],[3,3],[2,2]
    # ]
    # #
    xy1 = xy[0]
    xyf = xy[-1]
    return xy,xy1,xyf

'''
This function finds the distance between all the points in the list.
It adds the total distance to a sum to find the overall total distance
'''
def distance(xy):
        x_dist = 0
        y_dist = 0
        total_dist = (x_dist**2 + y_dist**2)**.5
        for idx in range(len(xy)-1):
            x_dist = xy[idx+1][0] - xy[idx][0]
            y_dist = xy[idx+1][1] - xy[idx][1]
            total_dist += (x_dist**2 + y_dist**2)**.5
            # print(total_dist)
        return total_dist

'''
The pathing function randomizes the order of elements in the list.
This method of path finding is very inefficient, but for small amounts of points (< 25) it returns a solution relatively quick.

A copy of the xy list is made so we can delete elements from the copied list. If we only referenced it, deleting would globally delete.
a new xy list is made so we can add the first and last elements so they are not randomized and are always in the same location
'''
def pathing(xy,xy1,xyf):
    old_xy = copy(xy[1:-2])
    new_xy = [xy1]
    while True:
        rand_idx = random.randint(0,len(old_xy)-1)    
        new_xy.append(old_xy[rand_idx])
        del old_xy[rand_idx]
        if len(old_xy) == 0:
            break
    new_xy.append(xyf)
    dist = distance(new_xy)
    return new_xy, dist

'''
The plotting function plots (duh!)
It uses list comprehension to separate the x's and y's from the new_xy list.
Other than that it plots normally, using scatter and plot to plot dots as well as a semi transparent line.
'''
def plotting(new_xy,total_dist):
    xs = [x[0] for x in new_xy]
    ys = [y[1] for y in new_xy]
    
    plt.cla()
    plt.scatter(xs[1:-1],ys[1:-1],c='k')
    plt.scatter(xs[0],ys[0],c='g')
    plt.scatter(xs[-1],ys[-1],c='r')
    plt.plot(xs,ys,c='b',alpha=.2)
    plt.title('Path Length: ' + str(round(total_dist,3)))
    plt.draw()
    plt.pause(.5)

'''
Main function uses all the functions described above.

It also has a while loop that allows us to repeat the pathing scripts points^4 times. 
I chose the points^4 number somewhat arbitrarily, but from values ranging from 5-25 it works and finds a reasonable solution.

In the while loop some lists are created including ALL xy's, ALL distances, and the xy's in string form.
I selected these all and put them in a pandas dataframe in order to use the sort_values() and drop_duplicates() function.
The drop_dupes() function does not work on lists, so that is why the "xy string" is included in the dataframe.

Finally we plot the 20 lowest values using plt.pause() to hold on each plot for half a second.
'''
def main():
    start_time = time.time()

    points = 10
    random.seed(0)
    xsys,x1y1,xfyf = point_locations(points)

    loops = 0
    full_xy = []
    full_xy_str = []
    full_dist = []
    while loops < points**4:
        new_xy,dist_xy = pathing(xsys,x1y1,xfyf)
        loops += 1
        full_xy.append(new_xy)
        full_xy_str.append(str(new_xy))
        full_dist.append(dist_xy)
    print(loops)
    df = pd.DataFrame()
    df["distance"] = full_dist
    df["xy"] = full_xy
    df["xy_str"] = full_xy_str

    df = df.sort_values(by=["distance"],ascending=False).drop_duplicates(subset=["xy_str"]).reset_index(drop=True)
    print(df)
    print("Ran in {} seconds".format(round(time.time()-start_time,3)))

    for i in df.index[-20:]:
        plotting(df.iloc[i]["xy"],df.iloc[i]["distance"])

'''
Not really useful in this case as I wont run this main script outside this file, but I have it ingrained in me to always use this format to run main()
'''

if __name__ == "__main__":
    os.system('cls')
    main()
    plt.show()
