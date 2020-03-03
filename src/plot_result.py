import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D

file_o = open("out", 'r')
lines = file_o.readlines()
counter = 0
condition = True

while condition:
    line = lines[counter]
    if line.strip() != "Node Num - x - y":
        counter += 1
        continue
    counter += 1
    line = lines[counter]
    line_data = line.strip().split()
    coords = np.array([float(line_data[3]), float(line_data[5])])
    counter +=1
    break
while condition:
    line = lines[counter]
    if line.strip() == "":
        counter += 1
        break
    line_data = line.strip().split()
    coords = np.vstack((coords, np.array([float(line_data[3]),
        float(line_data[5])])))
    counter += 1
while condition:
    line = lines[counter]
    if line.strip() != "Element Num - node 1 - node 2 - node 3":
        counter += 1
        continue
    counter += 1
    line = lines[counter]
    line_data = line.strip().split()
    elements = np.array([int(line_data[3]), int(line_data[5]), 
        int(line_data[7])])
    counter +=1
    break
while condition:
    line = lines[counter]
    if line.strip() == "":
        counter += 1
        break
    line_data = line.strip().split()
    elements = np.vstack((elements, np.array([int(line_data[3]), int(line_data[5]), 
        int(line_data[7])])))
    counter += 1
while condition:
    line = lines[counter]
    if line.strip() != "Boundary Num - node 1 - node 2":
        counter += 1
        continue
    counter += 1
    line = lines[counter]
    line_data = line.strip().split()
    boundaries = np.array([int(line_data[3]), int(line_data[5])])
    counter +=1
    break
while condition:
    line = lines[counter]
    if line.strip() == "":
        counter += 1
        break
    line_data = line.strip().split()
    boundaries = np.vstack((boundaries, np.array([int(line_data[3]),
        int(line_data[5])]))) 
    counter += 1
while condition:
    line = lines[counter]
    if line.strip() != "Concentration Coefficients":
        counter += 1
        continue
    counter += 1
    line = lines[counter]
    line_data = line.strip().split()
    line_data[0] = line_data[0].split("[")[1]
    line_data = line_data[:-1]
    coefficients = np.array(list(map(float, line_data)))
    break

x = coords[:,0]
y = coords[:,1]
triang1 = mtri.Triangulation(x,y)
triang2 = mtri.Triangulation(x+60,y)
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
ax.plot_trisurf(triang1, coefficients[42:], cmap='jet')
ax.plot_trisurf(triang2, coefficients[:43], cmap='jet')
ax.view_init(elev=90, azim=-90)
ax.set_xlim(0,120)
ax.set_ylim(0,120)
plt.show()
