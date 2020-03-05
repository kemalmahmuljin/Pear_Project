import matplotlib.pyplot as plt
import numpy as np

file_name = "mesh_test_out"
file_o = open(file_name, 'r')
fig = plt.figure()
line = file_o.readline()
line_data = file_o.readline().strip().split()
coords = np.array([float(line_data[3]), float(line_data[5])])
line = file_o.readline()
while line.strip() != "":
    line_data = line.strip().split()
    coords = np.vstack((coords, 
        np.array([float(line_data[3]), float(line_data[5])])))
    line = file_o.readline()

line = file_o.readline()
line_data = line.strip().split()
nodes = np.array([int(line_data[3]), int(line_data[5]), int(line_data[7])])
line = file_o.readline()
while line.strip() != "":
    line_data = line.strip().split()
    nodes = np.vstack((nodes, 
        np.array([int(line_data[3]), int(line_data[5]), int(line_data[7])])))
    line = file_o.readline()

for elem in nodes:
    x = [coords[elem[0]][0], coords[elem[1]][0], coords[elem[2]][0]
            , coords[elem[0]][0]]
    y = [coords[elem[0]][1], coords[elem[1]][1], coords[elem[2]][1]
            , coords[elem[0]][1]]
    plt.plot(x, y)

plt.title("Working importer 23/02/20")
plt.show()
