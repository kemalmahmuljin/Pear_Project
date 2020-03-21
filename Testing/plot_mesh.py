import matplotlib.pyplot as plt
import numpy as np

file_name = "elements"
coords = np.loadtxt("coords", skiprows=1)
nodes = np.loadtxt("elements").astype(int)

for elem in nodes:
    x = [coords[elem[0]][0], coords[elem[1]][0], coords[elem[2]][0]
            , coords[elem[0]][0]]
    y = [coords[elem[0]][1], coords[elem[1]][1], coords[elem[2]][1]
            , coords[elem[0]][1]]
    plt.plot(x, y)

plt.title("Working importer 23/02/20")
plt.show()

count = 0
for elem in nodes:
    ax = plt.axes()
    for i in range(3):
        ax.arrow(coords[elem[i%3]][0], coords[elem[i%3]][1],
                coords[elem[(i+ 1)%3]][0]-coords[elem[i%3]][0],
                coords[elem[(i + 1)%3]][1]-coords[elem[i%3]][1], 
                head_width=0.001, head_length=0.003, width=0.00001,
                length_includes_head=True)
    plt.xlim(-0.01,0.40)
    plt.ylim(-0.01,0.140)
    plt.show()
    count += 1
