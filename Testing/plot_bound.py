import matplotlib.pyplot as plt
import numpy as np
from io import StringIO 

coords_file = open("coords", 'r')
coords_str = coords_file.read()
coords_str = coords_str[coords_str.find('\n') + 1:]
coords_stream = StringIO(coords_str)
coords = np.loadtxt(coords_stream)

boundaries_file = open("boundaries", 'r')
bound_str = boundaries_file.read()
bound_stream = StringIO(bound_str)
boundaries = np.loadtxt(bound_stream)

ax = plt.axes()
for elem in boundaries:
    x = [coords[int(elem[0])][0], coords[int(elem[1])][0]]
    y = [coords[int(elem[0])][1], coords[int(elem[1])][1]]
    if elem[2] == 0:
        color = 'r'
    else:
        color = 'b'
    ax.arrow(coords[int(elem[0])][0],coords[int(elem[0])][1],
            coords[int(elem[1])][0]-coords[int(elem[0])][0],
            coords[int(elem[1])][1]-coords[int(elem[0])][1], 
            head_width=0.002, length_includes_head=True, lw=0.00001, 
            color=color)
plt.xlim(0,0.050)
plt.ylim(0,0.120)    

plt.title("Working boundaries 05/03/20")
plt.show()
