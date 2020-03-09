import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from io import StringIO 

element_file = open("elements", 'r')
elem_str = element_file.read()
elem_stream = StringIO(elem_str)
elements = np.loadtxt(elem_stream)

boundaries_file = open("boundaries", 'r')
bound_str = boundaries_file.read()
bound_stream = StringIO(bound_str)
boundaries = np.loadtxt(bound_stream)

coords_file = open("coords", 'r')
coords_str = coords_file.read()
coords_str = coords_str[coords_str.find('\n') + 1:]
coords_stream = StringIO(coords_str)
coords = np.loadtxt(coords_stream)

coeff_file = open("initial_coeff", 'r')
coeff_str = coeff_file.read().strip()
coeff_str = coeff_str[coeff_str.find('[') + 1:-2]
coeff_stream = StringIO(coeff_str)
coeff = np.loadtxt(coeff_stream)

f_file = open("f_vector", 'r')
f_str = f_file.read().strip()
f_str = f_str[f_str.find('[') + 1:-2]
f_stream = StringIO(f_str)
f_vector = np.loadtxt(f_stream)

stiff = np.loadtxt("stiff_2", skiprows=1)
def plot_coeff(coords, coeff):
    x = coords[:,0]
    y = coords[:,1]
    triang1 = mtri.Triangulation(x,y)
    triang2 = mtri.Triangulation(x+0.060,y)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.plot_trisurf(triang1, coeff[coeff.shape[0]//2:], cmap='jet')
    ax.plot_trisurf(triang2, coeff[:coeff.shape[0]//2+1], cmap='jet')
    ax.view_init(elev=90, azim=-90)
    ax.set_xlim(0,0.120)
    ax.set_ylim(0,0.120)
    plt.show()

if __name__ == "__main__":
    plot_coeff(coords, coeff)
