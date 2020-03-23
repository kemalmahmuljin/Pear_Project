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

coeff_file = open("final_coeff", 'r')
coeff_str = coeff_file.read()
coeff_str = coeff_str[coeff_str.find('\n') + 1:]
coeff_stream = StringIO(coeff_str)
coeff = np.loadtxt(coeff_stream)

f_file = open("f_vector", 'r')
f_str = f_file.read()
f_str = f_str[f_str.find('\n') + 1:]
f_stream = StringIO(f_str)
f_vector = np.loadtxt(f_stream)

f2_file = open("f_vector_lin", 'r')
f2_str = f2_file.read()
f2_str = f2_str[f2_str.find('\n') + 1:]
f2_stream = StringIO(f2_str)
f2_vector = np.loadtxt(f2_stream)

#stiff = np.loadtxt("stiff_2", skiprows=1)
def plot_coeff(coords, coeff):
    x = coords[:,0]
    y = coords[:,1]
    triang1 = mtri.Triangulation(x,y, elements.astype(int))
    triang2 = mtri.Triangulation(x+0.060,y, elements.astype(int))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.plot_trisurf(triang1, coeff[coeff.shape[0]//2:], cmap='jet')
    ax.plot_trisurf(triang2, coeff[:coeff.shape[0]//2], cmap='jet')
    ax.view_init(elev=90, azim=-90)
    ax.set_xlim(0,0.120)
    ax.set_ylim(0,0.120)
    plt.show()

if __name__ == "__main__":
    plot_coeff(coords, coeff)
