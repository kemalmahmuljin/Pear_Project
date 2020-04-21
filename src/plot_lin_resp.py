import numpy as np
#import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from io import StringIO 
import argparse

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
    #plt.savefig('figure.jpg')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('mesh_to_plot', type=str, nargs='+',
            help=("Number of mesh to plot"))
    
    args = parser.parse_args()
    
    element_file = open("../output/elements_" + args.mesh_to_plot[0], 'r')
    elem_str = element_file.read()
    elem_stream = StringIO(elem_str)
    elements = np.loadtxt(elem_stream)
    
    boundaries_file = open("../output/boundaries_" + args.mesh_to_plot[0], 'r')
    bound_str = boundaries_file.read()
    bound_stream = StringIO(bound_str)
    boundaries = np.loadtxt(bound_stream)
    
    coords_file = open("../output/coords_" + args.mesh_to_plot[0], 'r')
    coords_str = coords_file.read()
    coords_str = coords_str[coords_str.find('\n') + 1:]
    coords_stream = StringIO(coords_str)
    coords = np.loadtxt(coords_stream)
    
    i_coeff_file = open("../output/linear_resp_coeff_" + args.mesh_to_plot[0], 'r')
    i_coeff_str = i_coeff_file.read()
    i_coeff_str = i_coeff_str[i_coeff_str.find('\n') + 1:]
    i_coeff_stream = StringIO(i_coeff_str)
    i_coeff = np.loadtxt(i_coeff_stream)
   
    to_plot = i_coeff
    plot_coeff(coords, to_plot)
