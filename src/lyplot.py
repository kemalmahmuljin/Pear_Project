import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from io import StringIO 
import argparse
import plotly.offline as py
import plotly.figure_factory as FF

#stiff = np.loadtxt("stiff_2", skiprows=1)
def plot_coeff(elements,coords, coeff):	
	fig_lin_CO2  = FF.create_trisurf(x=coords[:,0], y=coords[:,1], z=coeff[coeff.shape[0]//2:] ,simplices=elements.astype(int),title="Linearized CO2" )
	fig_lin_O2 = FF.create_trisurf(x=coords[:,0], y=coords[:,1], z=coeff[:coeff.shape[0]//2] ,simplices=elements.astype(int),title="Linearized O2" )

	fig_lin_CO2.show()
	fig_lin_O2.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('function_to_plot', type=str, nargs='+',
            help=("Function to plot: \n initial_c: Initial coefficients \n final_c: "
            " Final coefficients \n  f_vector: Constant vector \n f_vector_lin: "
            " Constant vector + linearized term \n f_non_lin: Linearized non "
            " linear addition to f vector \n MC_lin: Linearized Matlab coefficient"))
    
    args = parser.parse_args()
    
    element_file = open("../output/elements", 'r')
    elem_str = element_file.read()
    elem_stream = StringIO(elem_str)
    elements = np.loadtxt(elem_stream)
    
    boundaries_file = open("../output/boundaries", 'r')
    bound_str = boundaries_file.read()
    bound_stream = StringIO(bound_str)
    boundaries = np.loadtxt(bound_stream)
    
    coords_file = open("../output/coords", 'r')
    coords_str = coords_file.read()
    coords_str = coords_str[coords_str.find('\n') + 1:]
    coords_stream = StringIO(coords_str)
    coords = np.loadtxt(coords_stream)
    
    f_coeff_file = open("../output/final_coeff", 'r')
    f_coeff_str = f_coeff_file.read()
    f_coeff_str = f_coeff_str[f_coeff_str.find('\n') + 1:]
    f_coeff_stream = StringIO(f_coeff_str)
    f_coeff = np.loadtxt(f_coeff_stream)
    
    i_coeff_file = open("../output/initial_coeff", 'r')
    i_coeff_str = i_coeff_file.read()
    i_coeff_str = i_coeff_str[i_coeff_str.find('\n') + 1:]
    i_coeff_stream = StringIO(i_coeff_str)
    i_coeff = np.loadtxt(i_coeff_stream)
   
    f_file = open("../output/f_vector", 'r')
    f_str = f_file.read()
    f_str = f_str[f_str.find('\n') + 1:]
    f_stream = StringIO(f_str)
    f_vector = np.loadtxt(f_stream)
    
    f2_file = open("../output/f_vector_lin", 'r')
    f2_str = f2_file.read()
    f2_str = f2_str[f2_str.find('\n') + 1:]
    f2_stream = StringIO(f2_str)
    f2_vector = np.loadtxt(f2_stream)
    
    #cc_vector = np.loadtxt("../output/calculated_coeff")
    
    #fc_vector = np.loadtxt("../output/f_calc")

    MC_file = open("../output/MC_lin", 'r')
    MC_str = MC_file.read()
    MC_str = MC_str[MC_str.find('\n') + 1:]
    MC_stream = StringIO(MC_str)
    MC_vector = np.loadtxt(MC_stream)

    cond = args.function_to_plot[0]
    print(cond)
    if cond == "initial_c":
        to_plot = i_coeff
    elif cond == "final_c":
        to_plot = f_coeff
    elif cond == "f_vector":
        to_plot = f_vector
    elif cond == "f_vector_lin":
        to_plot = f2_vector
    elif cond == "f_non_lin":
        to_plot = f2_vector-f_vector
    elif cond == "f_calc":
        to_plot = fc_vector
    elif cond == "calc_coef":
        to_plot = cc_vector
    elif cond == "MC_lin":
        to_plot = MC_vector 
    else:
        print("Error in the input argument use -h for help")
        exit()
    plot_coeff(elements, coords, to_plot)
