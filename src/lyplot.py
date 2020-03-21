import numpy as np
import scipy as sp
from io import StringIO 
######
import plotly.offline as py
import plotly.figure_factory as FF

######


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

fcoeff_file = open("final_coeff", 'r')
fcoeff_str = fcoeff_file.read().strip()
fcoeff_str = fcoeff_str[fcoeff_str.find('[') + 1:-2]
fcoeff_stream = StringIO(fcoeff_str)
fcoeff = np.loadtxt(fcoeff_stream)


f_file = open("f_vector", 'r')
f_str = f_file.read().strip()
f_str = f_str[f_str.find('[') + 1:-2]
f_stream = StringIO(f_str)
f_vector = np.loadtxt(f_stream)

stiff = np.loadtxt("stiff_2", skiprows=1)

def plott(elements,coords, coeff):

	#fig = make_subplots(rows=2, cols=2)
	
		
	fig_lin_CO2  = FF.create_trisurf(x=coords[:,0], y=coords[:,1], z=coeff[coeff.shape[0]//2:] ,simplices=elements.astype(int),title="Linearized CO2" )
	fig_lin_O2 = FF.create_trisurf(x=coords[:,0], y=coords[:,1], z=coeff[:coeff.shape[0]//2] ,simplices=elements.astype(int),title="Linearized O2" )

	fig_fin_O2 = FF.create_trisurf(x=coords[:,0], y=coords[:,1], z=fcoeff[fcoeff.shape[0]//2:] ,simplices=elements.astype(int),title="Final O2" )
	fig_fin_CO2 = FF.create_trisurf(x=coords[:,0], y=coords[:,1], z=fcoeff[:fcoeff.shape[0]//2] ,simplices=elements.astype(int) ,title="Final CO2" )


	fig_lin_CO2.show()
	fig_lin_O2.show()
	fig_fin_CO2.show()
	fig_fin_O2.show()

if __name__ == "__main__":

    plott(elements, coords, coeff)
    #plot_coeff(coords, fcoeff)
