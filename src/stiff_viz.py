import numpy as np
import scipy.sparse as sparse
import matplotlib.pylab as plt
a = np.loadtxt("stiff_2", skiprows=1)
b = sparse.csr_matrix(a)
plt.spy(b, markersize=1)
plt.show()
