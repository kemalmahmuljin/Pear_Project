import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../output/lin_mod_runtime")
plt.plot(data[:,0], data[:,1])
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.title("Runtime for Linear Respiration")
plt.xlabel('Number of Nodes []')
plt.ylabel('Time [s]')
plt.xlim(1e2, 1e4)
plt.show()
