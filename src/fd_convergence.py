import numpy as np
import matplotlib.pyplot as plt

fd_data = np.loadtxt('./fd_convergence')
plt.plot(fd_data[:, 0], np.sqrt(fd_data[:, 1]))
plt.xscale('log')
plt.yscale('log')
plt.title("Finite Difference Convergence")
plt.xlabel(r'$Finite \; Differences \; \epsilon $')
plt.ylabel(r'$\frac{\|J_{fd} - J_{a}\|_2}{\|J_{a}\|_2}$')
plt.grid()
plt.tight_layout()
plt.show()
