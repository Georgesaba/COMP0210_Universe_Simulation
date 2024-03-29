import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("test.txt").T

plt.figure(figsize = (8,6))
plt.title("Comparison of analytic potential vs particle mesh potential")
plt.plot(data[0], label = "potential function")
plt.plot(data[1], label = "expected potential function")
plt.xlabel("X index")
plt.ylabel("Potential Value")
plt.legend(loc="best")
plt.grid(True)
plt.savefig("test_potential/plot.png")