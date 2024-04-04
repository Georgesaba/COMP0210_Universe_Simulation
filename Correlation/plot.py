import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

filename = "Comparison_4_1._1.04"

data = pd.read_csv("Correlation/" + filename + ".csv").iloc[1:, :] #due to inf value
step = 0.5/len(data)
x = np.arange(0, 0.5, step)
x += (x[1] - x[0])/2

plt.figure()
plt.title("Correlations for Different Expansion Factors")
plt.xlabel("Radial Bin")
plt.ylabel("Radial Correlation Function")
for i in data.columns:
   plt.plot(x,data[i],label = f"Expansion Factor: {i}")
plt.legend(loc='best')
plt.grid(True)
plt.savefig("Correlation/" + filename + ".png")