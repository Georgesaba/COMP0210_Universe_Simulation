import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

filename = "Comparison_4_1.03_1.04"

data = pd.read_csv("Correlation/" + filename + ".csv").iloc[1:, :] #due to inf value

plt.figure()
plt.title("Correlation/Correlation Plots for Different Expansion Factors")
for i in data.columns:
   plt.plot(data[i],label = f"Expansion Factor: {i}")
plt.legend(loc='best')
plt.grid(True)
plt.savefig("Correlation/" + filename + ".png")