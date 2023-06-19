import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

CP_file = pd.read_csv("./Pressure distribution.csv", header=None)

x = CP_file.iloc[:,0]
CP = CP_file.iloc[:,1]
print(CP)
CP_array = np.array(CP)
print(CP_array)

# p_tot = 0.5 * CP * 0.770816 * 110^2 + 57182.0
# print()

def total_pressure(V, C_P, rho_inf):
    p_tot = 0.5 * C_P * rho_inf * V**2
    return p_tot

y = total_pressure(110, CP, 0.770816)



fig, ax = plt.subplots()

# Possible colours: 1) blue, 2) red, 3) green, 4) yellow, 5)
ax.plot(x, y, linestyle='-', label="Label", color='blue')

ax.set_xlabel('Location [x/c]')
ax.set_ylabel('Pressure [Pa]')
ax.minorticks_on()
ax.grid(which='major', color='grey', linestyle='-')
ax.grid(which='minor', color='lightgrey', linestyle='dotted')
ax.axhline(y=0, xmin=0, xmax=1, linewidth=2, color='black')
# ax.invert_yaxis()
plt.gca().invert_yaxis()
#ax.legend()
plt.title("Pressure distribution")
plt.show()