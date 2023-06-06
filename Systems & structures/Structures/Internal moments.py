########## Determine internal moments and shear forces ##########
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import take-off and cruise force distributions
CruiseDistribution = 'MainWing_a=6.25_v=110.00ms.csv'
TODistribution = 'MainWing_a=18.00_v=63.60ms.csv'

data_cruise = pd.read_csv(CruiseDistribution, index_col=None)
data_TO = pd.read_csv(TODistribution, index_col=None)

data_cruise = data_cruise.to_dict(orient='list')
data_TO = data_TO.to_dict(orient='list')

# print(data_TO['y-span'])

# Aircraft parameters
Sweep05 = np.deg2rad(35.36)
b = data_cruise['y-span'][-1] * 2 # span in meters
n = 100 # amount of elements
tc = np.ones(n) * 0.15 # thickness to chord ratio of airfoil
S = 81/2 # m^2 surface
T = [1000,1000] # Thrust of engines in Newton
T = np.array(T) / ( b /( 2 * n)) # normalised point Thrust in N/m
ZT = [2,3,4,5,6] # position of engines along z in meters

PointWeights = [1000,1000] # point weights in Newton
PointWeights = np.array(PointWeights) /( b /( 2 * n)) # normalised point weights in N/m

ZPointWeights = [3,1] # point weight positions in meters
StrucWeight = 7000 # Newtons/m

y_points = np.linspace(data_cruise['y-span'][0],data_cruise['y-span'][-1],n)
Chord_mapped = np.interp(y_points, data_cruise['y-span'] , data_cruise['Chord'])

def GenerateLoadingCruise(n, b, T, ZT, StrucWeight, ZPointWeights, PointWeights, data_cruise):

    y_points = np.linspace(data_cruise['y-span'][0],data_cruise['y-span'][-1],n) # data points position along y
    Cl_mapped = np.interp(y_points, data_cruise['y-span'] , data_cruise['Cl'])
    Chord_mapped = np.interp(y_points, data_cruise['y-span'] , data_cruise['Chord'])

    L = Cl_mapped * Chord_mapped * b/(2 *n) * 0.5 * 1.225 * 110**2

    StrucWeightDistribution = np.ones(n) * StrucWeight

    D = 0.1 * L

    for i in range(len(T)):
        posT = ZT[i]
        xp = int(np.round(posT / (b/2) * n, 0))
        D[xp] = D[xp] - T[i]

    W = StrucWeightDistribution

    for i in range(len(ZPointWeights)):
        ZPtWt = ZPointWeights[i]
        xp = int(np.round(ZPtWt / (b/2) * n, 0))
        W[xp] += PointWeights[i]

    M = np.zeros(n) * 1000
    return D,L,W,M, y_points

D,L,W,M, y_points = GenerateLoadingCruise(n, b, T, ZT, StrucWeight, ZPointWeights, PointWeights, data_cruise)

def InternalLoads(D,L,W,M,b,n,Sweep05):
    Vx = np.cumsum(np.flip(D*b/(2*n)))[::-1]
    Vy = np.cumsum(np.flip((-L+W)*b/(2*n)))[::-1]


    cumload = np.cumsum(np.flip((L-W)* np.tan(Sweep05)*b/(2*n)))[::-1]
    shear = np.cumsum(np.flip(cumload*b/(2*n)))[::-1]
    Mz = np.cumsum(np.flip(M*b/(2*n)))[::-1] + np.cumsum(np.flip(shear*b/(2*n)))[::-1]

    cumload = np.cumsum(np.flip((-L+W)*b/(2*n)))[::-1]
    shear = np.cumsum(np.flip(cumload * b / (2 * n)))[::-1]
    Mx = np.cumsum(shear * b / (2 * n))[::-1]

    cumload = np.cumsum(np.flip((-D)* b / (2 * n)))[::-1]
    shear = np.cumsum(np.flip(cumload * b / (2 * n)))[::-1]
    My = np.cumsum(np.flip(shear * b / (2 * n)))[::-1]
    return Vx, Vy, Mz, Mx, My

Vx, Vy, Mz, Mx, My = InternalLoads(D,L,W,M,b,n,Sweep05)

t = tc * np.interp(y_points, data_cruise['y-span'] , data_cruise['Chord']) # thickness of airfoil distribution along span in meters
sigma_max = 324000000 # Pa yield strength 2024 Al
tau_max = 283000000 # Pa yield shear strength 2024 Al

def necessary_values(sigma_max, tau_max, t, Vx, Vy, Mz, Mx, My):
    Iyy = Mx * t / (2 * sigma_max) # necessary second moment of area to withstand normal forces a
    # assuming two area points at the thickness extremities:
    A = Iyy * 4 / t**2
    # determining the necessary spar thickness for shear stresses due to internal shear force:
    Q = 2 * A * 0.5 * t
    ts = Vy * Q / (Iyy * tau_max) # spar thickness needed to withstand shear stress due to bending
    return Iyy, ts

Iyy, ts = necessary_values(sigma_max, tau_max, t, Vx, Vy, Mz, Mx, My)

plt.plot(y_points, L-W)
plt.show()