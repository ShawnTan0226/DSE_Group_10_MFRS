import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import integrate

n = 100 # resolution of data along half span

Sweep05 = np.deg2rad(35.36)

# Weights in the aircraft
BatteryMass = 2903 + 10600 # kg from 2nd weight estimation
StructuralMass = 1500 + 2813 + 180 + 558 + 126 # kg from 2nd weight estimation
EngineMass = 1224 # kg for 7 engines, evenly spread
EngineMass = EngineMass / 7
ThrustTO = 6000 # N thrust of one engine at Lift Off
ThrustCruise = 4288.32 # N thrust of one engine at cruise
# Take off force at hook
acceleration = 2.84

mass_total = BatteryMass + StructuralMass + EngineMass
ForceHookTO = acceleration * mass_total - ThrustTO * 7
print(ForceHookTO)

##### Import planform and aerodynamic data #####

def importdat(n):
    CruiseDistribution = 'MainWing_a=6.25_v=110.00ms.csv'
    TODistribution = 'MainWing_a=18.00_v=63.60ms.csv'
    data_cruise = pd.read_csv(CruiseDistribution, index_col=None)
    data_TO = pd.read_csv(TODistribution, index_col=None)
    data_cruise = data_cruise.to_dict(orient='list')
    data_TO = data_TO.to_dict(orient='list')
    y_old =  data_cruise['y-span']

    # Interpolating to n data points with a new y axis

    y_points = np.linspace(data_cruise['y-span'][0], data_cruise['y-span'][-1], n)
    for i in data_cruise.keys():
        data_cruise[i] = np.interp(y_points, y_old, data_cruise[i])
    for i in data_TO.keys():
        data_TO[i] = np.interp(y_points, y_old, data_TO[i])

    # Test plots
    #plt.plot(y_points,data_cruise['Chord'])
    #plt.show()

    return data_cruise, data_TO, y_points

##### Define the force distributions along the beam (half span) #####

def forces(CruiseDistribution, TODistribution, y_points, n):

    # Lift and drag at cruise :
    Cl = CruiseDistribution['Cl']
    Cdi = CruiseDistribution['ICd']
    c = CruiseDistribution['Chord']
    V = 110 # m/s
    rho = 0.7361 # kg/m^3
    dy = (CruiseDistribution['y-span'][-1] - CruiseDistribution['y-span'][0]) / n
    Lcruise = Cl * 0.5 * rho * V**2 * c # N/m
    Dcruise = Cdi * 0.5 * rho * V**2 * c # N/m

    # Lift and drag at Lift off :
    Cl = TODistribution['Cl']
    Cdi = TODistribution['ICd']
    c = TODistribution['Chord']
    V = 63.6  # m/s
    rho = 1.225  # kg/m^3
    dy = (TODistribution['y-span'][-1] - TODistribution['y-span'][0]) / n
    Ltakeoff = Cl * 0.5 * rho * V**2 * c # N/m
    Dtakeoff = Cdi * 0.5 * rho * V**2 * c  # N/m


    # Weight distribution along span
    # distributed weights of batteries, payload and structure
    Wave = (StructuralMass + BatteryMass)/2 * 9.81 / TODistribution['y-span'][-1] # N/m
    W = Wave * TODistribution['Chord'] / ((TODistribution['Chord'][0] + TODistribution['Chord'][-1]) / 2) # N/m mass spread normalised with chord distribution.
    W[0] += EngineMass * 9.81 / (CruiseDistribution['y-span'][-1] / n)
    W[-1] += EngineMass*2 * 9.81 / (CruiseDistribution['y-span'][-1] / n)
    W[int(round(n/3))] += EngineMass * 9.81 / (CruiseDistribution['y-span'][-1] / n)
    W[int(round(n*2/3))] += EngineMass * 9.81 / (CruiseDistribution['y-span'][-1] / n)

    # Lower weight:
    #W = W * 0.3

    # Drag force distribution :
    # Thrust force :
    T = np.zeros(n)
    T[0] += ThrustTO + ForceHookTO
    T[-1] += ThrustTO * 2
    T[int(round(n/3))] += ThrustTO
    T[int(round(n*2/3))] += ThrustTO
    TTO = T / (CruiseDistribution['y-span'][-1] / n)

    T = np.zeros(n)
    T[0] += ThrustCruise
    T[-1] += ThrustCruise * 2
    T[int(round(n / 3))] += ThrustCruise
    T[int(round(n * 2 / 3))] += ThrustCruise
    TC = T / (CruiseDistribution['y-span'][-1] / n)

    # Moment distribution :
    V = 110 # m/s
    rho = 0.7361 # kg/m^3
    M_cruise = CruiseDistribution['CmAirf@chord/4'] * 1/2 * rho * V**2 * CruiseDistribution['Chord']

    V = 63.6  # m/s
    rho = 1.225  # kg/m^3
    M_TO = TODistribution['CmAirf@chord/4'] * 1 / 2 * rho * V**2 * TODistribution['Chord']

    return Lcruise, Ltakeoff, TC, TTO, W, M_cruise, M_TO, Dcruise, Dtakeoff

def TestForces(Lcruise, Ltakeoff, Tcruise, W, TODistribution):
    dy = (TODistribution['y-span'][-1] - TODistribution['y-span'][0]) / n
    c = TODistribution['Chord']
    #print('Ltot = ',np.trapz(Lcruise/9.81, TODistribution['y-span']),' kg')
    #print('Wtot = ',np.trapz(W/9.81,TODistribution['y-span']),' kg')

def InternalLoads(L, T, W, D, M, n, y_points, halfspan,sweep):
    b = halfspan * 2
    Dtot = T - D # drag and thrust act on the x axis
    Vx = integrate.cumtrapz(np.flip(Dtot * b / (2 * n)))[::-1]
    Vy = integrate.cumtrapz(np.flip((-L + W) * b / (2 * n)))[::-1]
    Vx = np.append(Vx,[0])
    Vy = np.append(Vy,[0])

    My = -integrate.cumtrapz(np.flip(Vx * b / (2 * n)))[::-1]
    Mx = -integrate.cumtrapz(np.flip(Vy * b / (2 * n)))[::-1]
    My = np.append(My,[0])
    Mx = np.append(Mx,[0])
    #Mz = integrate.cumtrapz(np.flip(((L-W) * np.tan(Sweep05) * y_points - M) * b / (2 * n)))[::-1
    # Internal torque calculation: first due to the sweep angle then added on due to aerodynamic moment
    Ml = []
    Mw =[]
    yp = y_points
    count = 2

    # due to lift and weight moment arm
    for i in yp[1:]:
        Mli = IntegrateTorqueFromLift(count, yp, -L, sweep)
        Mwi = IntegrateTorqueFromLift(count, yp, W, sweep)
        Ml.append(Mli)
        Mw.append(Mwi)
        count += 1

    # due to aerodynamic moment:
    Mzm = integrate.cumtrapz(np.flip(M), np.flip(yp))
    Mz = (np.array(Ml) + np.array(Mw) + np.array(Mzm))[::-1]
    Mz = np.append(Mz, [0])

    #plt.plot(halfspan, np.array(Mw)/200)
    #plt.plot(halfspan, np.array(Ml)/200 )
    #plt.plot(data['y-span'][:-1], np.array(Ml) + np.array(Mw))
    # plt.plot(data['y-span'], )
    #plt.plot(data['y-span'], My)
    #plt.show()

    return Vx, Vy, Mx, My, Mz

def IntegrateTorqueFromLift(c, axis, data, sweep):
    data = np.flip(data)
    data = data[:c]
    axis = np.flip(axis)
    axis = axis[:c]
    data = data * np.tan(sweep) * (axis - axis[c-1])

    M = np.trapz(data, axis)
    return M

def PrelimSizing(Vx, Vy, Mx, My, Mz, chords):
    # Material properties
    sigma = 324000000 # Pa yield strength 2024 Al
    tau = 283000000 # Pa yield shear strength 2024 Al
    tc = 0.15 # thickness over chord ratio of airfoil
    t = chords * tc

    Ixx = abs(Mx) * t / (2 * sigma)
    Iyy = abs(My) * chords / (2 * sigma)
    # assuming two area points at the thickness extremities:
    A = Ixx * 4 / t ** 2
    # determining the necessary spar thickness for shear stresses due to internal shear force:
    Q = 2 * A * 0.5 * t
    ts = abs(Vy) * Q / (Iyy * tau) # spar thickness needed to withstand shear stress due to bending

    # thickness of wingbox for torque
    twb = abs(Mz)/(2 * chords**2 * (0.4*0.15) * tau)

    #plt.plot(y_points, A)
    #plt.plot(y_points, Ixx)
    #plt.show()
    return Ixx, Iyy, ts, twb

nl = -1 # loadfactor
CruiseDistribution, TODistribution, y_points = importdat(n)
Lcruise, Ltakeoff, TC, TTO, W, M_cruise, M_TO, Dcruise, Dtakeoff = forces(CruiseDistribution, TODistribution, y_points, n)
Vx, Vy, Mx, My, Mz = InternalLoads(nl*Lcruise, TC, W, abs(nl)*Dcruise, nl*M_cruise, n, y_points, CruiseDistribution['y-span'][-1], Sweep05)
#plt.plot(y_points, Vy, label='Internal shear force on the y axis, take off')
Ixx, Iyy, ts, twb = PrelimSizing(Vx, Vy, Mx, My, Mz, CruiseDistribution['Chord'])

#print(M_cruise)
# TestForces(Lcruise, Ltakeoff, Tcruise, W, TODistribution)

#plt.xlabel("z")
#plt.ylabel("Vy [N]")
#plt.grid()
#plt.legend(
#plt.show()

def plots():
    plt.plot(y_points, ts, label='ts')
    plt.plot(y_points, twb, label='twb')
    #plt.plot(y_points, -W, label='Ixx')
    #plt.plot(y_points, Mx, label='Mx')
    #plt.plot(y_points, , label='tot')
    #plt.plot(y_points, My, label='My')
    #plt.plot(y_points, Mz, label='Mz')
    plt.xlabel("z [m]")
    plt.ylabel("Ixx [m^4]")
    plt.grid()
    plt.legend()
    plt.show()

plots()