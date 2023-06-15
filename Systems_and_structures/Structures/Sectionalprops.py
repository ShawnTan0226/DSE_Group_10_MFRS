import numpy as np
from matplotlib import pyplot as plt
import numpy as np

class Sectional_Properties:
    def __init__(self):
        self.spar_pos = np.array([0.15, 0.55])
        self.Ixx_stringer = 0.0000000001
        self.Izz_stringer = 0.0000000001
        self.run()


    def DefineAirfoil(self):
        airfoil = {'x': [], 'y': []}

        with open('MH 91  14.98%.dat', 'r') as file:
             for line in file:
                values = line.split()
                airfoil['x'].append(float(values[0]))
                airfoil['y'].append(float(values[1]))
        return airfoil

    def CellGen(self):
        A_cell = np.array([0])
        for x_spar in self.spar_pos:
            y_spar = np.interp(x_spar,self.airfoil['x'],self.airfoil['y'])
            A_this_cell = np.trapz(self.airfoil['y'])
            A_cell = np.concatenate((A_cell, np.array([A_this_cell])))
        return A_cell

    def run(self):
        self.airfoil = self.DefineAirfoil()
        self.a_cell = self.CellGen()
        plt.plot(self.airfoil['x'], self.airfoil['y'])
        plt.plot([0.15,0.15],[-0.08,0.08])
        plt.plot([0.55, 0.55], [-0.08, 0.08])

        plt.axis('equal')
        plt.show()

SectionalProperties = Sectional_Properties()
SectionalProperties.run()