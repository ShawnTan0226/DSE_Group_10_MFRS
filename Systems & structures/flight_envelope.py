import numpy as np
import matplotlib.pyplot as plt

class Flight_Envelope:
    def __init__(self, MTOW, CL_max, v_cruise, v_dive, rho, wing_loading):
        self.MTOW = MTOW
        self.CL_max = CL_max
        self.v_cruise = v_cruise
        self.v_dive = v_dive
        self.rho = rho
        self.wing_loading = wing_loading

        # Set array of velocities
        self.v_array = np.arange(0, v_dive, 10)

    def manoeuver_loads(self):
        # Positive part
        # Initialise load factor array
        # self.manoeuver_array = []
        # Define max manoeuver load factor from STANAG 4671
        self.max_load_factor = min(2.1 + 10900/(self.MTOW+4536), 3.8)
        # Calculate load factors up until max manoeuver load
        n = (0.5 * self.rho * np.square(self.v_array) * self.CL_max) / self.wing_loading
        self.manoeuver_array = np.clip(n, 0, self.max_load_factor)
        # self.manoeuver_array.append(n)

        # Negative part
        # Initialise load factor array
        # self.neg_manoeuver_array = []
        # Define max negative manoeuver load factor from STANAG 4671
        self.min_load_factor = -0.4 * self.max_load_factor
        self.v_array_negative = self.v_array[self.v_array <= self.v_cruise]
        # Calculate load factors up until max negative manoeuver load
        n_neg = -(0.5 * self.rho * np.square(self.v_array_negative) * self.CL_max) / self.wing_loading
        self.neg_manoeuver_array = np.clip(n_neg, 0, self.min_load_factor)
        # self.neg_manoeuver_array.append(n_neg)

    def plot(self):
        fig, ax = plt.subplots()

        ax.plot(self.v_array, self.manoeuver_array, linestyle='-',  label="Yeh")
        ax.plot(self.v_array_negative, self.neg_manoeuver_array, linestyle='-', label="Faka")

        ax.set_xlabel('Velocity (m/s)')
        ax.set_ylabel('Load factor')
        ax.minorticks_on()
        ax.grid(which='major', color='grey', linestyle='-')
        ax.grid(which='minor', color='lightgrey', linestyle='dotted')
        ax.axhline(y=0, xmin=0, xmax=1, linewidth=2, color='black')
        # plt.gca().invert_yaxis()
        ax.legend()
        plt.title("Flight envelope")
        plt.show()

flight_envelope = Flight_Envelope(19000, 1.7, 320, 500, 1.225, 2200)
flight_envelope.manoeuver_loads()
flight_envelope.plot()
# print((flight_envelope.manoeuver_array, flight_envelope.neg_manoeuver_array))





