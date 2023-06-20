import numpy as np
import matplotlib.pyplot as plt

class Flight_Envelope:
    def __init__(self, MTOW, CL_max, CL_alpha, v_cruise, v_dive, rho, wing_loading):
        self.MTOW = MTOW
        self.CL_max = CL_max
        self.CL_alpha = CL_alpha
        self.v_cruise = v_cruise
        self.v_dive = v_dive
        self.rho = rho
        self.wing_loading = wing_loading

        # Set array of velocities
        self.v_array = np.arange(0, v_dive + 0.5, 0.5)

    def manoeuver_loads(self):
        # Positive part
        # Define max manoeuver load factor from STANAG 4671
        self.max_load_factor = min(2.1 + 10900/(self.MTOW+4536), 3.8)
        # Calculate load factors up until max manoeuver load
        n = (0.5 * self.rho * np.square(self.v_array) * self.CL_max) / self.wing_loading
        self.manoeuver_array = np.clip(n, 0, self.max_load_factor)
        self.v_array = np.append(self.v_array, self.v_dive)
        self.manoeuver_array = np.append(self.manoeuver_array, 0)

        # Negative part
        # Define max negative manoeuver load factor from STANAG 4671
        self.min_load_factor = -0.4 * self.max_load_factor
        self.v_array_negative = self.v_array[self.v_array <= self.v_cruise]
        # Calculate load factors up until max negative manoeuver load
        n_neg = -(0.5 * self.rho * np.square(self.v_array_negative) * self.CL_max) / self.wing_loading
        self.neg_manoeuver_array = np.clip(n_neg, self.min_load_factor, 0)
        self.v_array_negative = np.append(self.v_array_negative, self.v_dive)
        self.neg_manoeuver_array = np.append(self.neg_manoeuver_array, 0)

    def gust_loads(self):
        # Gust speeds from STANAG 4671
        v_gust_cruise   = 15.2  # m/s
        v_gust_dive     = 7.6   # m/s
        # Find load factor differential
        def n_diff(v_gust, V):
            dn = (self.rho * V * self.CL_alpha * v_gust) / (2 * self.wing_loading)
            return dn
        self.vc_vd_array = np.array([0,self.v_cruise, self.v_dive, self.v_dive])
        # v_gust_array = [v_gust_cruise, -v_gust_cruise, v_gust_dive, -v_gust_dive]
        # for v_gust in v_gust_array:
        self.gust_array_up = np.array([1, 1+n_diff(v_gust_cruise, self.v_cruise), 1+n_diff(v_gust_dive, self.v_dive), 1])
        self.gust_array_down = np.array([1, 1+n_diff(-v_gust_cruise, self.v_cruise), 1+n_diff(-v_gust_dive, self.v_dive), 1])

    def plot(self):
        fig, ax = plt.subplots()

        ax.plot(self.v_array, self.manoeuver_array, linestyle='-', label='Manoeuver loading', color='blue')
        ax.plot(self.v_array_negative, self.neg_manoeuver_array, linestyle='-', color='blue')
        ax.plot(self.vc_vd_array, self.gust_array_up, linestyle='-', label='Gust loading', color='red')
        ax.plot(self.vc_vd_array, self.gust_array_down, linestyle='-', color='red')


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

flight_envelope = Flight_Envelope(6521.24, 1.4, 3.724, 110, 1.4*110, 0.770816, 1412.107)
flight_envelope.manoeuver_loads()
flight_envelope.gust_loads()
flight_envelope.plot()


