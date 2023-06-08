import numpy as np

class LandingGear:
    def __init__(self, xcg, span, sweep_QC, MTOW, chord_tip, chord_root, MAC, dihedral):
        self.xcg = xcg
        self.span = span
        self.sweep = sweep_QC
        self.MTOW = MTOW
        self.chord_tip = chord_tip
        self.chord_root = chord_root
        self.MAC = MAC
        self.dihedral = dihedral

    def positioning(self):
        # Assumed initial MLG position
        self.pos_x_MLG = self.xcg + 0.15 * self.MAC
        # NLG distance flowing from this
        load_percentage_NLG = 0.1   # Set percentage of MTOW NLG supports
        self.pos_x_NLG = (self.xcg - self.pos_x_MLG*(1-load_percentage_NLG))/load_percentage_NLG

        # Find wingtip distance
        self.dis_x_cg_wingtip = 0.25 * self.chord_root + np.tan(self.sweep) * self.span/2 + 0.75 * self.chord_tip - self.xcg
        # Find lm and ln
        self.dis_x_cg_NLG = self.xcg - self.pos_x_NLG
        self.dis_x_cg_MLG = self.pos_x_MLG - self.xcg

    def loads(self):
        # Set number of struts
        n_s = 1
        # Loads on the landing gears
        self.load_NLG = (self.MTOW * self.dis_x_cg_MLG) / (self.dis_x_cg_NLG + self.dis_x_cg_MLG)
        self.load_MLG = (self.MTOW * self.dis_x_cg_NLG) / (n_s * (self.dis_x_cg_NLG + self.dis_x_cg_MLG))

    def clearance(self, pusher = False):
        # Longitudinal clearance
        self.ang_lon = np.radians(15)   # degrees

        # Assess the clearance in case of TE pusher engines, most outer engines then critical
        pos_y_engine = self.span / 2
        pos_x_engine = self.xcg + self.dis_x_cg_wingtip
        prop_radius = 0.5
        if pusher:
            # Option if wingtip is constraining
            self.height_MLG_tip = np.tan(self.ang_lon) * (self.dis_x_cg_wingtip - self.dis_x_cg_MLG) - np.tan(self.dihedral) * self.span/2
            # Option if pusher engine is constraining
            self.height_MLG_pusher = np.tan(self.ang_lon) * (pos_x_engine - self.pos_x_MLG) + prop_radius - np.tan(self.dihedral) * pos_y_engine
            self.height_MLG = max(self.height_MLG_pusher,self.height_MLG_tip)
        else:
            self.height_MLG = np.tan(self.ang_lon) * (self.dis_x_cg_wingtip - self.dis_x_cg_MLG) - np.tan(self.dihedral) * self.span/2

        # Lateral tip-over
        self.ang_tip = np.radians(55)  # degrees
        # Perpendicular distance
        w = self.height_MLG / np.tan(self.ang_tip)
        # Angle from NLG to half track width
        beta = np.arcsin(w / (self.xcg - self.pos_x_NLG))
        self.track_width_MLG = 2 * np.tan(beta) * (self.pos_x_MLG - self.pos_x_NLG)

        # Lateral clearance
        if pusher:
            self.ang_lat = np.arctan((self.height_MLG + np.tan(self.dihedral) * pos_y_engine - prop_radius) / (pos_y_engine - self.track_width_MLG / 2))
        else:
            self.ang_lat = np.arctan((self.height_MLG + np.tan(self.dihedral) * self.span/2) / (self.span / 2 - self.track_width_MLG / 2))

        if self.ang_lat <= np.radians(5):
            print("Lateral clearance criterion not fulfilled!")



drone = LandingGear(4.342, 16.4869908, np.radians(38), 6521.24, 1.550632326, 5.800, 3, 0)
drone.positioning()
drone.clearance()
drone.loads()