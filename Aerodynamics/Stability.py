import numpy as np
class Stab:
    def __init__(self,plane,Coefficients,MTOW):
        self.plane=plane
        self.Coeff=Coefficients
        self.V0=plane.V

        self.muc=MTOW/(Coefficients.rho*plane.S*plane.MAC)
        self.mub=MTOW/(Coefficients.rho*plane.S*plane.b_tot)

    def get_asymm_eigen(self):
        
        self.A_asymm = 16*(self.mub**3)*(self.Coeff.KX2*self.Coeff.KZ2-self.Coeff.KXZ**2)  #squared KXZ

        self.B_asymm = -4*(self.mub**2)*(2*self.Coeff.C_Y_b*(self.Coeff.KX2*self.Coeff.KZ2-self.Coeff.KXZ**2)+self.Coeff.C_n_r*self.Coeff.KX2+self.Coeff.C_l_p*self.Coeff.KZ2+(
            self.Coeff.C_l_r+self.Coeff.C_n_p)*self.Coeff.KXZ)  # squared KXZ

        self.C_asymm = 2*self.mub*((self.Coeff.C_Y_b*self.Coeff.C_n_r-self.Coeff.C_Y_r*self.Coeff.C_n_beta)*self.Coeff.KX2+(
            self.Coeff.C_Y_b*self.Coeff.C_l_p-self.Coeff.C_l_beta*self.Coeff.C_Y_p)*self.Coeff.KZ2+((
            self.Coeff.C_Y_b*self.Coeff.C_n_p-self.Coeff.C_n_beta*self.Coeff.C_Y_p)+(
            self.Coeff.C_Y_b*self.Coeff.C_l_r-self.Coeff.C_l_beta*self.Coeff.C_Y_r))*self.Coeff.KXZ+4*self.mub*self.Coeff.C_n_beta*self.Coeff.KX2+4*self.mub*self.Coeff.C_l_beta*self.Coeff.KXZ+0.5*(self.Coeff.C_l_p*self.Coeff.C_n_r-self.Coeff.C_n_p*self.Coeff.C_l_r))

        self.D_asymm = -4*self.mub*self.Coeff.C_L*(self.Coeff.C_l_beta*self.Coeff.KZ2+self.Coeff.C_n_beta*self.Coeff.KXZ)+2*self.mub*(
            self.Coeff.C_l_beta*self.Coeff.C_n_p-self.Coeff.C_n_beta*self.Coeff.C_l_p)+0.5*self.Coeff.C_Y_b*(
            self.Coeff.C_l_r*self.Coeff.C_n_p-self.Coeff.C_n_r*self.Coeff.C_l_p)+0.5*self.Coeff.C_Y_p*(
            self.Coeff.C_l_beta*self.Coeff.C_n_r-self.Coeff.C_n_beta*self.Coeff.C_l_r)+0.5*self.Coeff.C_Y_r*(
            self.Coeff.C_l_p*self.Coeff.C_n_beta-self.Coeff.C_n_p*self.Coeff.C_l_beta) # Added minus sign to the front

        self.E_asymm = self.Coeff.C_Y_b*(self.Coeff.C_l_beta*self.Coeff.C_n_r-self.Coeff.C_n_beta*self.Coeff.C_l_r)

        self.eigenvalues_asymm=np.roots([self.A_asymm, self.B_asymm, self.C_asymm, self.D_asymm, self.E_asymm])

    def get_symm_eigen(self):
        self.A_symm = 4 * (self.muc ** 2) * self.Coeff.KY2 * (self.Coeff.C_Z_alphadot - 2 * self.muc)

        self.B_symm = self.Coeff.C_m_alphadot * 2 * self.muc * (self.Coeff.C_Z_q + 2 * self.muc) - self.Coeff.C_m_q * 2 * self.muc * (
                    self.Coeff.C_Z_alphadot - 2 * self.muc) - 2 * self.muc * self.Coeff.KY2 * (
                             self.Coeff.C_X_u * (self.Coeff.C_Z_alphadot - 2 * self.muc) - 2 * self.muc * self.Coeff.C_Z_alpha)

        self.C_symm = self.Coeff.C_m_alpha * 2 * self.muc * (self.Coeff.C_Z_q + 2 * self.muc) - self.Coeff.C_m_alphadot * (
                    2 * self.muc * self.Coeff.C_X_0 + self.Coeff.C_X_u * (self.Coeff.C_Z_q + 2 * self.muc)) + self.Coeff.C_m_q * (
                             self.Coeff.C_X_u * (
                                 self.Coeff.C_Z_alphadot - 2 * self.muc) - 2 * self.muc * self.Coeff.C_Z_alpha) + 2 * self.muc * (
                     self.Coeff.KY2) * (self.Coeff.C_X_alpha * self.Coeff.C_Z_u - self.Coeff.C_Z_alpha * self.Coeff.C_X_u)

        self.D_symm = self.Coeff.C_m_u * (self.Coeff.C_X_alpha * (self.Coeff.C_Z_q + 2 * self.muc) - self.Coeff.C_Z_0 * (
                    self.Coeff.C_Z_alphadot - 2 * self.muc)) - self.Coeff.C_m_alpha * (2 * self.muc * self.Coeff.C_X_0 + self.Coeff.C_X_u * (
                    self.Coeff.C_Z_q + 2 * self.muc)) + self.Coeff.C_m_alphadot * (
                             self.Coeff.C_X_0 * self.Coeff.C_X_u - self.Coeff.C_Z_0 * self.Coeff.C_Z_u) + self.Coeff.C_m_q * (
                             self.Coeff.C_X_u * self.Coeff.C_Z_alpha - self.Coeff.C_Z_u * self.Coeff.C_X_alpha)

        self.E_symm = -self.Coeff.C_m_u * (self.Coeff.C_X_0 * self.Coeff.C_X_alpha + self.Coeff.C_Z_0 * self.Coeff.C_Z_alpha) + self.Coeff.C_m_alpha * (
                    self.Coeff.C_X_0 * self.Coeff.C_X_u + self.Coeff.C_Z_0 * self.Coeff.C_Z_u)
    
        self.eigenvalues_symm=np.roots([self.A_symm, self.B_symm, self.C_symm, self.D_symm, self.E_symm])

    def Routh_discriminant(self):
        self.Routh=self.B_asymm*self.C_asymm*self.D_asymm-self.D_asymm**2*self.A_asymm-self.B_asymm**2*self.E_asymm
        return self.Routh
    
    def add_text_to_file(self,file_path, text):
        with open(file_path, 'a') as file:
            file.write(text)

    def record_stability(self):
        text='Symmetric eigenvalues: '+str(self.eigenvalues_symm)+'\nAsymmetric eigenvalues: '+str(self.eigenvalues_asymm)+'\nCoefficients'+str(self.Coeff.coefficients)+'\n\n'
        print(text)
        self.add_text_to_file('./Record/Stability record.txt', text)
        