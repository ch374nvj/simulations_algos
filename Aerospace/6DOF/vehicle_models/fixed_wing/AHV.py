"""
This library returns an object which contains the following intertial 
properties of the called object.
m_kg        = vmod['m_kg']
Jxz_b_kgm2  = vmod['Jxz_b_kgm2']
Jxx_b_kgm2  = vmod['Jxx_b_kgm2']
Jyy_b_kgm2  = vmod['Jyy_b_kgm2']
Jzz_b_kgm2  = vmod['Jzz_b_kgm2']
CD_approx = vmod['CD_approx']
CL_approx = vmod['CL_approx']
CY_approx = vmod['CY_approx']
Aref_m2 = vmod['Aref_m2]

References: 
    J.D. Shaughnessy, S.Z. Pinckney, J.D. McMinn, C.I. Cruz, M.-L. Kelley, Hypersonic vehicle
        simulation model: winged-cone configuration (1990)
    K. Sachan, R. Padhi, State-constrained robust adaptive cruise control design for air-breathing
        hypersonic vehicles, in 2018 AIAA Guidance, Navigation, and Control Conference (2018)
    K. Sachan, R. Padhi, Nonlinear robust neuro-adaptive flight control for hypersonic 
        vehicles with state constraints, Control Engineering Practice (2020)
"""

class AHV_model:
    def __init__(self) -> None:
        self.m_kg = 136080
        self.Aref_m2 = 334.73
        self.b_m = 18.29
        self.c_m = 24.38
        self.T_max_N = 1467900
        self.Jxx_b_kgm2 = -7.8809e-5 * self.m_kg**2 + 25.88576 * self.m_kg - 6.9683e+5
        self.Jyy_b_kgm2 = -8.2890e-4 * self.m_kg**2 + 265.9889 * self.m_kg - 7.3048e+6
        self.Jzz_b_kgm2 = self.Jyy_b_kgm2
        self.Jxz_b_kgm2 = 0
        self.M = 1 # Mach number

        self.vmod = {
                    "m_kg" : self.m_kg,
                    "Jxx_b_kgm2" : self.Jxx_b_kgm2,
                    "Jyy_b_kgm2" : self.Jyy_b_kgm2,
                    "Jzz_b_kgm2" : self.Jzz_b_kgm2,
                    "Jxz_b_kgm2" : self.Jxz_b_kgm2,
                    # "Vterm_mps"  : self.Vterm_mps,
                    # "CD_approx"  : self.CD_approx,
                    # "CL_approx"  : 0.0,
                    # "CY_approx"  : 0.0,
                    "Aref_m2"    : self.Aref_m2   
                }

    def model_update(self, TAS, c):
        self.M = TAS/c

    def CD(self, alpha_deg, u_deg) ->  float:
        Cd_alpha = 8.717e-2 - 3.307e-2 * self.M + 3.179e-2 * alpha_deg + 5.036e-3 * self.M ** 2
        Cd_del_e = 0
        Cd_del_a = 0
        Cd_del_r = 0

        Cd = Cd_alpha + Cd_del_e + Cd_del_a + Cd_del_r
        return Cd

    def CL(self, alpha_deg, u_deg) -> float:
        Cl_alpha = -8.19e-2 + 4.7e-2 * self.M + 1.86e-2 * alpha_deg - 9.19e-3 * self.M ** 2
        Cl_del_e = -4.14e-4 * u_deg[0]
        Cl_del_a = -4.14e-4 * u_deg[1]
        Cl_del_r = 0

        Cl = Cl_alpha + Cl_del_e + Cl_del_a + Cl_del_r
        return Cl

    def CY(self, beta_deg, u_deg):
        Cy_beta = -0.292*self.M + 5.48e-2 * self.M**2 - 4.32e-3 * self.M**3
        Cy_del_e = 0
        Cy_del_a = 0
        Cy_del_r = 3.84e-4 * u_deg[2]

        Cy = Cy_beta * beta_deg + Cy_del_e + Cy_del_a + Cy_del_r
        return Cy

    def Clm(self, alpha_deg, beta_deg, u_deg, TAS, p, r):
        p *= 57.3
        r *= 57.3
        Cl_beta = -0.14 + 3.32e-2 * self.M - 7.59e-4 * alpha_deg - 3.79e-3 * self.M ** 2
        Cl_del_e = -1.17e-4 * u_deg[0]
        Cl_del_a = 1.17e-4  * u_deg[1]
        Cl_del_r = 1.144e-4 * u_deg[2]

        Cl_p = -0.299 + 7.47e-2 * self.M + 1.38e-3 * alpha_deg - 9.13e-3 * self.M**2
        Cl_r = 0.382 - 0.106*self.M + 1.94e-3 * alpha_deg + 1.45e-2 * self.M**2 - 1.02e-3 * self.M**3
        
        Clm = Cl_beta * beta_deg + Cl_del_e + Cl_del_a + Cl_del_r +  Cl_p * (p * self.b_m / (2*TAS)) + Cl_r * (r * self.b_m / (2*TAS))
        return Clm

    def Cm(self, alpha_deg, u_deg, TAS, q):
        q*= 57.3
        Cm_alpha = -2.19e-2 + 7.73e-3 * self.M - 2.26e-3 * alpha_deg
        Cm_del_e = 2.89e-4 * u_deg[0]
        Cm_del_a = 2.89e-4 * u_deg[1]
        Cm_del_r = 1.58e-3 * self.M**2

        Cm_q = -1.36 + 0.38 * self.M - 5.42e-2 * self.M**2 + 3.8e-3 * self.M ** 3

        
        Cm = Cm_alpha + Cm_del_e + Cm_del_a + Cm_del_r + Cm_q * (q * self.c_m / (2*TAS))
        return Cm

    def Cn(self, alpha_deg, beta_deg, u_deg, TAS, p, r):
        p *= 57.3
        r *= 57.3
        Cn_beta = 5.91e-2 * self.M - 1.48e-2 * self.M**2 + 1.27e-3 * self.M ** 3
        Cn_del_e = 0
        Cn_del_a = 0
        Cn_del_r = -5.28e-4 * u_deg[2]

        Cn_p = 0.368 - 9.79e-2 * self.M + 1.24e-2 * self.M**2
        Cn_r = -2.41 + 0.596 * self.M - 2.74e-3 * alpha_deg - 7.57e-2 * self.M**2 + 4.9e-3 * self.M**3

        Cn = Cn_beta * beta_deg + Cn_del_e + Cn_del_a + Cn_del_r + Cn_p * (p * self.b_m / (2 * TAS)) + Cn_r * (r * self.b_m / (2 * TAS))
        return Cn
