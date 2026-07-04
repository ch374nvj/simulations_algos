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
"""

class AHV:
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

    def model_update(self, TAS, c):
        self.M = TAS/c

    def CD(self, alpha_deg, u) ->  float:
        Cd_alpha = 8.717e-2 - 3.307e-2 * self.M + 3.179e-2 * alpha_deg + 5.036e-3 * self.M ** 2
        Cd_del_e = 0
        Cd_del_a = 0
        Cd_del_r = 0

        Cd = Cd_alpha + Cd_del_e + Cd_del_a + Cd_del_r
        return Cd

    def CL(self, alpha_deg, u) -> float:
        Cl_alpha = -8.19e-2 + 4.7e-2 * self.M + 1.86e-2 * alpha_deg - 9.19e-3 * self.M ** 2
        Cl_del_e = -4.14e-4 * u[0]
        Cl_del_a = -4.14e-4 * u[1]
        Cl_del_r = 0

        Cl = Cl_alpha + Cl_del_e + Cl_del_a + Cl_del_r
        return Cl

    def CY(self):
        pass

    def Clm(self):
        pass

    def Cm(self):
        pass

    def Cn(self):
        pass
