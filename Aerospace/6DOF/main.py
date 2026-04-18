import math
import ussa1976
import matplotlib.pyplot as plt
import numpy as np

from governing_equations import flat_earth_eom
from numerical_integrators import numerical_integrator
from tools.Interpolators import interp
from vehicle_models.sphere import spheres

from tools.profutils import profile

@profile(output_dir='logs')
def main():
    # ##########################################################
    # Part1: Initialization
    # ##########################################################

    # Atmospheric data
    atmosphere = ussa1976.compute()

    # Essential gravity and atmospheric data
    amod_alt_m       = atmosphere['z'].values
    amod_rho_kgpm3   = atmosphere['rho'].values
    amod_c_mps       = atmosphere['cs'].values
    amod_g_mps2      = ussa1976.core.compute_gravity(amod_alt_m)

    # amod: atmospheric model with gravity
    amod = {
        "alt_m"     : amod_alt_m,
        "rho_kgpm3" : amod_rho_kgpm3,
        "c_mps"     : amod_c_mps,
        "g_mps2"    : amod_g_mps2
    }

    # vmod: vehicle model
    vmod = spheres.BowlingBall()

    # Setting init conditions or trim conditions
    u0_b_mps   = 0.001  
    v0_b_mps   = 0
    w0_b_mps   = 0
    p0_b_rps   = 0
    q0_b_rps   = 0
    r0_b_rps   = 0
    phi0_rad   = 0*math.pi/180
    theta0_rad = -90*math.pi/180
    psi0_rad   = 0
    p10_n_m    = 0
    p20_n_m    = 0
    p30_n_m    = -30000

    # Assignin' init condtions to the init state array
    x0 = np.array([
        u0_b_mps,
        v0_b_mps,
        w0_b_mps,
        p0_b_rps,
        q0_b_rps,
        r0_b_rps,
        phi0_rad,
        theta0_rad,
        psi0_rad,
        p10_n_m,
        p20_n_m,
        p30_n_m,
    ])

    x0 = x0.transpose()
    nx0 = x0.size

    # Set time conditions
    t0_s = 0.0
    tf_s = 10.0
    h_s  = 0.01

    # ##########################################################
    # Part2: Numerically approx the governing eqs
    # ##########################################################

    # Preallocating soln array
    t_s = np.arange(t0_s, tf_s + h_s, h_s)
    nt_s = t_s.size
    x = np.empty((nx0, nt_s), dtype=float)

    # Assign init condition to soln array
    x[:, 0] = x0

    t_s, x = numerical_integrator.forward_euler(flat_earth_eom.flat_earth_eom, t_s, x, h_s, vmod=vmod, amod=amod)

    # Data pre-alloc & post processing

    # True Airspeed calc 
    TAS_mps = np.zeros((nt_s,1))
    for i, _ in enumerate(t_s):
        TAS_mps[i,0] = math.sqrt(x[0,i]**2 + x[1,i]**2 + x[2,i]**2)

    # Altitude, speed of sound, air density
    Altitude_m  = np.zeros((nt_s, 1))
    Cs_mps      = np.zeros((nt_s, 1))
    Rho_kgm3    = np.zeros((nt_s, 1))

    for i, _ in enumerate(t_s):
        Altitude_m[i, 0] = -x[11,i]
        Cs_mps[i, 0]     = interp(amod['alt_m'], amod['c_mps'], Altitude_m[i,0])
        Rho_kgm3[i, 0]   = interp(amod['alt_m'], amod['rho_kgpm3'], Altitude_m[i,0])

    # Angle of attack
    Alpha_rad = np.zeros((nt_s, 1))
    for i, _ in enumerate(t_s):
        if x[2,i] == 0 and x[0,i] == 0:
            w_over_u = 0
        else:
            w_over_u = x[2,i]/x[0,i]
        
        Alpha_rad[i,0] = math.atan(w_over_u)

    # Angle of sideslip
    Beta_rad = np.zeros((nt_s, 1))
    for i, _ in enumerate(t_s):
        if x[1, i] == 0 and TAS_mps[i,0] == 0:
            v_over_VT = 0
        else:
            v_over_VT = x[1, i]/TAS_mps[i,0]
        
        Beta_rad[i, 0] = math.asin(v_over_VT)

    # Mach number
    Mach = np.zeros((nt_s, 1))
    for i, _ in enumerate(t_s):
        Mach[i,0] = TAS_mps[i,0]/Cs_mps[i,0]

    # ##########################################################
    # Part3: Plot Data
    # ##########################################################

    # Subplots
    fig, axes = plt.subplots(2, 4, figsize=(10,6))
    fig.set_facecolor('black')

    # x-axis velocity
    axes[0,0].plot(t_s, x[0,:], label='u', color='yellow')
    axes[0,0].set_xlabel('Time (sec)', color='white')
    axes[0,0].set_ylabel('u [m/s]', color='white')
    axes[0,0].grid(True)
    axes[0,0].set_facecolor('black')
    axes[0,0].tick_params(colors='white')

    # y-axis velocity
    axes[0,1].plot(t_s, x[1,:], label='v', color='yellow')
    axes[0,1].set_xlabel('Time (sec)', color='white')
    axes[0,1].set_ylabel('v [m/s]', color='white')
    axes[0,1].grid(True)
    axes[0,1].set_facecolor('black')
    axes[0,1].tick_params(colors='white')

    # z-axis velocity
    axes[0,2].plot(t_s, x[2,:], label='w', color='yellow')
    axes[0,2].set_xlabel('Time (sec)', color='white')
    axes[0,2].set_ylabel('w [m/s]', color='white')
    axes[0,2].grid(True)
    axes[0,2].set_facecolor('black')
    axes[0,2].tick_params(colors='white')

    # phi
    axes[0,3].plot(t_s, x[6,:], label='phi', color='yellow')
    axes[0,3].set_xlabel('Time (sec)', color='white')
    axes[0,3].set_ylabel('phi [rad]', color='white')
    axes[0,3].grid(True)
    axes[0,3].set_facecolor('black')
    axes[0,3].tick_params(colors='white')

    # roll rate, p
    axes[1,0].plot(t_s, x[3,:], label='p', color='yellow')
    axes[1,0].set_xlabel('Time (sec)', color='white')
    axes[1,0].set_ylabel('p [r/s]', color='white')
    axes[1,0].grid(True)
    axes[1,0].set_facecolor('black')
    axes[1,0].tick_params(colors='white')

    # pitch rate, q
    axes[1,1].plot(t_s, x[4,:], label='q', color='yellow')
    axes[1,1].set_xlabel('Time (sec)', color='white')
    axes[1,1].set_ylabel('q [r/s]', color='white')
    axes[1,1].grid(True)
    axes[1,1].set_facecolor('black')
    axes[1,1].tick_params(colors='white')

    # yaw rate, r
    axes[1,2].plot(t_s, x[5,:], label='r', color='yellow')
    axes[1,2].set_xlabel('Time (sec)', color='white')
    axes[1,2].set_ylabel('r [r/s]', color='white')
    axes[1,2].grid(True)
    axes[1,2].set_facecolor('black')
    axes[1,2].tick_params(colors='white')

    # pitch, theta
    axes[1,3].plot(t_s, x[7,:], label='theta', color='yellow')
    axes[1,3].set_xlabel('Time (sec)', color='white')
    axes[1,3].set_ylabel('theta [m/s]', color='white')
    axes[1,3].grid(True)
    axes[1,3].set_facecolor('black')
    axes[1,3].tick_params(colors='white')

    plt.tight_layout()
    plt.savefig('savefig/sphere_drop_test_1.png')
    plt.show()

if __name__ == "__main__":
    main()