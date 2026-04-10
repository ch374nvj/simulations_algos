import math

import matplotlib.pyplot as plt
import numpy as np

from governing_equations import flat_earth_eom
from numerical_integrators import numerical_integrator

# ##########################################################
# Part1: Initialization
# ##########################################################

# Vehicle Defn
r_sphere_m = 0.08
m_sphere_kg = 5
J_sphere_kgm2 = 0.4*m_sphere_kg*r_sphere_m

amod = {
    "m_kg" : 1,
    "Jxz_b_kgm2" : 0,
    "Jxx_b_kgm2" : J_sphere_kgm2,
    "Jyy_b_kgm2" : J_sphere_kgm2,
    "Jzz_b_kgm2" : J_sphere_kgm2
}

# Setting init conditions or trim conditions
u0_b_mps   = 0  
v0_b_mps   = 0
w0_b_mps   = 0
p0_b_rps   = 0
q0_b_rps   = 0
r0_b_rps   = 0
phi0_rad   = 90*math.pi/180
theta0_rad = 0*math.pi/180
psi0_rad   = 0
p10_n_m    = 0
p20_n_m    = 0
p30_n_m    = 0

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

t_s, x = numerical_integrator.forward_euler(flat_earth_eom.flat_earth_eom, t_s, x, h_s, amod)

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