import numpy as np

def forward_euler(f, t_s, x, u, h_s, vmod, amod):
    """_summary_

    Forward Euler integrator to approximate solution of differential eqn

    Args:
        f   (func) : function representing RHS of a ODE (x_dot = f(x,t))
        t_s (list) : vector of time points at which numerical solution will be approximated
        x   (list) : numerically approximated solns of f(x,t)
        h_s (float): step size in sec
        vmod(dict) : vehicle (aircraft) model
        amod(dict) : atmosphere model

    Returns:
        t_s (list) : vector of time points at which numerical solution was be approximated
        x   (list) : numerically approximated solns of f(x,t)
    """
    
    # Fwd Euler numerical integration
    dx = np.empty((12, len(t_s)), dtype=float)

    for i in range(1, len(t_s)):
        dx[:, i] = f(t_s[i-1], x[:, i-1], u, amod, vmod)
        x[:, i]  = x[:, i-1] + h_s * dx[:, i] 

    return t_s, x, dx

def rk4(f, t_s, x, u, h_s, vmod, amod):
    dx = np.empty((12, len(t_s)), dtype=float)

    for i in range(1,len(t_s)):
        ti = t_s[i-1]
        yi = x[:,i-1]

        # 4th order RK
        k1 = f(ti, yi, u, amod, vmod)
        k2 = f(ti + h_s/2, yi + k1*h_s/2, u, amod, vmod)
        k3 = f(ti + h_s/2, yi + k2*h_s/2, u, amod, vmod)
        k4 = f(ti + h_s, yi + k3*h_s, u, amod, vmod)

        y = yi + (h_s/6)*(k1 + 2*k2 + 2*k3 + k4)
        x[:, i] = y

        dx[:,i] = k1 #unintegrated state rates

    return t_s, x, dx
