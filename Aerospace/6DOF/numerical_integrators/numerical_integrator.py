import numpy as np

def forward_euler(f, t_s, x, h_s, amod):
    """_summary_

    Forward Euler integrator to approximate solution of differential eqn

    Args:
        f   (func) : function representing RHS of a ODE (x_dot = f(x,t))
        t_s (list) : vector of time points at which numerical solution will be approximated
        x   (list) : numerically approximated solns of f(x,t)
        h_s (float): step size in sec
        amod(dict) : aircraft model

    Returns:
        t_s (list) : vector of time points at which numerical solution was be approximated
        x   (list) : numerically approximated solns of f(x,t)
    """
    
    # Fwd Euler numerical integration
    for i in range(1, len(t_s)):
        x[:, i] = x[:, i-1] + h_s * f(t_s[i-1], x[:, i-1], amod)

    return t_s, x;
