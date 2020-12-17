import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
mV = 1e-3; ms = 1e-3

class IonChannel(object):
    def __init__(self):
        self.nr_variables = -1

    def compute_derivatives(self):
        pass

    def find_steady_state(self):
        result_pre = integrate.odeint(self.compute_derivatives, np.zeros(shape=(self.nr_variables,)), np.linspace(0,10000,100000), args=(-80e-3,))
        print("Steady state:")
        return (result_pre[-1,:])
