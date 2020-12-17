from IonChannel import IonChannel
from IonChannel import mV, ms
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

#TODO not complete yet
#I_x = g_x*(V-V_x) = g_x * m_p_x * h_q_x * (V-V_x)
class KQT3(IonChannel):
    def __init__(self):
        super()
        self.nr_variables = 2

    def compute_derivatives(self, y, t, voltage):
        def m_KQT3_inf(V):
            # updated version
            m_v_05 = -12.8 * mV
            k_a = 15.8 * mV
            return (1)/(1+np.exp(-(V-m_v_05)/k_a))

        def tau_m_KQT3(V):
            a = 30.0 * ms
            b = 18.1 * mV
            c = 20 * mV
            d = 1.0 * ms
            return a/(np.exp(-(V-b)/c)) + d

        def h_KQT3_inf(V):
            h_v_05 = 47.3 * mV
            k_i = 11.1 * mV
            return (1)/(1+np.exp((V-h_v_05)/k_i))

        def tau_h_KQT3(V):
            a = 88.5 * ms
            b = 50.0 * mV
            c = 15.0 * mV
            d = 53.4 * ms
            return a/(np.exp(-(V-b)/c)) + d

        m_KQT3 = y[0]
        h_KQT3 = y[1]

        dm_KQT3_dt = (m_KQT3_inf(voltage)-m_KQT3)/(tau_m_KQT3(voltage))

        dh_KQT3_dt = (h_KQT3_inf(voltage)-h_KQT3)/tau_h_KQT3(voltage)

        #I_KQT3 = g_KQT3 * m_3_KQT3 * (0.7 * h_f_KQT3 + 0.3 * h_s_KQT3) * (V - V_K)

        return np.asarray([dm_KQT3_dt, dh_KQT3_dt])


channel = KQT3()
init_values = channel.find_steady_state()

def run_experiment(V_pre,V_post):
    t0 = 0e-3
    t1 = 100e-3
    t2 = 700e-3

    tsteps1 = int((t1-t0)*1e6)
    tsteps2 = int((t2-t1)*1e6)

    result_pre = integrate.odeint(channel.compute_derivatives, np.asarray(init_values), np.linspace(t0,t1,tsteps1), args=(V_pre,))
    print(result_pre[-1,:])
    result_post = integrate.odeint(channel.compute_derivatives, result_pre[-1,:], np.linspace(t1,t2,tsteps2),args=(V_post,))

    m_KQT3_pre = result_pre[:,0]
    h_KQT3_pre = result_pre[:,1]

    m_KQT3_post = result_post[:,0]
    h_KQT3_post = result_post[:,1]

    if True:
        g_KQT3 = 1.8 * 1e-9
        V_K = -80*1e-3
        I_KQT3_pre = np.transpose(g_KQT3 * np.power(m_KQT3_pre,3) * (1.0 * h_KQT3_pre) * (V_pre - V_K))
        I_KQT3_post = np.transpose(g_KQT3 * np.power(m_KQT3_post,3) * (1.0 * h_KQT3_post) * (V_post - V_K))

        print(I_KQT3_pre)
        print()
        print(I_KQT3_post)
        plt.plot(np.linspace(t0,t2,tsteps1+tsteps2),np.hstack((I_KQT3_pre, I_KQT3_post)))
        plt.xlabel("Time (ms)")
        plt.ylabel("I_KQT3 (arbitrary units)")

run_experiment(-80e-3, 40e-3)
run_experiment(-80e-3, 30e-3)
run_experiment(-80e-3, 20e-3)
run_experiment(-80e-3, 10e-3)
run_experiment(-80e-3, 0e-3)
run_experiment(-80e-3, -10e-3)
run_experiment(-80e-3, -20e-3)
run_experiment(-80e-3, -30e-3)
run_experiment(-80e-3, -40e-3)

plt.show()
