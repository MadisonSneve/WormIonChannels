from IonChannel import IonChannel
from IonChannel import mV, ms
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

#I_x = g_x*(V-V_x) = g_x * m_p_x * h_q_x * (V-V_x)
class SHK1(IonChannel):
    def __init__(self):
        super()
        self.nr_variables = 2

    def compute_derivatives(self, y, t, voltage):
        def m_SHK1_inf(V):
            # updated version
            m_v_05 = 20.4 * mV
            k_a = 7.7 * mV
            return (1)/(1+np.exp(-(V-m_v_05)/k_a))

        def tau_m_SHK1(V):
            a = 26.6 * ms
            b = -33.7 * mV
            c = 15.8 * mV
            d = -33.7 * mV
            e = 15.4 * mV
            f = 2.0 * ms
            return a/(np.exp(-(V-b)/c) + np.exp((V-d)/e)) + f

        def h_SHK1_inf(V):
            h_v_05 = -7.0 * mV
            k_i = 5.8 * mV
            return (1)/(1+np.exp((V-h_v_05)/k_i))

        def tau_h_SHK1(V):
            a = 539.2 * ms
            return a

        m_SHK1 = y[0]
        h_SHK1 = y[1]

        dm_SHK1_dt = (m_SHK1_inf(voltage)-m_SHK1)/(tau_m_SHK1(voltage))

        dh_SHK1_dt = (h_SHK1_inf(voltage)-h_SHK1)/tau_h_SHK1(voltage)

        #I_SHK1 = g_SHK1 * m_3_SHK1 * (0.7 * h_f_SHK1 + 0.3 * h_s_SHK1) * (V - V_K)

        return np.asarray([dm_SHK1_dt, dh_SHK1_dt])


channel = SHK1()
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

    m_SHK1_pre = result_pre[:,0]
    h_SHK1_pre = result_pre[:,1]

    m_SHK1_post = result_post[:,0]
    h_SHK1_post = result_post[:,1]

    if True:
        g_SHK1 = 1.8 * 1e-9
        V_K = -80*1e-3
        I_SHK1_pre = np.transpose(g_SHK1 * np.power(m_SHK1_pre,3) * (1.0 * h_SHK1_pre) * (V_pre - V_K))
        I_SHK1_post = np.transpose(g_SHK1 * np.power(m_SHK1_post,3) * (1.0 * h_SHK1_post) * (V_post - V_K))

        print(I_SHK1_pre)
        print()
        print(I_SHK1_post)
        plt.plot(np.linspace(t0,t2,tsteps1+tsteps2),np.hstack((I_SHK1_pre, I_SHK1_post)))
        plt.xlabel("Time (ms)")
        plt.ylabel("I_SHK1 (arbitrary units)")

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
