from IonChannel import IonChannel
from IonChannel import mV, ms
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

#I_x = g_x*(V-V_x) = g_x * m_p_x * h_q_x * (V-V_x)
class KVS1(IonChannel):
    def __init__(self):
        super()
        self.nr_variables = 2

    def compute_derivatives(self, y, t, voltage):
        def m_KVS1_inf(V):
            # updated version
            m_v_05 = 57.1 * mV
            k_a = 25.0 * mV
            return (1)/(1+np.exp(-(V-m_v_05)/k_a))

        def tau_m_KVS1(V):
            a = 30.0 * ms
            b = 18.1 * mV
            c = 20 * mV
            d = 1.0 * ms
            return a/(np.exp(-(V-b)/c)) + d

        def h_KVS1_inf(V):
            h_v_05 = 47.3 * mV
            k_i = 11.1 * mV
            return (1)/(1+np.exp((V-h_v_05)/k_i))

        def tau_h_KVS1(V):
            a = 88.5 * ms
            b = 50.0 * mV
            c = 15.0 * mV
            d = 53.4 * ms
            return a/(np.exp(-(V-b)/c)) + d

        m_KVS1 = y[0]
        h_KVS1 = y[1]

        dm_KVS1_dt = (m_KVS1_inf(voltage)-m_KVS1)/(tau_m_KVS1(voltage))

        dh_KVS1_dt = (h_KVS1_inf(voltage)-h_KVS1)/tau_h_KVS1(voltage)

        #I_KVS1 = g_KVS1 * m_3_KVS1 * (0.7 * h_f_KVS1 + 0.3 * h_s_KVS1) * (V - V_K)

        return np.asarray([dm_KVS1_dt, dh_KVS1_dt])


channel = KVS1()
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

    m_KVS1_pre = result_pre[:,0]
    h_KVS1_pre = result_pre[:,1]

    m_KVS1_post = result_post[:,0]
    h_KVS1_post = result_post[:,1]

    if True:
        g_KVS1 = 1.8 * 1e-9
        V_K = -80*1e-3
        I_KVS1_pre = np.transpose(g_KVS1 * np.power(m_KVS1_pre,3) * (1.0 * h_KVS1_pre) * (V_pre - V_K))
        I_KVS1_post = np.transpose(g_KVS1 * np.power(m_KVS1_post,3) * (1.0 * h_KVS1_post) * (V_post - V_K))

        print(I_KVS1_pre)
        print()
        print(I_KVS1_post)
        plt.plot(np.linspace(t0,t2,tsteps1+tsteps2),np.hstack((I_KVS1_pre, I_KVS1_post)))
        plt.xlabel("Time (ms)")
        plt.ylabel("I_KVS1 (arbitrary units)")

run_experiment(-80e-3, 40e-3)
run_experiment(-80e-3, 30e-3)
run_experiment(-80e-3, 20e-3)
run_experiment(-80e-3, 10e-3)
run_experiment(-80e-3, 0e-3)
run_experiment(-80e-3, -10e-3)
run_experiment(-80e-3, -20e-3)
run_experiment(-80e-3, -30e-3)
run_experiment(-80e-3, -40e-3)

plt.savefig('KVS1.png')
plt.show()
