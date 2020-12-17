from IonChannel import IonChannel
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

#I_x = g_x*(V-V_x) = g_x * m_p_x * h_q_x * (V-V_x)
class SHL1(IonChannel):
    def compute_derivatives(y, t, voltage):
        def m_SHL1_inf(V):
            # updated version
            m_v_05 = 11.2 * mV #potentially -6.8
            k_a = 14.1 * mV
            return (1)/(1+np.exp(-(V-m_v_05)/k_a))

        def tau_m_SHL1(V):
            a = 13.8 * ms #potentially 1.4
            a = 1.4 * ms
            b = -17.5 * mV
            c = 12.9 * mV
            d = -3.7 * mV
            e = 6.5 * mV
            f = 1.9 * ms #potentially 0.2
            f = 0.2 * ms
            return a/(np.exp(-(V-b)/c) + np.exp((V-d)/e)) + f

        def h_f_SHL1_inf(V):
            h_v_05 = -33.1 * mV
            k_i = 8.3 * mV
            return (1)/(1+np.exp((V-h_v_05)/k_i))

        def h_s_SHL1_inf(V):
            h_v_05 = -33.1 * mV
            k_i = 8.3 * mV
            return (1)/(1+np.exp((V-h_v_05)/k_i))

        def tau_h_f_SHL1(V):
            a = 539.2 * ms #potentially 53.9
            b = -28.2 * mV
            c = 4.9 * mV
            d = 27.3 * ms #potentially 2.7
            return ((a)/(1 + np.exp((V-b)/c))) + d

        def tau_h_s_SHL1(V):
            a = 8422.0 * ms #potentially 842.2
            b = -37.7 * mV
            c = 6.4 * mV
            d = 118.9 * ms #potentially 11.9
            return ((a)/(1 + np.exp((V-b)/c))) + d

        #I_SHL1 = g_SHL1 * (m_SHL1)**3 * (0.7 * h_f_SHL1 + 0.3 * h_s_SHL1) * (V - V_K)

        m_SHL1 = y[0]
        h_f_SHL1 = y[1]
        h_s_SHL1 = y[2]

        dm_SHL1_dt = (m_SHL1_inf(voltage)-m_SHL1)/(tau_m_SHL1(voltage))

        dh_f_SHL1_dt = (h_f_SHL1_inf(voltage)-h_f_SHL1)/tau_h_f_SHL1(voltage)
        dh_s_SHL1_dt = (h_s_SHL1_inf(voltage)-h_s_SHL1)/tau_h_s_SHL1(voltage)

        #I_SHL1 = g_SHL1 * m_3_SHL1 * (0.7 * h_f_SHL1 + 0.3 * h_s_SHL1) * (V - V_K)

        return np.asarray([dm_SHL1_dt, dh_f_SHL1_dt, dh_s_SHL1_dt])


init_values = [0.00154979, 0.99648898, 0.99648898]

def run_experiment(V_pre,V_post):
    t0 = 0e-3
    t1 = 100e-3
    t2 = 700e-3

    tsteps1 = int((t1-t0)*1e6)
    tsteps2 = int((t2-t1)*1e6)

    result_pre = integrate.odeint(SHL1.compute_derivatives, np.asarray(init_values), np.linspace(t0,t1,tsteps1), args=(V_pre,))
    print("Steady state:")
    print(result_pre[-1,:])
    result_post = integrate.odeint(SHL1.compute_derivatives, result_pre[-1,:], np.linspace(t1,t2,tsteps2),args=(V_post,))

    m_SHL1_pre = result_pre[:,0]
    h_f_SHL1_pre = result_pre[:,1]
    h_s_SHL1_pre = result_pre[:,2]

    m_SHL1_post = result_post[:,0]
    h_f_SHL1_post = result_post[:,1]
    h_s_SHL1_post = result_post[:,2]

    if True:
        g_SHL1 = 1.8 * 1e-9
        V_K = -80*1e-3
        I_SHL1_pre = np.transpose(g_SHL1 * np.power(m_SHL1_pre,3) * (0.7 * h_f_SHL1_pre + 0.3 * h_s_SHL1_pre) * (V_pre - V_K))
        I_SHL1_post = np.transpose(g_SHL1 * np.power(m_SHL1_post,3) * (0.7 * h_f_SHL1_post + 0.3 * h_s_SHL1_post) * (V_post - V_K))

        print(I_SHL1_pre)
        print()
        print(I_SHL1_post)
        plt.plot(np.linspace(t0,t2,tsteps1+tsteps2),np.hstack((I_SHL1_pre, I_SHL1_post)))
        plt.xlabel("Time (ms)")
        plt.ylabel("I_SHL1 (arbitrary units)")

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
