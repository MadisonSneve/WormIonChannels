import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

#I_x = g_x*(V-V_x) = g_x * m_p_x * h_q_x * (V-V_x)
def SHL1(y, t, voltage):

    mV = 1e-3; ms = 1e-3
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

    old_eqs = '''# m-related equations
    m_k_a = 14.1 #mV
    m_v_05 = 11.2 #mV # or -6.8 calibrated
    #m_v_05 = -6.8
    m_SHL1_inf = lambda V : 1/(1+np.exp(-(V-m_v_05)/m_k_a))
    m_a = 13.8; m_b = -17.5; m_c = 12.9; m_d = -3.7; m_e = 6.5; m_f = 1.9 #m_f could be 0.2 calibrated
    tau_m_SHL1 = lambda V : ((m_a)/((np.exp(-(V-m_b)/m_c)) + (np.exp(-(V-m_d)/m_e)))) + m_f

    print(tau_m_SHL1(-10*1e-3))

    # h-related equations
    tau_h_a_f = 539.2/10; tau_h_b_f = -28.2; tau_h_c_f = 4.9; tau_h_d_f = 27.3 # tau_f_h_a could be 53.9 calibrated
    tau_h_a_s = 8422.0/10; tau_h_b_s = -37.7; tau_h_c_s = 6.4; tau_h_d_s = 118.9/10 # a could be 842.2; d could be 11.9 calibrated
    h_v_05 = -33.1
    h_k_i  = 8.3
    h_f_SHL1_inf = lambda V : (1)/(1+np.exp((V-h_v_05)/h_k_i))
    h_s_SHL1_inf = lambda V : (1)/(1+np.exp((V-h_v_05)/h_k_i))
    tau_h_f_SHL1 = lambda V : ((tau_h_a_f)/(1 + np.exp((V-tau_h_b_f)/tau_h_c_f))) + tau_h_d_f
    tau_h_s_SHL1 = lambda V : ((tau_h_a_s)/(1 + np.exp((V-tau_h_b_s)/tau_h_c_s))) + tau_h_d_s
    '''

    #print(tau_m_SHL1(-10*1e-3), "expecting 18ms")
    #print(tau_h_f_SHL1(-10*1e-3),"expecting 40ms")
    #print(tau_h_s_SHL1(-10*1e-3),"expecting 200ms")

    dm_SHL1_dt = (m_SHL1_inf(voltage)-m_SHL1)/(tau_m_SHL1(voltage))

    dh_f_SHL1_dt = (h_f_SHL1_inf(voltage)-h_f_SHL1)/tau_h_f_SHL1(voltage)
    dh_s_SHL1_dt = (h_s_SHL1_inf(voltage)-h_s_SHL1)/tau_h_s_SHL1(voltage)

    #vv = np.linspace(-100*1e-3,70*1e-3,1000)
    #hinf = h_f_SHL1_inf(vv)
    #taum = tau_m_SHL1(vv)
    #minf = m_SHL1_inf(vv)
    #plt.plot(vv, hinf)
    #plt.plot(vv, minf)
    #plt.plot(vv, taum)
    #plt.show()


    #I_SHL1 = g_SHL1 * m_3_SHL1 * (0.7 * h_f_SHL1 + 0.3 * h_s_SHL1) * (V - V_K)

    return np.asarray([dm_SHL1_dt, dh_f_SHL1_dt, dh_s_SHL1_dt])


init_values = [0.00154979, 0.99648898, 0.99648898]

def run_experiment(V_pre,V_post):
    t0 = 0e-3
    t1 = 100e-3
    t2 = 700e-3

    tsteps1 = int((t1-t0)*1e6)
    tsteps2 = int((t2-t1)*1e6)

    result_pre = integrate.odeint(SHL1, np.asarray(init_values), np.linspace(t0,t1,tsteps1), args=(V_pre,))
    print("Steady state:")
    print(result_pre[-1,:])
    result_post = integrate.odeint(SHL1, result_pre[-1,:], np.linspace(t1,t2,tsteps2),args=(V_post,))

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
