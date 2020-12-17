import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

#I_x = g_x*(V-V_x) = g_x * m_p_x * h_q_x * (V-V_x)
def KVS1(y, t, voltage):
    '''
        mV = 1e-3; ms = 1e-3
        # updated version
        m_v_05 = 11.2 * mV #potentially -6.8
        k_a = 14.1 * mV
        m_KVS1_inf = lambda V : (1)/(1+np.exp(-(V-m_v_05)/k_a))

        a = 13.8 * ms #potentially 1.4
        b = -17.5 * mV
        c = 12.9 * mV
        d = -3.7 * mV
        tau_m_KVS1 = lambda V: a/(1 + np.exp((V-b)/c)) + d

        h_v_05 = -33.1 * mV
        k_i = 8.3 * mV
        h_KVS1_inf = lambda V: (1)/(1+np.exp((V-h_v_05)/k_i))

        a = 539.2 * ms #potentially 53.9
        b = -28.2 * mV
        c = 4.9 * mV
        d = 27.3 * ms #potentially 2.7
        tau_h_KVS1 = lambda V: ((a)/(1 + np.exp((V-b)/c))) + d

        m_KVS1 = y[0]
        h_f_KVS1 = y[1]
        h_s_KVS1 = y[2]
    '''
    V_0_5 = 57.1 * mV #potentially 27.1
    k_a = 25.0 * mV
    m_KVS1_inf = (1)/(1-np.exp(-(V-V_0_5)/k_a))

    V_0_5 = 47.3 * mV #potentially 17.3
    k_i = 11.1 * mV
    h_KVS1_inf = (1)/(1-np.exp((V-V_0_5)/k_i))

    a = 30.0 * ms #potentially 3.0
    b = 18.1 * mV
    c = -20 * mV
    d = 1.0 * mV #potentially 0.1
    tau_m_KVS1 = ((a)/(1 + e^((V-b)/c))) + d

    a = 88.5 * ms #potentially  8.9
    b = 50.0 * mV
    c = -15.0 * mV
    d = 53.4 * ms #potentially 5.3
    tau_h_KVS1 = ((a)/(1 + e^((V-b)/c))) + d
    old_eqs = '''# m-related equations
    m_k_a = 14.1 #mV
    m_v_05 = 11.2 #mV # or -6.8 calibrated
    #m_v_05 = -6.8
    m_KVS1_inf = lambda V : 1/(1+np.exp(-(V-m_v_05)/m_k_a))
    m_a = 13.8; m_b = -17.5; m_c = 12.9; m_d = -3.7; m_e = 6.5; m_f = 1.9 #m_f could be 0.2 calibrated
    tau_m_KVS1 = lambda V : ((m_a)/((np.exp(-(V-m_b)/m_c)) + (np.exp(-(V-m_d)/m_e)))) + m_f

    print(tau_m_KVS1(-10*1e-3))

    # h-related equations
    tau_h_a_f = 539.2/10; tau_h_b_f = -28.2; tau_h_c_f = 4.9; tau_h_d_f = 27.3 # tau_f_h_a could be 53.9 calibrated
    tau_h_a_s = 8422.0/10; tau_h_b_s = -37.7; tau_h_c_s = 6.4; tau_h_d_s = 118.9/10 # a could be 842.2; d could be 11.9 calibrated
    h_v_05 = -33.1
    h_k_i  = 8.3
    h_f_KVS1_inf = lambda V : (1)/(1+np.exp((V-h_v_05)/h_k_i))
    h_s_KVS1_inf = lambda V : (1)/(1+np.exp((V-h_v_05)/h_k_i))
    tau_h_f_KVS1 = lambda V : ((tau_h_a_f)/(1 + np.exp((V-tau_h_b_f)/tau_h_c_f))) + tau_h_d_f
    tau_h_s_KVS1 = lambda V : ((tau_h_a_s)/(1 + np.exp((V-tau_h_b_s)/tau_h_c_s))) + tau_h_d_s
    '''

    #print(tau_m_KVS1(-10*1e-3), "expecting 18ms")
    #print(tau_h_f_KVS1(-10*1e-3),"expecting 40ms")
    #print(tau_h_s_KVS1(-10*1e-3),"expecting 200ms")

    dm_KVS1_dt = (m_KVS1_inf(voltage)-m_KVS1)/(tau_m_KVS1(voltage))

    dh_f_KVS1_dt = (h_f_KVS1_inf(voltage)-h_f_KVS1)/tau_h_f_KVS1(voltage)
    dh_s_KVS1_dt = (h_s_KVS1_inf(voltage)-h_s_KVS1)/tau_h_s_KVS1(voltage)

    vv = np.linspace(-100*1e-3,70*1e-3,1000)
    hinf = h_f_KVS1_inf(vv)
    taum = tau_m_KVS1(vv)
    minf = m_KVS1_inf(vv)
    #plt.plot(vv, hinf)
    #plt.plot(vv, minf)
    plt.plot(vv, taum)
    plt.show()


    #I_KVS1 = g_KVS1 * m_3_KVS1 * (0.7 * h_f_KVS1 + 0.3 * h_s_KVS1) * (V - V_K)

    return np.asarray([dm_KVS1_dt, dh_f_KVS1_dt, dh_s_KVS1_dt])


init_values = [0.00154979, 0.99648898, 0.99648898]

t0 = 0
t1 = 100

V_pre = -80e-3
V_post = 0e-3

result_pre = integrate.odeint(KVS1, np.asarray(init_values), np.linspace(t0,t1,(t1-t0)*1000), args=(V_pre,))
print("Steady state:")
print(result_pre[-1,:])
t2 = 200
result_post = integrate.odeint(KVS1, result_pre[-1,:], np.linspace(t1,t2,(t2-t1)*1000),args=(V_post,))

m_KVS1_pre = result_pre[:,0]
h_f_KVS1_pre = result_pre[:,1]
h_s_KVS1_pre = result_pre[:,2]

m_KVS1_post = result_post[:,0]
h_f_KVS1_post = result_post[:,1]
h_s_KVS1_post = result_post[:,2]

if True:
    g_KVS1 = 1.8 * 1e-9
    V_K = -80*1e-3
    I_KVS1_pre = np.transpose(g_KVS1 * np.power(m_KVS1_pre,3) * (0.7 * h_f_KVS1_pre + 0.3 * h_s_KVS1_pre) * (-80*1e-3 - V_K))
    I_KVS1_post = np.transpose(g_KVS1 * np.power(m_KVS1_post,3) * (0.7 * h_f_KVS1_post + 0.3 * h_s_KVS1_post) * (40*1e-3 - V_K))

    print(I_KVS1_pre)
    print()
    print(I_KVS1_post)
    plt.plot(np.linspace(t0,t2,(t2-t0)*1000),np.hstack((I_KVS1_pre, I_KVS1_post)))
    plt.xlabel("Time (ms)")
    plt.ylabel("I_KVS1 (arbitrary units)")
    plt.show()
