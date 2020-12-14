# Test of one ion channel:


eqs = '''
dV_i/dt  = (-g_units*(V_i-E_cell)-I_gap-I_syn + I_Ext)/C_units : volt
I_gap : amp
I_syn : amp
I_Ext : amp
Vj : volt
Vth : volt
ds_j/dt = a_r * phi * (1-s_j)-a_d*s_j : 1
phi = 1/(1+exp(-beta*(V_i-Vth))) : 1
'''

#voltage-gated potassium channels: SHL-1, KVS-1, SHK-1, IRK1-3, KQT-3, EGL-36, EGL-2
#calcium-regulated channels: BK (SLO-1, SLO-2), SK(KCNL-1/4)

#generic forms



#use calibrated numbers in parentheses when available

I_x = g_x*(V-V_x) = g_x * m_p_x * h_q_x * (V-V_x)
dm_SHL1/dt = (m_SHL1_inf(V)-m_SHL1)/(tau_m_SHL1(V))
dh_x/dt = (h_x,inf(V)-h_x)/tau_h_x(V))
m_SHL1 = (1)/(1-e^(-(V-V_0.5)/k_a))
h_f_SHL1 = h_s_SHL1,inf(V) = (1)/(1+e^((V-V_0.5)/k_i))
tau_m_SHL1(V) = ((a)/((e^(-(V-b)/c)) + (e^(-(V-b)/e)))) + f
tau_h_f_SHL1(V) = ((a_f)/(1 + e^((V-b_f)/c_f))) + d_f
tau_h_s_SHL1(V) = ((a_s)/(1 + e^((V-b_s)/c_s))) + d_s
C*dV/dt = -I_ion + I_ext
I_ion = I_SHL1
I_SHL1 = g_SHL1 * m_3_SHL1 * (0.7 * h_f_SHL1 + 0.3 * h_s_SHL1) * (V - V_K)



#SHL1
m_SHL1_inf = lambda V : (1)/(1-e^(-(V-V_0.5)/k_a))
tau_m_SHL1 = lambda V : ((a)/((e^(-(V-b)/c)) + (e^(-(V-b)/e)))) + f
h_f_SHL1 = lambda V : (1)/(1+e^((V-V_0.5)/k_i))
h_s_SHL1 = lambda V : (1)/(1+e^((V-V_0.5)/k_i))
tau_h_f_SHL1 = lambda V : ((a_f)/(1 + e^((V-b_f)/c_f))) + d_f
tau_h_s_SHL1 = lambda V : ((a_s)/(1 + e^((V-b_s)/c_s))) + d_s
I_SHL1 = g_SHL1 * (m_SHL1)**3 * (0.7 * h_f_SHL1 + 0.3 * h_s_SHL1) * (V - V_K)



# updated version

V_0.5 = 11.2 * mV #potentially -6.8
k_a = 14.1 * mV
m_SHL1_inf = lambda V : (1)/(1-e^(-(V-V_0.5)/k_a))

a =   13.8 * ms #potentially 1.4
b = -17.5 * mV
c = 12.9 * mV
d = -3.7 * mV
e = 6.5 * mV
f = 1.9 * ms #potentially 0.2
tau_m_SHL1 = lambda V: ((a)/((e^(-(V-b)/c)) + (e^(-(V-b)/e)))) + f

V_0.5 = -33.1 * mV
k_i = 8.3 * mV
h_f_SHL1 = lambda V: (1)/(1+e^((V-V_0.5)/k_i))
h_s_SHL1 = lambda V: (1)/(1+e^((V-V_0.5)/k_i))

a = 539.2 * ms #potentially 53.9
b = -28.2 * mV
c = 4.9 mV
d = 27.3 * ms #potentially 2.7
tau_h_f_SHL1 = lambda V: ((a)/(1 + e^((V-b)/c))) + d

a = 8422.0 * ms #potentially 842.2
b = -37.7 * mV
c = 6.4 * mV
d = 118.9 * ms #potentially 11.9
tau_h_s_SHL1 = lambda V: ((a)/(1 + e^((V-b)/c))) + d
I_SHL1 = g_SHL1 * (m_SHL1)**3 * (0.7 * h_f_SHL1 + 0.3 * h_s_SHL1) * (V - V_K)
