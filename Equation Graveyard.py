
#generic forms
I_x = g_x*(V-V_x) = g_x * m_p_x * h_q_x * (V-V_x)
dm_x/dt = (m_x,inf(V)-m_x)/(tau_m_x(V))
dh_x/dt = (h_x,inf(V)-h_x)/tau_h_x(V))
m_x,inf(V) = (1)/(1-e^(-(V-V_0.5)/k_a))
h_x,inf(V) = (1)/(1-e^(V-V_0.5)/k_i))
tau_m_x(V) = tau_h_x(V) = (a)/(1+e^((V-b)/c))+d
C*dV/dt = -I_ion + I_ext
I_ion = (I_SHL1 + I_KVS1 + I_SHK1+ I_IRK + I_KQT3 + I_EGL36 + I_EGL2) + (I_EGL19 + I_UNC2 + I_CCA1) + (I_SLO1/EGL19 + I_SLO1/




#BK Channels Modeling
dm_BK_s/dt = (m_BK_s,inf(V,Ca)-m_BK_s)/(tau_m_BK_s(V,Ca))

m_BK_s,inf = k_plus/(k_plus + k_minus)
tau_m_BK_s = (1)/(k_plus + k_minus)
k_minus = w_minus(V)*f_minus(Ca)
k_plus = w_plus(V)*f_plus(Ca)

w_minus(V) = w_minus_0*e^(-w_yx*V)
w_plus(V) = w_plus_0*e^(-w_xy*V)
f_minus(Ca) = (1)/(1+(Ca/K_yx)^(n_yx))
f_plus(Ca) = (1)/(1+(K_xy/Ca)^(n_xy))

m_BK_s,inf(V, Ca) = (1)/(1+e^(-((V-V_0.5_a)/k_a)))
V_0.5_a = k_a * ((log(w_0_minus/w_0_plus) + log(1 + (K_xy/Ca)^n_xy) - log(1 + (Ca/K_yx)^n_yx)))
k_a - (1)/(w_yx - w_xy)

tau_m_BK_s(V, Ca) = ((e^(w_xy*y))/(w_0_plus)) * (1 + (K_xy/Ca)^(n_xy)) * ((1)/(1 + e^(-(V-V_0.5_a)/k_a)))


#isolated SLO1 and SLO2 currents
I_BK_s = g_BK * m_BK_s * (V-V_k)

m_BK_inf(V, Ca) = ((m_CaV * k_plus_0) * (alpha + beta + k_minus_c))/((k_plus_0 + k_minus_0)*(k_minus_c + alpha) + (beta*k_minus_c))
tau_m_BK(V, Ca) = ((alpha + beta + k_minus_c))/((k_plus_0 + k_minus_0)*(k_minus_c + alpha) + (beta*k_minus_c))
alpha = (m_CaV,inf)/(tau_m_CaV)
beta = tau_-1_m_CaV - alpha

#total current flowing across BK channels within the BK-CaV complex is given by:
I_BK = g_BK * m_BK * h_CaV * (V-V_K)



#SK Channels Modeling
m_KCNL,inf(Ca) = (Ca)/(K_Ca + Ca) #where KCa = 0.33 uM, Hill coefficient is assumed equal to 1
I_KCNL = g_KCNL * m_KCNL * (V-V_K) #total current flowing through SK channels


#Intracellular calcium modeling
[CA2+]_n_i,o = ((i_Ca)/(8*pi*r*D_Ca*F) * exp((-r)/(sqrt((D_Ca)/(k_plus_B * [B]_tot))))) #where D_Ca = 250 um^2s^-1, r = 13 nm, F = 96485 C*mol^-1, k_plus_B = 500 uM^-1s^-1, [B]_tot = 30 uM
#calcium current through a singe open voltage-gated calcium channel:
i_Ca = g_sc * (V-V_Ca) #g_sc = 40 pS, V_Ca = 60 mV
#when calcium channel is closed, [Ca]_n_i,c = 0.05 uM

d[Ca2+]_m_i/dt = (1 - H(V - V_Ca)) * ((-f * alpha * I_Ca) - (([Ca2+]_m_i - [Ca2+]_m_eq)/(tau_Ca)) + H(V-V_Ca)*(-([Ca2+]_m_i - [Ca2+]_m_eq)/(tau_Ca)))
alpha = (1)/(2*V_cell * F)
I_Ca = I_EGL19 + I_UNC2 + I_CCA1
#V_Cell is cell volume, taken from neuromorpho.org
#f = 0.001, tau_Ca = 33 ms, H(V) = heaviside function, V_Ca = 60 mV, [CA2+]_m_eq = [Ca2+]_n_i,c = 0.05 uM




#SHL1
V_0.5 = 11.2 * mV #potentially -6.8
k_a = 14.1 * mV
m_SHL1_inf = lambda V : (1)/(1-e^(-(V-V_0.5)/k_a))

a =   13.8 * ms #potentially 1.4
b = -17.5 * mV
c = 12.9 * mV
d = -3.7 * mV
e = 6.5 * mV
f = 1.9 * ms #potentially 0.2
tau_m_SHL1(V) = ((a)/((e^(-(V-b)/c)) + (e^(-(V-b)/e)))) + f

V_0.5 = -33.1 * mV
k_i = 8.3 * mV
h_f_SHL1,inf(V) = h_s_SHL1,inf(V) = (1)/(1+e^((V-V_0.5)/k_i))

a = 539.2 * ms #potentially 53.9
b = -28.2 * mV
c = 4.9 mV
d = 27.3 * ms #potentially 2.7
tau_f_h_SHL1(V) = ((a)/(1 + e^((V-b)/c))) + d

a = 8422.0 * ms #potentially 842.2
b = -37.7 * mV
c = 6.4 * mV
d = 118.9 * ms #potentially 11.9
tau_s_h_SHL1(V) = ((a)/(1 + e^((V-b)/c))) + d
I_SHL1 = g_SHL1 * (m_SHL1)**3 * (0.7 * h_f_SHL1 + 0.3 * h_s_SHL1) * (V - V_K)


#KVS1
V_0.5 = 57.1 * mV #potentially 27.1
k_a = 25.0 * mV
m_KVS1_inf = (1)/(1-e^(-(V-V_0.5)/k_a))

V_0.5 = 47.3 * mV #potentially 17.3
k_i = 11.1 * mV
h_KVS1_inf = (1)/(1-e^((V-V_0.5)/k_i))

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

I_KVS1 = g_KVS1 * m_KVS1 * h_KVS1 * (V - V_K)


#SHK1

V_0.5 = 20.4 * mV
k_a = 7.7 * mV
m_SHK1_inf = (1)/(1-e^(-(V-V_0.5)/k_a))

a = 26.6 * ms
b = -37.7 * mV
c = 15.8 * mV
d = -37.7 * mV
e = 11.2 * mV
f = 3.8 * ms
tau_m_SHK1 = ((a)/((e^(-(V-b)/c)) + (e^((V-d)/e)))) + f

V_0.5 = -7.0 * mV
k_i = 5.8 * mV
h_SHK1_inf = (1)/(1-e^((V-V_0.5)/k_i))

a = 1400 * ms
tau_h_SHK1 = a
I_SHK1 = g_SHK1 * m_SHK1 * h_SHK1 * (V - V_K)


#KQT3

V_0.5 = -12.8 * mV (#potentially 7.7)
k_a = 15.8 * mV
m_f_KQT3_inf = m_s_KQT3_inf = (1)/(1-e^(-(V-V_0.5)/k_a))

a = 395.3 * ms #potentially 39.5
b = -38.1 * mV
c = 33.6 * mV
tau_f_m_KQT3 = (a)/(1 + ((V+b)/c)**2)

a = 5503.0 * ms #potentially 550.3
b = 5343.3 * ms #potentially 534.5
c = -0.0283 * mV**-1
d = -23.9 * mV
e = 4590 * ms #potentially 459.1
f = -0.0357 * mV**-1
g = 14.2 * mV
tau_s_m_KQT3 = a + (b)/(1 + (10)**(-c * (d-V))) + (e)/(1 + (10)**(-f*(g+V)))

V_0.5 =  -1.1 * mV
k_i = 28.8 * mV
a = 0.5
b = 0.5
w_KQT3_inf = a + (b)/(1 + e^((V-V_0.5)/k_i))

V_0.5 = -45.3 * mV
k_i = 12.3 * mV
a = 0.3
b = 0.7
s_KQT3_inf = a + (b)/(1 + e^((V-V_0.5)/k_i))

a = 0.5 * ms
b = 2.9 * ms
c = -48.1 * mV
d = 48.8 * mV
tau_w_KQT3 = a + (b)/(1 + ((V-c)/d)**2)

a = 500 * ms
tau_s_KQT3 = a

I_KQT3 = g_KQT3 * (0.7*m_f_KQT3 + 0.3 * m_s_KQT3) * w_KQT3 * s_KQT3 * (V-V_K)


#EGL2
V_0.5 = 6.9 * mV
k_a = 14.9 * mV
m_EGL2_inf = (1)/(1-e^(-(V-V_0.5)/k_a))

a = 16.8 * ms #potentially 8.4
b = -122.6 * mV
c = -13.8 * mV
d = 8.1 * ms #potentially 4.1
tau_m_EGL2 = (a)/(1 + e^((V-b)/c)) + d

I_EGL2 = g_EGL2 * m_EGL2 * (V-V_K)


#EGL36
V_0.5 = 63 * mV
k_a = 28.5 * mV
m_f_EGL36_inf = m_m_ELG36_inf = m_s_EGL36_inf = (1)/(1-e^(-(V-V_0.5)/k_a))

a = 13.0 * ms
tau_f_m_EGL36 = a
a = 63.0 * ms
tau_m_m_EGL36 = a
a = 355.0 * ms
tau_s_m_EGL36 = a

I_EGL36 = g_EGL36 * (0.33 * m_f_EGL36 + 0.36 * m_m_EGL36 + 0.39 * m_s_EGL36) * (V-V_K)


#IRK
V_0.5 = -86.5 * mV
k_a = -28.0 * mV
m_IRK_inf = (1)/(1-e^((V-V_0.5)/k_a))

a = 17.1 * ms
b = -17.8 * mV
c = 20.3 * mV
d = -43.4 * mV
e = 11.2 * mV
f = 3.8 * ms
tau_m_IRK = ((a)/((e^(-(V-b)/c)) + (e^((V-d)/e)))) + f

I_IRK = g_IRK * m_IRK * (V-V_K)







#VOLTAGE-GATED CALCIUM CURRENTS


#EGL19
V_0.5 = 5.6 * mV #potentially -4.4
k_a = 7.5 * mV
m_EGL19_inf = (1)/(1-e^(-(V-V_0.5)/k_a))

a = 2.9 * ms
b = 5.2 * mV #potentially -4.8
c = 6.0 * mV
d = 1.9 * ms
e = 1.4 * mV #potentially -8.6
f = 30.0 * mV
g = 2.3 * ms
tau_m_EGL19 = (a * e^(-(((V-b)/c))**2)) + (d * e^(-(((V-e)/f))**2)) + g

V_0.5 = 24.9 * mV #potentially 14.9
k_i = 12.0 * mV
k_b_i = -10.5 * mV #potentially -20.5
V_b_0.5 = 8.1 * mV
a = 1.4
b = 0.1
c = 6.0
d = 0.6
h_EGL19_inf = ((a)/((1+e^(-(V-V_0.5)/k_i))) + b) + ((c)/((1+e^((V-V_b_0.5)/k_b_i))) + d)

a = 0.4
b = 44.6 * ms
c = -23.0 * mV #potentially -33.0
d = 5.0 * mV
e = 36.4 * ms
f = 28.7 * mV #potentially 18.7
g = 3.7 * mV
h = 43.1 * ms
tau_h_EGL19 = (a) * ((b)/((1+e^((V-c)/d))) + (e)/((1+e^((V-f)/g))) + h)
I_EGL19 = g_EGL19 * m_EGL19 * h_EGL19 * (V-V_Ca)


#UNC2
V_0.5 = -12.2 * mV #potentially -37.2
k_a = 4.0 * mV
m_UNC2_inf = (1)/(1-e^(-(V-V_0.5)/k_a))

a = 4.5 * ms
b = -8.2 * mV #potentially -38.2
c = 9.1 * mV
d = 15.4 * mV
e = 0.3 * ms
tau_m_UNC2 = ((a)/((e^(-(V-b)/c)) + (e^((V-b)/d)))) + e

V_0.5 = -52.5 * mV #potentially -77.5
k_i = 5.6 * mV
h_UNC2_inf = (1)/(1-e^((V-V_0.5)/k_i))

a = 83.8 * ms #potentially 142.5
b = 52.9 * mV #potentially 22.9
c = 3.5 * mV
d = 72.1 * ms #potentially 122.6
e = 23.9 * mV #potentially -6.1
f = 3.6 * mV
tau_h_UNC2 = (a)/(1 + e^(-(V-b)/c)) + (d)/(1 + e^((V-e)/f))
I_UNC2 = g_UNC2 * m_UNC2 * h_UNC2 * (V-V_Ca)


#CCA1

V_0.5 = -43.32 * mV #potentially -57.7
k_a = 7.6 * mV #potentially 2.4
m_CCA1_inf = (1)/(1-e^(-(V-V_0.5)/k_a))

V_0.5 = -58.0 * mV #potentially -73.0
k_i = 7.0 * mV #potentially 8.1
h_CCA1_inf = (1)/(1-e^((V-V_0.5)/k_i))

a = 40.0 * ms #potentially 20
b = -62.5 * mV #potentially -92.5
c = -12.6 * mV #potentially 21.1
d = 0.7 * ms #potentially 0.4
tau_m_CCA1 = (a)/(1 + e^((V-b)/c)) + d

a = 280 * ms #potentially 22.4
b = -60.7 * mV #potentially -75.7
c = 8.5 * mV #potentially 9.4
d = 19.8 * ms #potentially 1.6
t_h_CCA1 = (a)/(1 + e^((V-b)/c)) + d


I_CCA1 = g_CCA1 * (m_CCA1)**2 * h_CCA1 * (V-V_Ca)







#CALCIUM-REGULATED POTASSIUM CURRENTS

#SLO1-SLO2
#CaV = EGL19, UNC2

#SLO1
w_yx = 0.013 * mV**-1
w_xy = -0.028 mV**-1
w_minus_0 = 3.15 * ms**-1
w_plus_0 = 0.16 * ms**-1
K_xy = 55.73 * uM
n_xy = 1.30
K_yx = 0.034 * uM
n_yx = 10**-4

#SLO2
w_yx = 0.019 * mV**-1
w_xy = -0.024 * mV**-1
w_minus_0 - 0.87 * ms**-1
w_plus_0 = 0.028 * ms**-1
K_xy = 93.45 * uM
n_xy = 1.84
K_yx = 3294.55 * uM
n_yx = 10**-5


m_BK_inf(V, Ca) = ((m_CaV * k_plus_0) * (alpha + beta + k_minus_c))/((k_plus_0 + k_minus_0)*(k_minus_c + alpha) + (beta*k_minus_c))
tau_m_BK(V, Ca) = ((alpha + beta + k_minus_c))/((k_plus_0 + k_minus_0)*(k_minus_c + alpha) + (beta*k_minus_c))
alpha = (m_CaV,inf)/(tau_m_CaV)
beta = tau_-1_m_CaV - alpha
I_BK = g_BK * m_BK * h_CaV * (V-V_K)


#KCNL
K_Ca = 0.33 * uM
m_KCNL_inf = (Ca)/(K_Ca + Ca)

a = 6.3 * ms
tau_m_KCNL = a

I_KCNL = g_KCNL * m_KCNL * (V-V_K)


#NCA
I_NCA = g_NCA * (V - V_Na)


#LEAK
I_LEAK = g_LEAK * (V - V_L)


#Intracellular calcium calculation
g_sc = 40 * pS
V_Ca = 60 * mV
r = 13 * nm
F = 96485 * C * (mol**-1)
D_Ca = 250 * (mu)*2 * m * (s)**-1
k_plus_B = 500 * (mu) * (M)**-1 * (s)**-1
[B]_tot = 30 * uM
[Ca2+]_n_c_i = 0.05 * uM
V_cell_AWC = 31.16 * uM**3
V_cell_RMD = 5.65 * uM**3
f = 0.001
tau_Ca = 50 * ms
[Ca2+]_m_eq = 0.05 * uM
