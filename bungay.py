# 0  IIf
# 1  IIL
# 2  mIIaf
# 3  mIIaL
# 4  Vf
# 5  VL
# 6  Vaf
# 7  VaL
# 8  VIIf
# 9  VIIL
# 10 VIIaf
# 11 VIIaL
# 12 VIIIf
# 13 VIIIL
# 14 VIIIaf
# 15 VIIIaL
# 16 IXf
# 17 IXL
# 18 IXaf
# 19 IXaL
# 20 Xf
# 21 XL
# 22 Xaf
# 23 XaL
# 24 APCf
# 25 APCL
# 26 PSf
# 27 PSL
# 28 VIIIaif
# 29 VIIIaiL
# 30 Vaif
# 31 VaiL
# 32 PCf
# 33 PCL
# 34 TFL
# 35 TFVIIaL
# 36 TFVIIL
# 37 TFVIIaIXL
# 38 TFVIIaXL
# 39 TFVIIaXaL
# 40 TFVIIXaL
# 41 IXaVIIIaL
# 42 XaVaL
# 43 IXaVIIIaXL
# 44 VXaL
# 45 VIIIXaL
# 46 IIaf
# 47 V_IIaL
# 48 VIII_IIaL
# 49 XaVaIIL
# 50 XaVamIIaL
# 51 XIf
# 52 XI_IIaf
# 53 XIaf
# 54 APCPSL
# 55 APCPSVIIIaL
# 56 TFPIf
# 57 ATf
# 58 IIaATf
# 59 TFPIXaf
# 60 TFPIXaTFVIIaL
# 61 APCPSVaL
# 62 IXaATf
# 63 XaATf
# 64 VIIXaL
# 65 VmIIaL
# 66 VIIImIIaL
# 67 TML
# 68 IIaTML
# 69 IIaTMPCL
# 70 mIIaATL
# 71 XIaIXL
# 72 LIPID


import numpy as np # The Numpy package is used frequently for this work
from scipy.integrate import solve_ivp # The solve_ivp method from Scipy is used to simulate the model
import matplotlib.pyplot as plt # matplotlib is used for plotting the thrombin generation curves in the plotThr method

def setIC(ICvector = np.array([10e-12, 1.45e-6, 2e-8, 1e-8, 1e-8/100, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.45e-6, 2.5e-9]), lipid=6.79e-3):
    """
    Converts a list of factor levels into an initial condition vector for the Bungay model assuming 0 concnetration for all remaining species. \n
    Inputs: \n
    ICvector - List of factor levels in the order TF, II, V, VII, VIIa, VIII, IX, X, XI, AT, TFPI. Leave black for default initial concentrations (10pM of TF) \n
    Outputs: \n
    IC - Vector of initial conditions for the Bungay model
    """
    IC = np.zeros(73);
    IC[[34,0,4,8,10,12,16,20,51,57,56,72]] = np.append(ICvector,lipid); # np.append adds the concentration of lipids onto the input. Indices correspond to factors in initial condition vector
    return IC

def getRates():
    """
    Get the reaction rates for simulating the Bungay model. \n
    Inputs: \n
    None \n
    Outputs: \n
    List of reaction rates
    """
    k1 = np.array([0.5*1e9, 0.005, 0.005*1e9, 0.005, 0.01*1e9, 2.09, 0.34, 0.1*1e9, 32.5, 1.5, 0.05*1e9, 44.8, 15.2, 0.1*1e9, 0.2, 1.0*1e9, 1.0, 0.1*1e9, 10.7, 8.3, 0.1*1e9, 1.0, 0.043, 0.1*1e9, 2.1, 0.023, 0.1*1e9, 6.94, 0.23, 0.1*1e9, 13.8, 0.9, 0.1*1e9, 100.0, 0.1*1e9, 66.0, 13.0, 15.0, 0.05*1e9, 44.8, 15.2, 0.1*1e9, 10.0, 1.43, 0.1*1e9, 1.6, 0.4, 0.1*1e9, 1.6, 0.4, 0.016*1e9, 3.3E-4, 0.01*1e9, 0.0011, 4.9E-7*1e9, 2.3E-6*1e9, 6.83E-5*1e9, 0.1*1e9, 6.94, 1.035, 0.1*1e9, 13.8, 0.9, 1.0*1e9, 0.5, 0.1*1e9, 6.4, 3.6, 6.83E-6*1e9, 0.1*1e9, 0.5, 0.01*1e9, 1.417, 0.183, 1.0]);
    kon = np.array([0.0043*1e9, 0.05*1e9, 0.05*1e9, 0.057*1e9, 0.05*1e9, 0.05*1e9, 0.05*1e9, 0.05*1e9, 0.05*1e9, 0.05*1e9, 0.01*1e9, 0.029*1e9, 0.05*1e9, 0.05*1e9, 0.05*1e9, 0.057*1e9, 0.05*1e9]);
    koff = np.array([1.0, 0.475, 0.145, 0.17, 0.66, 0.227, 0.1, 0.335, 0.115, 0.115, 1.9, 3.3, 3.5, 0.2, 0.335, 0.17, 11.5]);
    nva = np.array([100.0]);
    k=np.concatenate((k1,kon,koff,nva));
    return k

def getThr(k,y,maxt):
    """
    Simulate the Bungay model and return the thrombin generation curve. \n
    Inputs: \n
    k - List of reaction rates, \n
    y - Initial condition vector from setIC, \n
    maxt - Simulation time specified in seconds (1200s = 20min) \n
    Outputs: \n
    Tuple containing a vector of time points and a vector of thrombin concentrations
    """
    sol = ODESolver(k,y,maxt);
    return (sol.t,sol.y[46])

def plotThr(k,y,maxt):
    """
    Simulate the Bungay model and plots the thrombin generation curve. \n
    Inputs: \n
    k - List of reaction rates, \n
    y - Initial condition vector from setIC, \n
    maxt - Simulation time specified in seconds (1200s = 20min) \n
    Outputs: \n
    None, plot of thrombin generation curve is displayed.
    """
    (t,thr) = getThr(k,y,maxt);
    plt.plot(t/60,thr/1e-9)
    plt.xlabel('Time (min)')
    plt.ylabel('Thrombin Concentration (nM)')
    plt.title("Bungay Model")
    plt.show()

def test():
    """
    Simulate the Bungay model and plots the thrombin generation curve for default initial concentrations (Requires no inputs). \n
    Inputs: \n
    None \n
    Outputs: \n
    None, plot of thrombin generation curve is displayed.
    """
    y = setIC();
    k = getRates();
    plotThr(k,y,1200)

def ODESolver(k,c0,maxt):
    """
    Simulate the Bungay model and returns the ODE solution object. \n
    Inputs: \n
    k - List of reaction rates, \n
    y - Initial condition vector from setIC, \n
    maxt - Simulation time specified in seconds (1200s = 20min) \n
    Outputs: \n
    sol - ODE solution object containing timecourse data of all species
    """
    sol = solve_ivp(ODE, (0,maxt), c0, args=(k, ), method='BDF', rtol=1e-3,atol=1e-14);
    return sol

def ODE(t,y,k):
    """
    Evaluate the ODE of the Bungay model. \n
    Inputs: \n
    t - Current simulation time (not needed for evaluation but required by solve_ivp), \n
    y - Current concentration vector, \n
    k - List of reaction rates \n
    Outputs: \n
    dy - Gradient of the concentration vector
    """
    
    kon = k[75:92];
    koff = k[92:109];
    nva = k[109];
    
    dy = np.zeros(73);
    lipidChange = 0;
    
    # II lipid binding
    rate=(kon[0]*y[0]*y[72]/nva-koff[0]*y[1]);
    dy[[0,1]] += [-rate, rate];
    lipidChange += rate;
    
    # mIIa lipid binding
    rate=(kon[1]*y[2]*y[72]/nva-koff[1]*y[3]);
    dy[[2,3]] += [-rate, rate];
    lipidChange += rate;
    
    # V lipid binding
    rate=(kon[2]*y[4]*y[72]/nva-koff[2]*y[5]);
    dy[[4,5]] += [-rate, rate];
    lipidChange += rate;
    
    # Va lipid binding
    rate=(kon[3]*y[6]*y[72]/nva-koff[3]*y[7]);
    dy[[6,7]] += [-rate, rate];
    lipidChange += rate;
    
    # VII lipid binding
    rate=(kon[4]*y[8]*y[72]/nva-koff[4]*y[9]);
    dy[[8,9]] += [-rate, rate];
    lipidChange += rate;
    
    # VIIa lipid binding
    rate=(kon[5]*y[10]*y[72]/nva-koff[5]*y[11]);
    dy[[10,11]] += [-rate, rate];
    lipidChange += rate;
    
    # VIII lipid binding
    rate=(kon[6]*y[12]*y[72]/nva-koff[6]*y[13]);
    dy[[12,13]] += [-rate, rate];
    lipidChange += rate;
    
    # VIIIa lipid binding
    rate=(kon[7]*y[14]*y[72]/nva-koff[7]*y[15]);
    dy[[14,15]] += [-rate, rate];
    lipidChange += rate;
    
    # IX lipid binding
    rate=(kon[8]*y[16]*y[72]/nva-koff[8]*y[17]);
    dy[[16,17]] += [-rate, rate];
    lipidChange += rate;
    
    # IXa lipid binding
    rate=(kon[9]*y[18]*y[72]/nva-koff[9]*y[19]);
    dy[[18,19]] += [-rate, rate];
    lipidChange += rate;
    
    # X lipid binding
    rate=(kon[10]*y[20]*y[72]/nva-koff[10]*y[21]);
    dy[[20,21]] += [-rate, rate];
    lipidChange += rate;
    
    # Xa lipid binding
    rate=kon[11]*y[22]*y[72]/nva-koff[11]*y[23];
    dy[[22,23]] += [-rate, rate];
    lipidChange += rate;
    
    # APC lipid binding
    rate=(kon[12]*y[24]*y[72]/nva-koff[12]*y[25]);
    dy[[24,25]] += [-rate, rate];
    lipidChange += rate;
    
    # PS lipid binding
    rate=(kon[13]*y[26]*y[72]/nva-koff[13]*y[27]);
    dy[[26,27]] += [-rate, rate];
    lipidChange += rate;
    
    # VIIIai lipid binding
    rate=(kon[14]*y[28]*y[72]/nva-koff[14]*y[29]);
    dy[[28,29]] += [-rate, rate];
    lipidChange += rate;
    
    # Vai lipid binding
    rate=(kon[15]*y[30]*y[72]/nva-koff[15]*y[31]);
    dy[[30,31]] += [-rate, rate];
    lipidChange += rate;
    
    # PC lipid binding
    rate=(kon[16]*y[32]*y[72]/nva-koff[16]*y[33]);
    dy[[32,33]] += [-rate, rate];
    lipidChange += rate;
    
    dy[72] -= nva*lipidChange;
    
    # TF + VIIa <-> TF:VIIa
    rate=(k[0]*y[34]*y[11]-k[1]*y[35]);
    dy[[34,11,35]] += [-rate, -rate, rate];
    
    # TF + VII <-> TF:VII
    rate=(k[2]*y[34]*y[9]-k[3]*y[36]);
    dy[[9,34,36]] += [-rate, -rate, rate];
    
    # IX + TF:VIIa <-> IX:TF:VIIa
    rate=(k[4]*y[35]*y[17]-k[5]*y[37]);
    dy[[17,35,37]] += [-rate, -rate, rate];
    
    # IX:TF:VIIa -> IXa + TF:VIIa
    rate=k[6]*y[37];
    dy[[37,35,19]] += [-rate, rate, rate];
    
    # X + TF:VIIa <-> X:TF:VIIa
    rate=(k[7]*y[35]*y[21]-k[8]*y[38]);
    dy[[21,35,38]] += [-rate, -rate, rate];
    
    # X:TF:VIIa -> Xa:TF:VIIa
    rate=k[9]*y[38];
    dy[[38,39]] += [-rate, rate];
    
    # Xa:TF:VIIa -> Xa + TF:VIIa
    rate=k[74]*y[39];
    dy[[39,35,23]] += [-rate, rate, rate];
    
    # TF:VII + Xa <-> Xa:TF:VII
    rate=(k[10]*y[36]*y[23]-k[11]*y[40]);
    dy[[23,36,40]] += [-rate, -rate, rate];
    
    # Xa:TF:VII -> Xa + TF:VIIa
    rate=k[12]*y[40];
    dy[[40,23,35]] += [-rate, rate, rate];
    
    # IXa + VIIIa <-> IXa:VIIIa
    rate=(k[13]*y[19]*y[15]-k[14]*y[41]);
    dy[[15,19,41]] += [-rate, -rate, rate];
    
    # Xa + Va <-> Xa:Va
    rate=(k[15]*y[23]*y[7]-k[16]*y[42]);
    dy[[7,23,42]] += [-rate, -rate, rate];
    
    # X + IXa:VIIIa <-> X:IXa:VIIIa
    rate=(k[17]*y[41]*y[21]-k[18]*y[43]);
    dy[[21,41,43]] += [-rate, -rate, rate];
    
    # X:IXa:VIIIa -> Xa + IXa:VIIIa
    rate=k[19]*y[43];
    dy[[43,23,41]] += [-rate, rate, rate];
    
    # V + Xa <-> V:Xa
    rate=(k[20]*y[5]*y[23]-k[21]*y[44]);
    dy[[5,23,44]] += [-rate, -rate, rate];
    
    # V:Xa -> Va + Xa
    rate=k[22]*y[44];
    dy[[44,7,23]] += [-rate, rate, rate];
    
    # Xa + VIII <-> Xa:VIII
    rate=k[23]*y[13]*y[23]-k[24]*y[45];
    dy[[13,23,45]] += [-rate, -rate, rate];
    
    # Xa:VIII -> Xa + VIIIa
    rate=k[25]*y[45];
    dy[[45,15,23]] += [-rate, rate, rate];
    
    # V + IIa <-> V:IIa
    rate=(k[26]*y[5]*y[46]-k[27]*y[47]);
    dy[[5,46,47]] += [-rate, -rate, rate];
    
    # V:IIa -> Va:IIa
    rate=k[28]*y[47];
    dy[[47,7,46]] += [-rate, rate, rate];
    
    # VIII + IIa <-> VIII:IIa
    rate=(k[29]*y[13]*y[46]-k[30]*y[48]);
    dy[[13,46,48]] += [-rate, -rate, rate];
    
    # VIII:IIa -> VIIIa + IIa
    rate=k[31]*y[48];
    dy[[48,15,46]] += [-rate, rate, rate];
    
    # Xa:Va + II <-> Xa:Va:II
    rate=(k[32]*y[42]*y[1]-k[33]*y[49]);
    dy[[1,42,49]] += [-rate, -rate, rate];
    
    # Xa:Va + mIIa <-> Xa:Va:mIIa
    rate=(k[34]*y[42]*y[3]-k[35]*y[50]);
    dy[[3,42,50]] += [-rate, -rate, rate];
    
    # Xa:Va:II -> Xa:Va:mIIa
    rate=k[36]*y[49];
    dy[[49,50]] += [-rate, rate];
    
    # Xa:Va:mIIa -> Xa:Va + IIa
    rate=k[37]*y[50];
    dy[[50,42,46]] += [-rate, rate, rate];
    
    # VII + Xa <-> VII:Xa
    rate=(k[38]*y[9]*y[23]-k[39]*y[64]);
    dy[[9,23,64]] += [-rate, -rate, rate];
    
    # VII:Xa -> VIIa + Xa
    rate=k[40]*y[64];
    dy[[64,11,23]] += [-rate, rate, rate];
    
    # XI + IIa <-> XI:IIa
    rate=(k[41]*y[51]*y[46]-k[42]*y[52]);
    dy[[46,51,52]] += [-rate, -rate, rate];
    
    # XI:IIa -> XIa + IIa
    rate=k[43]*y[52];
    dy[[52,46,53]] += [-rate, rate, rate];
    
    # PCa:PS + VIIIa <-> PCa:PS:VIIIa
    rate=(k[44]*y[54]*y[15]-k[45]*y[55]);
    dy[[15,54,55]] += [-rate, -rate, rate];
    
    # PCa:PS:VIIIa -> PCa:PS + VIIIai
    rate=k[46]*y[55];
    dy[[55,29,54]] += [-rate, rate, rate];
    
    # PCa:PS + Va <-> PCa:PS:Va
    rate=(k[47]*y[54]*y[7]-k[48]*y[61]);
    dy[[7,54,61]] += [-rate, -rate, rate];
    
    # PCa:PS:Va -> PCa:PS + Vai
    rate=k[49]*y[61];
    dy[[61,31,54]] += [-rate, rate, rate];
    
    # TFPI + Xa <-> Xa:TFPI
    rate=(k[50]*y[56]*y[22]-k[51]*y[59]);
    dy[[22,56,59]] += [-rate, -rate, rate];
    
    # Xa:TFPI + TF:VIIa <-> Xa:TFPI:TF:VIIa
    rate=(k[52]*y[59]*y[35]-k[53]*y[60]);
    dy[[35,59,60]] += [-rate, -rate, rate];
    
    # IXa + AT -> IXa:AT
    rate=k[54]*y[18]*y[57];
    dy[[18,57,62]] += [-rate, -rate, rate];
    
    # Xa + AT -> Xa:AT
    rate=k[55]*y[22]*y[57];
    dy[[22,57,63]] += [-rate, -rate, rate];
    
    # IIa + AT -> IIa:AT
    rate=k[56]*y[46]*y[57];
    dy[[46,57,58]] += [-rate, -rate, rate];
    
    # V + mIIa <-> V:mIIa
    rate=(k[57]*y[5]*y[3]-k[58]*y[65]);
    dy[[3,5,65]] += [-rate, -rate, rate];
    
    # V:mIIa -> Va + mIIa
    rate=k[59]*y[65];
    dy[[65,3,7]] += [-rate, rate, rate];
    
    # VIII + mIIa <-> VIII:mIIa
    rate=(k[60]*y[13]*y[3]-k[61]*y[66]);
    dy[[3,13,66]] += [-rate, -rate, rate];
    
    # VIII:mIIa -> VIIIa + mIIa
    rate=k[62]*y[66];
    dy[[66,3,15]] += [-rate, rate, rate];
    
    # IIa + TM <-> IIa:TM
    rate=(k[63]*y[46]*y[67]-k[64]*y[68]);
    dy[[46,67,68]] += [-rate, -rate, rate];
    
    # IIa:TM + PC <-> IIa:TM:PC
    rate=(k[65]*y[68]*y[33]-k[66]*y[69]);
    dy[[33,68,69]] += [-rate, -rate, rate];
    
    # IIa:TM:PC -> IIa:TM + PCa
    rate=k[67]*y[69];
    dy[[69,68,25]] += [-rate, rate, rate];
    
    # mIIa + AT -> mIIa:AT
    rate=k[68]*y[2]*y[57];
    dy[[2,57,70]] += [-rate, -rate, rate];
    
    # PCa + PS <-> PCa:PS
    rate=(k[69]*y[25]*y[27]-k[70]*y[54]);
    dy[[25,27,54]] += [-rate, -rate, rate];
    
    # XIa + IX <-> XIa:IX
    rate=(k[71]*y[53]*y[17]-k[72]*y[71]);
    dy[[17,53,71]] += [-rate, -rate, rate];
    
    # XIa:IX -> XIa + IXa
    rate=k[73]*y[71];
    dy[[71,19,53]] += [-rate, rate, rate];
    
    return dy













