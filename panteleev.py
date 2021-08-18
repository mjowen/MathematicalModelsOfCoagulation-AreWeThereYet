# 0  TFVIIa
# 1  TFVII
# 2  TF
# 3  VIIa
# 4  VII
# 5  IXa
# 6  IX
# 7  Xa
# 8  X
# 9  IIa
# 10 II
# 11 Ia
# 12 I
# 13 VIIIa
# 14 VIII
# 15 Va
# 16 V
# 17 XIa
# 18 XI
# 19 AT
# 20 TFPI
# 21 XaTFPI
# 22 APC
# 23 PC

import numpy as np # The Numpy package is used frequently for this work
from scipy.integrate import solve_ivp # The solve_ivp method from Scipy is used to simulate the model
import matplotlib.pyplot as plt # matplotlib is used for plotting the thrombin generation curves in the plotThr method

def setIC(ICvector = np.array([10e-12, 1.45e-6, 2e-8, 1e-8, 1e-8/100, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.45e-6, 2.5e-9])):
    """
    Converts a list of factor levels into an initial condition vector for the Panteleev model assuming 0 concnetration for all remaining species. \n
    Inputs: \n
    ICvector - List of factor levels in the order TF, II, V, VII, VIIa, VIII, IX, X, XI, AT, TFPI. Leave black for default initial concentrations (10pM of TF) \n
    Outputs: \n
    IC - Vector of initial conditions for the Panteleev model
    """
    IC = np.zeros(73);
    IC[[2,10,16,4,3,14,6,8,18,19,20]] = ICvector; # Indices correspond to factors in initial condition vector
    return IC

def getRates():
    """
    Get the reaction rates for simulating the Panteleev model. \n
    Inputs: \n
    None \n
    Outputs: \n
    List of reaction rates
    """
    k1 = np.array([4.2*1e9/60, 0.0014*1e9/60, 0.4*1e9/60, 15/60, 5.8/60, 435/60, 0.06/60, 6350/60, 0.052*1e9/60, 45*1e9*1e9/60, 1.44/60, 5040/60, 54/60, 14/60, 0.03*1e9*1e9/60, 2e-05*1e9/60]);
    km = np.array([1.1/60, 0.02/60, 770/60]);
    K = np.array([210/1e9, 200/1e9, 238/1e9, 230, 1216, 278, 1655, 7200/1e9, 147/1e9, 71.7/1e9, 2.57/1e9, 1.5/1e9, 150/1e9, 0.118/1e9, 200/1e9, 320/1e9, 470/1e9, 2.9/1e9]);
    n = np.array([260, 750, 16000, 2700]);
    h = np.array([0.44*1e9/60, 6*1e9/60, 8.2e-06*1e9/60, 0.00015*1e9/60, 4e-05*1e9/60, 1.36e-05*1e9/60, 0.0012*1e9/60, 2.2e-05*1e9/60, 0.00041*1e9/60, 0.0001*1e9/60, 3e-06*1e9/60, 0.00037*1e9/60, 6.3e-05*1e9/60, 0.35/60, 7.7*1e9/60, 1.9e-05*1e9/60, 2.6e-05*1e9/60, 6e-06*1e9/60, 0.0054*1e9/60, 0.00014*1e9/60, 6e-06*1e9/60, 6e-06*1e9/60, 7e-07*1e9/60, 0.00039*1e9/60]);
    i = np.array([3000/1e9, 40000/1e9, 975/1e9, 1400/1e9, 79/1e9, 1700/1e9, 323/1e9]);
    p = np.array([7.5e-5/1e9]);
    k = np.concatenate((k1, km, K, n, h, i, p))
    return k

def getThr(k,y,maxt):
    """
    Simulate the Panteleev model and return the thrombin generation curve. \n
    Inputs: \n
    k - List of reaction rates, \n
    y - Initial condition vector from setIC, \n
    maxt - Simulation time specified in seconds (1200s = 20min) \n
    Outputs: \n
    Tuple containing a vector of time points and a vector of thrombin concentrations
    """
    sol = ODESolver(k,y,maxt);
    return (sol.t,sol.y[9])

def plotThr(k,y,maxt):
    """
    Simulate the Panteleev model and plots the thrombin generation curve. \n
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
    plt.title("Panteleev Model")
    plt.show()

def test():
    """
    Simulate the Panteleev model and plots the thrombin generation curve for default initial concentrations (Requires no inputs). \n
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
    Simulate the Panteleev model and returns the ODE solution object. \n
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
    Evaluate the ODE of the Panteleev model. \n
    Inputs: \n
    t - Current simulation time (not needed for evaluation but required by solve_ivp), \n
    y - Current concentration vector, \n
    k - List of reaction rates \n
    Outputs: \n
    dy - Gradient of the concentration vector
    """
    
    p = k[72];
    i = k[65:72];
    h = k[41:65];
    n = k[37:41];
    K = k[19:37];
    km = k[16:19]; #k minus (reverse) rather than Km
    k = k[0:16];
    
    dy = np.zeros(73);
    
    
    IXaBF = y[5]*p*n[0]/(K[10]+y[5]);
    VIIaTFF = y[0]/(1+y[6]/K[0]+y[8]/K[2]);
    XaVIIaTF = k[5]/(K[2]*km[2])*y[8]*VIIaTFF;
    VaB = y[15]*p*n[3]/(K[17]+y[15]);
    XaVaB = y[7]*VaB/(K[13]*(1+i[6]/K[14]+y[7]/K[13])+VaB);
    VaBF = VaB-XaVaB;
    XaF = y[7]-XaVaB;
    IIB = y[10]*p*n[2]/(K[16]*(1+y[8]/K[15]+y[10]/K[16]));
    IIaF = y[9]/(1+(y[11]+y[12])/K[7]);
    XB = y[8]*p*n[2]/(K[15]*(1+y[8]/K[15]+y[10]/K[16]));
    VIIIaBF = y[13]*p*n[1]/((K[11]+y[13])*(1+XB/(p*K[6])*(1+i[6]/K[12])));
    
    # Mass Action Law Reactions
    # TF + VII <-> TF:VII
    rate = k[0]*y[2]*y[4]-km[0]*y[1];
    dy[[4,2,1]] += [-rate, -rate, rate];
    
    # TF + VIIa <-> TF:VIIaF
    rate = k[0]*y[2]*y[3]-km[0]*VIIaTFF;
    dy[[3,2,0]] += [-rate, -rate, rate];
    
    # TF:VII + IIaF -> TF:VIIa + IIaF
    rate = k[1]*y[1]*IIaF;
    dy[[1,0]] += [-rate, rate];
    
    # TF:VII + XaF -> TF:VIIa + XaF
    rate = k[2]*y[1]*XaF;
    dy[[1,0]] += [-rate, rate];
    
    # VII + IIaF -> VIIa + IIaF
    rate = k[1]*y[4]*IIaF;
    dy[[4,3]] += [-rate, rate];
    
    # TF:VIIaF + Xa:TFPI -> TF:VIIa:Xa:TFPI
    rate = h[0]*VIIaTFF*y[21];
    dy[[0,21]] += [-rate, -rate];
    
    # TF:VIIa:Xa + TFPI -> TF:VIIa:Xa:TFPI
    rate = h[1]*XaVIIaTF*y[20];
    dy[[0,20]] += [-rate, rate];
    
    # IX + TF:VIIaF -> IXa + TF:VIIa
    rate = k[3]/K[0]*y[6]*VIIaTFF;
    dy[[6,5]] += [-rate, rate];
    
    # IXa + AT -> IXa:AT
    rate = h[2]*y[5]*y[19];
    dy[[5,19]] += [-rate, -rate];
    
    # X + TF:VIIaF -> Xa + TF:VIIa
    rate = k[5]/K[2]*y[8]*VIIaTFF;
    dy[[8,7]] += [-rate, rate];
    
    # XaF + TFPI <-> Xa:TFPI
    rate = k[8]*XaF*y[20]-km[1]*y[21];
    dy[[7,20,21]] += [-rate, -rate, rate];
    
    # XaF + AT -> Xa:AT
    rate = h[3]*XaF*y[19];
    dy[[7,19]] += [-rate, rate];
    
    # XaF + a2M -> Xai + a2M
    rate = h[4]*XaF*i[0];
    dy[7] -= rate; 
    
    # XaF + a1AT -> Xai + a1AT
    rate = h[5]*XaF*i[1];
    dy[7] -= rate; 
    
    # XaF + PCI -> Xai + PCI
    rate = h[6]*XaF*i[4];
    dy[7] -= rate;
    
    # XaVaB + AT -> Xa:AT + Va
    rate = h[7]*XaVaB*y[19];
    dy[[7,19]] += [-rate, -rate];
    
    # II + XaF -> IIa + Xa
    rate = k[9]*p*y[10]*XaF;
    dy[[10,9]] += [-rate, rate];
    
    # IIB + XaVaB -> IIa + XaVa
    rate = k[10]/p*IIB*XaVaB;
    dy[[10,9]] += [-rate, rate];
    
    # IIaF + AT -> IIa:AT
    rate = h[8]*IIaF*y[19];
    dy[[9,19]] += [-rate, -rate];
    
    # IIaF + a2M -> IIai + a2M
    rate = h[9]*IIaF*i[0];
    dy[9] -= rate;
    
    # IIaF + a1AT -> IIai + a1AT
    rate = h[10]*IIaF*i[1];
    dy[9] -= rate;
    
    # IIaF + PCI -> IIai + PCI
    rate = h[11]*IIaF*i[4];
    dy[9] -= rate;
    
    # IIaF + hep -> IIai + hep
    rate = h[12]*IIaF*i[3];
    dy[9] -= rate;
    
    # I + IIaF -> Ia + IIa
    rate = k[11]/K[7]*y[12]*IIaF;
    dy[[12,11]] += [-rate, rate];
    
    # VIIIa -> VIIIai
    rate = h[13]*y[13];
    dy[13] -= rate;
    
    # VaBF + PCa -> Vai + PCa
    rate = h[14]*VaBF*y[22];
    dy[15] -= rate;
    
    # XI + IIaF -> XIa + IIa
    rate = k[14]*p*y[18]*IIaF;
    dy[[18,17]] += [-rate, rate];
    
    # XIa + AT -> XIa:AT
    rate = h[15]*y[17]*y[19];
    dy[[17,19]] += [-rate, -rate];
    
    # XIa + a2AP -> XIa:a2AP
    rate = h[16]*y[17]*i[2];
    dy[17] -= rate;
    
    # XIa + a1AT -> XIa:a1AT
    rate = h[17]*y[17]*i[1];
    dy[17] -= rate;
    
    # XIa + PCI -> XIai + a2AP
    rate = h[18]*y[17]*i[4];
    dy[17] -= rate;
    
    # XIa + C1inh -> XIai + C1inh
    rate = h[19]*y[17]*i[5];
    dy[17] -= rate;
    
    # PC + IIaF -> PCa + IIa
    rate = k[15]*y[23]*IIaF;
    dy[[23,22]] += [-rate, rate];
    
    # PCa + a2M -> PCai + a2M
    rate = h[20]*y[22]*i[0];
    dy[22] -= rate;
    
    # PCa + a2AP -> PCai + a2AP
    rate = h[21]*y[22]*i[2];
    dy[22] -= rate;
    
    # PCa + a1AT -> PCai + a1AT
    rate = h[22]*y[22]*i[1];
    dy[22] -= rate;
    
    # PCa + PCI -> PCai + PCI
    rate = h[23]*y[22]*i[4];
    dy[22] -= rate;
    
    
    # Enzyme Reactions
    # IX by XIa
    rate = k[4]*y[6]*y[17]/(K[1]+y[6]);
    dy[[6,5]] += [-rate, rate];
    
    # XB by IXaBF
    rate = k[6]*IXaBF*XB/(p*K[3]);
    dy[[8,7]] += [-rate, rate];
         
    # XB by IXaBF*VIIIaBF
    rate = k[7]*IXaBF*VIIIaBF*XB/(p**2*K[4]*K[5]);
    dy[[8,7]] += [-rate, rate];
         
    # VIII by IIaF
    rate = k[12]*IIaF*y[14]/(K[8]+IIaF);
    dy[[14,13]] += [-rate, rate];
    
    # V by IIaF
    rate = k[13]*IIaF*y[16]/(K[9]+IIaF);
    dy[[16,15]] += [-rate, rate];
    
    return dy

















