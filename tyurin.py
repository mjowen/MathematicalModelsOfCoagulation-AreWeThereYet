# 0  TF
# 1  TFVII
# 2  TFVIIa
# 3  VIIa
# 4  VII
# 5  Xa
# 6  X
# 7  Va
# 8  V
# 9  XaVa
# 10 XIa
# 11 XI
# 12 IXa
# 13 IX
# 14 VIIIa
# 15 VIII
# 16 IXaVIIIa
# 17 IIa
# 18 II
# 19 TM
# 20 TMIIa
# 21 PCa
# 22 PC
# 23 TFPI
# 24 TFPIXa
# 25 ATIII
# 26 C1
# 27 PAI1
# 28 alpha1AT
# 29 alpha2AP
# 30 alpha2M
# 31 PCI
# 32 XIIa

import numpy as np # The Numpy package is used frequently for this work
from scipy.integrate import solve_ivp # The solve_ivp method from Scipy is used to simulate the model
import matplotlib.pyplot as plt # matplotlib is used for plotting the thrombin generation curves in the plotThr method

def setIC(ICvector = np.array([10e-12, 1.45e-6, 2e-8, 1e-8, 1e-8/100, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.45e-6, 2.5e-9]), includeExtras = False):
    """
    Converts a list of factor levels into an initial condition vector for the Tyurin model assuming 0 concentration for all remaining species. \n
    Inputs: \n
    ICvector - List of factor levels in the order TF, II, V, VII, VIIa, VIII, IX, X, XI, AT, TFPI. Leave black for default initial concentrations (10pM of TF), \n
    includeExtras - Boolean variable to determine whether or not the following inhibitors are included in the concentration vector (at their default concentrations): PCI, alpha1-AT, alpha2-M, alpha2-AP, PAI-1 and C1-inh. This variable should match the similar variable in getRates(). Default is False. \n
    Outputs: \n
    IC - Vector of initial conditions for the Tyurin model
    """
    IC = np.zeros(33);
    IC[[0,18,8,4,3,15,13,6,11,25,23]] = ICvector; # Indices correspond to factors in initial condition vector
    if includeExtras:
        IC[[31, 28, 30, 29, 27, 26]] = [7e-8, 4e-5, 3.25e-6, 9.5e-7, 4.6e-10, 2.1e-6];
    return IC

def getRates(includeExtras = False):
    """
    Get the reaction rates for simulating the Tyurin model. \n
    Inputs: \n
    includeExtras - Only relevant for Panteleev model. Boolean variable to determine whether or not the following inhibitors are included in the rate vector (at their default concentrations): PCI, alpha1-AT, alpha2-M, alpha2-AP, PAI-1 and C1-inh. This variable should match the similar variable in setIC(). Default is False. \n    
    Outputs: \n
    List of reaction rates
    """
    k = np.array([1e4*1e6/60, 1e4*1e6/60, 3*1e6/60, 2*1e6/60, 6.4e2*1e6/60, 0.027*1e6/60, 4.25e-1*1e6/60, 4.7e-3*1e6/60, 2.93e-2*1e6/60, 1*1e6/60, 1.88e-1*1e6/60, 1.57e-2*1e6/60, 9.6e2*1e6/60, 1.57e-2*1e6/60, 0.1*1e6/60, 2.94e-2*1e6/60, 3e-2*1e6/60, 1e-3*1e6/60, 4e-3*1e6/60, 1e-2*1e6/60, 3e-2*1e6/60, 12.6*1e6/60, 30*1e6/60, 60*1e6/60, 0.15*1e6/60, 6e-4*1e6/60, 2e-3/60, 2e-4/60, 21/60, 0.05/1e6, 86/60, 0.05/1e6, 8/60, 0.05/1e6, 75/60, 0.35/1e6, 0.011/60, 0.009/1e6, 42/60, 0.1/1e6, 4e-2/60, 1/1e6, 1500/60, 0.16/1e6, 0.147/60, 0.25/1e6, 108/60, 0.22/1e6, 2.25/60, 5.8e-2/1e6, 1700/60, 1.03/1e6, 14/60, 7.17e-2/1e6, 2.6/60, 1.04e-2/1e6, 3/60, 0.05/1e6, 39.6/60, 9.3e-3/1e6, 21.6/60, 2e-2/1e6, 5300/60, 5.9/1e6, 24/60, 2e-2/1e6]);
    return k

def getThr(k,y,maxt):
    """
    Simulate the Tyurin model and return the thrombin generation curve. \n
    Inputs: \n
    k - List of reaction rates, \n
    y - Initial condition vector from setIC, \n
    maxt - Simulation time specified in seconds (1200s = 20min) \n
    Outputs: \n
    Tuple containing a vector of time points and a vector of thrombin concentrations
    """
    sol = ODESolver(k,y,maxt);
    return (sol.t,sol.y[17])

def plotThr(k,y,maxt):
    """
    Simulate the Tyurin model and plots the thrombin generation curve. \n
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
    plt.title("Tyurin Model")
    plt.show()

def test():
    """
    Simulate the Tyurin model and plots the thrombin generation curve for default initial concentrations (Requires no inputs). \n
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
    Simulate the Tyurin model and returns the ODE solution object. \n
    Inputs: \n
    k - List of reaction rates, \n
    y - Initial condition vector from setIC, \n
    maxt - Simulation time specified in seconds (1200s = 20min) \n
    Outputs: \n
    sol - ODE solution object containing timecourse data of all species
    """
    sol = solve_ivp(ODE, (0,maxt), c0, args=(k, ), method='BDF', rtol=1e-4,atol=1e-16);
    return sol

def ODE(t,y,k):
    """
    Evaluate the ODE of the Tyurin model. \n
    Inputs: \n
    t - Current simulation time (not needed for evaluation but required by solve_ivp), \n
    y - Current concentration vector, \n
    k - List of reaction rates \n
    Outputs: \n
    dy - Gradient of the concentration vector
    """
    dy = np.zeros(33);
    
    # Enzyme Reactions
    # 1 XIIa activating XI
    rate = k[28]*y[11]*y[32]/( k[29] +y[11] );
    dy[[11,10]] += [-rate, rate];
    
    # 2 IIa activating XI
    rate = k[30]*y[11]*y[17]/( k[31]* (1+y[8]/k[53]+y[15]/k[61]) +y[11] );
    dy[[11,10]] += [-rate, rate];
    
    # 3 XIa activating XI
    rate = k[32]*y[11]*y[10]/( k[33]* (1+y[13]/k[35]) +y[11] );
    dy[[11,10]] += [-rate, rate];
    
    # 4 XIa activating IX
    rate = k[34]*y[10]*y[13]/( k[35]* (1+y[11]/k[33]) +y[13] );
    dy[[13,12]] += [-rate, rate];
    
    # 5 VIIa activating IX
    rate = k[36]*y[3]*y[13]/( k[37]* (1+y[6]/k[45]) +y[13] );
    dy[[13,12]] += [-rate, rate];
    
    # 6 TFVIIa activating IX
    rate = k[38]*y[2]*y[13]/( k[39]* (1+y[6]/k[47]) +y[13] );
    dy[[13,12]] += [-rate, rate];
    
    # 7 IXa activating X
    rate = k[40]*y[6]*y[12]/( k[41] +y[6] );
    dy[[6,5]] += [-rate, rate];
    
    # 8 IXaVIIIa activating X
    rate = k[42]*y[6]*y[16]/( k[43] +y[6] );
    dy[[6,5]] += [-rate, rate];
    
    # 9 VIIa activating X
    rate = k[44]*y[3]*y[6]/( k[45]* (1+y[13]/k[37]) +y[6] );
    dy[[6,5]] += [-rate, rate];
    
    # 10 TFVIIa activating X
    rate = k[46]*y[2]*y[6]/( k[47]* (1+y[13]/k[39]) +y[6] );
    dy[[6,5]] += [-rate, rate];
    
    # 11 Xa activating II
    rate = k[48]*y[5]*y[18]/( k[49]* (1+y[8]/k[55] +y[4]/k[57] +y[1]/k[59]) +y[18] );
    dy[[18,17]] += [-rate, rate];
    
    # 12 XaVa activating II
    rate = k[50]*y[9]*y[18]/( k[51] +y[18] );
    dy[[18,17]] += [-rate, rate];
    
    # 13 IIa activating V
    rate = k[52]*y[17]*y[8]/( k[53]* (1+y[11]/k[31]+y[15]/k[61]) +y[8] );
    dy[[8,7]] += [-rate, rate];
    
    # 14 Xa activating V
    rate = k[54]*y[5]*y[8]/( k[55]* (1+y[18]/k[49]+y[4]/k[57]+y[1]/k[59]) +y[8] );
    dy[[8,7]] += [-rate, rate];
    
    # 15 Xa activating VII
    rate = k[56]*y[5]*y[4]/( k[57]* (1+y[8]/k[55]+y[18]/k[49]+y[1]/k[59]) +y[4] );
    dy[[4,3]] += [-rate, rate];
    
    # 16 Xa activating TFVII
    rate = k[58]*y[1]*y[5]/( k[59]* (1+y[8]/k[55]+y[18]/k[49]+y[4]/k[57]) +y[1] );
    dy[[1,2]] += [-rate, rate];
    
    # 17 IIa activating VIII
    rate = k[60]*y[15]*y[17]/( k[61]* (1+y[11]/k[31] +y[8]/k[53]) +y[15] );
    dy[[15,14]] += [-rate, rate];
    
    # 18 TMIIa activating PC
    rate = k[62]*y[20]*y[22]/( k[63] +y[22] );
    dy[[22,21]] += [-rate, rate];
    
    # 19 PC inhibiting Va
    rate = k[64]*y[21]*y[7]/( k[65] +y[14] +y[7] +y[9] +y[16] );
    dy[7] -= rate;
    
    # 20 PC inhibiting VIIIa
    rate = k[64]*y[21]*y[14]/( k[65] +y[14] +y[7] +y[9] +y[16] );
    dy[14] -= rate;
    
    # 21 PC inhibiting IXaVIIIa
    rate = k[64]*y[21]*y[16]/( k[65] +y[14] +y[7] +y[9] +y[16] );
    dy[[12,16]] += [rate, -rate];
    
    # 22 PC inhibiting XaVa
    rate = k[64]*y[21]*y[9]/( k[65] +y[14]+y[7]+y[9]+y[16] );
    dy[[9,5]] += [-rate, rate];
    
    
    # Mass Action Reactions
    # Va + Xa -> Xa==Va
    rate = k[0]*y[7]*y[5];
    dy[[5,7,9]] += [-rate, -rate, rate];
    
    # VIIIa + IXa -> VIIIa==IXa
    rate = k[1]*y[14]*y[12];
    dy[[12,14,16]] += [-rate, -rate, rate];
    
    # VIIa + TF -> TF==VIIa
    rate = k[2]*y[3]*y[0];
    dy[[3,0,2]] += [-rate, -rate, rate];
    
    # VII + TF -> TF==VII
    rate = k[3]*y[4]*y[0];
    dy[[4,0,1]] += [-rate, -rate, rate];
    
    # TF==VIIa + TFPI==Xa -> TF==VIIa==TFPI==Xa
    rate = k[4]*y[2]*y[24];
    dy[[2,24]] += [-rate, -rate];
    
    # TF==VIIa + ATIII -> TF==VIIa==ATIII
    rate = k[5]*y[2]*y[25];
    dy[[2,25]] += [-rate, -rate];
    
    # IIa + ATIII -> IIa==ATIII
    rate = k[6]*y[17]*y[25];
    dy[[17,25]] += [-rate, -rate];
    
    # IIa + a1AT -> IIa==a1AT
    rate = k[7]*y[17]*y[28];
    dy[[17,28]] += [-rate, -rate];
    
    # IIa + a2M -> IIa==a2M
    rate = k[8]*y[17]*y[30];
    dy[[17,30]] += [-rate, -rate];
    
    # IIa + PCI -> IIa==PCI
    rate = k[9]*y[17]*y[31];
    dy[[17,31]] += [-rate, -rate];
    
    # Xa + ATIII -> Xa==ATIII
    rate = k[10]*y[5]*y[25];
    dy[[5,25]] += [-rate, -rate];
    
    # Xa + a1AT -> Xa==a1AT
    rate = k[11]*y[5]*y[28];
    dy[[5,28]] += [-rate, -rate];
    
    # Xa + TFPI -> Xa==TFPI
    rate = k[12]*y[5]*y[23];
    dy[[5,23,24]] += [-rate, -rate, rate];
    
    # XaVa + a1AT -> Xa==a1AT + Va
    rate = k[13]*y[9]*y[28];
    dy[[9,28,7]] += [-rate, -rate, rate];
    
    # XaVa + ATIII -> Xa==ATIII + Va
    rate = k[14]*y[9]*y[25];
    dy[[9,25,7]] += [-rate, -rate, rate];
    
    # IXa + ATIII -> IXa==ATIII
    rate = k[15]*y[12]*y[25];
    dy[[12,25]] += [-rate, -rate];
    
    # VIIIa==IXa + ATIII -> IXa==ATIII + VIIIa
    rate = k[16]*y[16]*y[25];
    dy[[16,25,14]] += [-rate, -rate, rate];
    
    # XIa + C1-Inh -> XI==C1-Inh
    rate = k[17]*y[10]*y[26];
    dy[[10,26]] += [-rate, -rate];
    
    # XIa + a1AT -> XI==a1AT
    rate = k[18]*y[10]*y[28];
    dy[[10,28]] += [-rate, -rate];
    
    # XIa + ATIII -> XIa==ATIII
    rate = k[19]*y[10]*y[25];
    dy[[10,25]] += [-rate, -rate];
    
    # XIa + a2AP -> XIa==a2AP
    rate = k[20]*y[10]*y[29];
    dy[[10,29]] += [-rate, -rate];
    
    # XIa + PAI1 -> XIa==PAI1
    rate = k[21]*y[10]*y[27];
    dy[[10,27]] += [-rate, -rate];
    
    # IIa + TM -> IIa==TM
    rate = k[22]*y[17]*y[19];
    dy[[17,19,20]] += [-rate, -rate, rate];
    
    # IIa==TM + PCI -> IIa==TM==PCI
    rate = k[23]*y[20]*y[31];
    dy[[20,31]] += [-rate, -rate];
    
    # PCa + PCI -> PCa==PCI
    rate = k[24]*y[21]*y[31];
    dy[[21,31]] += [-rate, -rate];
    
    # PCa + a1AT -> PCa==a1AT
    rate = k[25]*y[21]*y[28];
    dy[[21,28]] += [-rate, -rate];
    
    # TF==VIIa -> TF + VIIa
    rate = k[26]*y[2];
    dy[[2,0,3]] += [-rate, rate, rate];
    
    # TF==VII -> TF + VII
    rate = k[27]*y[1];
    dy[[1,0,4]] += [-rate, rate, rate];
    
    return dy






