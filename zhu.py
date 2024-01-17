# 0  TF
# 1  TFVII
# 2  TFVIIa
# 3  VII
# 4  Xa
# 5  X
# 6  Va
# 7  V
# 8  XaVa
# 9  XIa
# 10 XI
# 11 IXa
# 12 IX
# 13 VIIIa
# 14 VIII
# 15 IXaVIIIa
# 16 IIa
# 17 II
# 18 TM
# 19 TMIIa
# 20 PCa
# 21 PC
# 22 TFPI
# 23 TFPIXa
# 24 ATIII
# 25 PAI1
# 26 alpha1AT
# 27 alpha2AP
# 28 alpha2M
# 29 XIIf
# 30 XIIa
# 31 XII
# 32 K
# 33 PK
# 34 C1
# 35 Ia
# 36 I

import numpy as np # The Numpy package is used frequently for this work
from scipy.integrate import solve_ivp # The solve_ivp method from Scipy is used to simulate the model
import matplotlib.pyplot as plt # matplotlib is used for plotting the thrombin generation curves in the plotThr method

def setIC(ICvector = np.array([10e-12, 1.4e-6, 2e-8, 1e-8, 1e-8/100, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.4e-6, 2.5e-9]), includeExtras = False):
    """
    Converts a list of factor levels into an initial condition vector for the Zhu model assuming 0 concentration for all remaining species. \n
    Inputs: \n
    ICvector - List of factor levels in the order TF, II, V, VII, VIIa, VIII, IX, X, XI, AT, TFPI. Leave black for default initial concentrations (10pM of TF), \n
    includeExtras - Boolean variable to determine whether or not the following inhibitors are included in the concentration vector (at their default concentrations): PCI, alpha1-AT, alpha2-M, alpha2-AP, PAI-1 and C1-inh. This variable should match the similar variable in getRates(). Default is False. \n
    Outputs: \n
    IC - Vector of initial conditions for the Zhu model
    """
    IC = np.zeros(37);
    IC[[0,17,7,3,2,14,12,5,10,24,22]] = ICvector; # Indices correspond to factors in initial condition vector
    if includeExtras:
        IC[[26, 28, 27, 25, 34]] = [4e-5, 3.25e-6, 9.5e-7, 4.6e-10, 2.1e-6];
    return IC

def getRates(includeExtras = False):
    """
    Get the reaction rates for simulating the Zhu model. \n
    Inputs: \n
    includeExtras - Only relevant for Panteleev model. Boolean variable to determine whether or not the following inhibitors are included in the rate vector (at their default concentrations): PCI, alpha1-AT, alpha2-M, alpha2-AP, PAI-1 and C1-inh. This variable should match the similar variable in setIC(). Default is False. \n
    Outputs: \n
    List of reaction rates
    """
    k = np.array([1.98/60, 216/60, 2400/60, 342/60, 0.34/60, 0.034/60, 34/60, 225/60, 0.04/60, 1740/60, 2.25/60, 1700/60, 14/60, 2.6/60, 60/60, 5040/60, 10000*1e6/60, 10000*1e6/60, 0.35*1e6/60, 0.0047*1e6/60, 0.0293*1e6/60, 0.11*1e6/60, 0.0157*1e6/60, 960*1e6/60, 0.0294*1e6/60, 0.001*1e6/60, 0.004*1e6/60, 0.01*1e6/60, 0.03*1e6/60, 12.6*1e6/60, 0.22*1e6/60, 0.011*1e6/60, 0.005*1e6/60, 0.0013*1e6/60, 0.96*1e6/60, 0.185*1e6/60, 0.0091*1e6/60, 0.0032*1e6/60, 1*1e6/60, 0.29*1e6/60, 3.6*1e6/60, 0.0096*1e6/60, 2*1e6/60, 0.027*1e6/60, 640*1e6/60, 39.6/60, 103/60, 34/60, 1200*1e6/60, 1200*1e6/60, 1200*1e6/60, 1200*1e6/60, 402*1e6/60, 19.8/60, 0.0078/60, 11/1e6, 0.091/1e6, 37/1e6, 0.51/1e6, 0.5/1e6, 2/1e6, 0.5/1e6, 0.35/1e6, 2/1e6, 0.19/1e6, 0.058/1e6, 1/1e6, 0.0717/1e6, 0.0104/1e6, 0.02/1e6, 7.2/1e6, 0.0093/1e6, 0.38/1e6, 0.133/1e6, 7.7/1e6]);
    return k

def getThr(k,y,maxt):
    """
    Simulate the Zhu model and return the thrombin generation curve. \n
    Inputs: \n
    k - List of reaction rates, \n
    y - Initial condition vector from setIC, \n
    maxt - Simulation time specified in seconds (1200s = 20min) \n
    Outputs: \n
    Tuple containing a vector of time points and a vector of thrombin concentrations
    """
    sol = ODESolver(k,y,maxt);
    return (sol.t,sol.y[16])

def plotThr(k,y,maxt):
    """
    Simulate the Zhu model and plots the thrombin generation curve. \n
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
    plt.title("Zhu Model")
    plt.show()

def test():
    """
    Simulate the Zhu model and plots the thrombin generation curve for default initial concentrations (Requires no inputs). \n
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
    Simulate the Zhu model and returns the ODE solution object. \n
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
    Evaluate the ODE of the Zhu model. \n
    Inputs: \n
    t - Current simulation time (not needed for evaluation but required by solve_ivp), \n
    y - Current concentration vector, \n
    k - List of reaction rates \n
    Outputs: \n
    dy - Gradient of the concentration vector
    """
    dy = np.zeros(37);
    
    # Va + Xa -> Va==Xa
    rate = k[16]*y[6]*y[4];
    dy[[6,4,8]] += [-rate, -rate, rate];
    
    # VIIIa + IXa -> VIIIa==IXa
    rate = k[17]*y[13]*y[11];
    dy[[13,11,15]] += [-rate, -rate, rate];
    
    # IIa + ATIII -> IIa==ATIII
    rate = k[18]*y[16]*y[24];
    dy[[16,24]] += [-rate, -rate];
    
    # IIa + a1AT -> IIa==a1AT
    rate = k[19]*y[16]*y[26];
    dy[[16,26]] += [-rate, -rate];
    
    # IIa + a2M -> IIa==a2M
    rate = k[20]*y[16]*y[28];
    dy[[16,28]] += [-rate, -rate];
    
    # Xa + ATIII -> Xa==ATIII
    rate = k[21]*y[4]*y[24];
    dy[[4,24]] += [-rate, -rate];
    
    # Xa + a1AT -> Xa==a1AT
    rate = k[22]*y[4]*y[26];
    dy[[4,26]] += [-rate, -rate];
    
    # Xa + TFPI -> Xa==TFPI
    rate = k[23]*y[4]*y[22];
    dy[[4,22,23]] += [-rate, -rate, rate];
    
    # IXa + ATIII -> IXa==ATIII
    rate = k[24]*y[11]*y[24];
    dy[[11,24]] += [-rate, -rate];
    
    # XIa + C1 -> XIa==C1
    rate = k[25]*y[9]*y[34];
    dy[[9,34]] += [-rate, -rate];
    
    # XIa + a1AT -> XIa==a1AT
    rate = k[26]*y[9]*y[26];
    dy[[9,26]] += [-rate, -rate];
    
    # XIa + ATIII -> XIa==ATIII
    rate = k[27]*y[9]*y[24];
    dy[[9,24]] += [-rate, -rate];
    
    # XIa + a2AP -> XIa==a2AP
    rate = k[28]*y[9]*y[27];
    dy[[9,27]] += [-rate, -rate];
    
    # XIa + PAI1 -> XIa==PAI1
    rate = k[29]*y[9]*y[25];
    dy[[9,25]] += [-rate, -rate];
    
    # XIIa + C1 -> XIIa==C1
    rate = k[30]*y[30]*y[34];
    dy[[30,34]] += [-rate, -rate];
    
    # XIIa + a2AP -> XIIa==a2AP
    rate = k[31]*y[30]*y[27];
    dy[[30,27]] += [-rate, -rate];
    
    # XIIa + a2M -> XIIa==a2M
    rate = k[32]*y[30]*y[28];
    dy[[30,28]] += [-rate, -rate];
    
    # XIIa + ATIII -> XIIa==ATIII
    rate = k[33]*y[30]*y[24];
    dy[[30,24]] += [-rate, -rate];
    
    # XIIa + PAI1 -> XIIa==PAI1
    rate = k[34]*y[30]*y[25];
    dy[[30,25]] += [-rate, -rate];
    
    # XIIf + C1 -> XIIf==C1
    rate = k[35]*y[29]*y[34];
    dy[[29,34]] += [-rate, -rate];
    
    # XIIf + a2AP -> XIIf==a2AP
    rate = k[36]*y[29]*y[27];
    dy[[29,27]] += [-rate, -rate];
    
    # XIIf + ATIII -> XIIf==ATIII
    rate = k[37]*y[29]*y[24];
    dy[[29,24]] += [-rate, -rate];
    
    # K + C1 -> K==C1
    rate = k[38]*y[32]*y[34];
    dy[[32,34]] += [-rate, -rate];
    
    # K + a2M -> K==a2M
    rate = k[39]*y[32]*y[28];
    dy[[32,28]] += [-rate, -rate];
    
    # K + PAI1 -> K==PAI1
    rate = k[40]*y[32]*y[25];
    dy[[32,25]] += [-rate, -rate];
    
    # K + ATIII -> K==ATIII
    rate = k[41]*y[32]*y[24];
    dy[[32,24]] += [-rate, -rate];
    
    # VII + TF -> VII==TF
    rate = k[42]*y[3]*y[0];
    dy[[3,0,1]] += [-rate, -rate, rate];
    
    # VIIa==TF + ATIII -> VIIa==TF==ATIII
    rate = k[43]*y[2]*y[24];
    dy[[2,24]] += [-rate, -rate];
    
    # VIIa==TF + Xa==TFPI -> VIIa==TF==Xa==TFPI
    rate = k[44]*y[2]*y[23];
    dy[[2,23]] += [-rate, -rate];
    
    # PCa + Va -> PCa==Va
    rate = k[48]*y[20]*y[6];
    dy[[20,6]] += [-rate, -rate];
    
    # PCa + VIIIa -> PCa==VIIIa
    rate = k[49]*y[20]*y[13];
    dy[[20,13]] += [-rate, -rate];
    
    # PCa + Va==Xa -> PCa==Va==Xa
    rate = k[50]*y[20]*y[8];
    dy[[20,8]] += [-rate, -rate];
    
    # PCa + VIIIa==IXa -> PCa==VIIIa==IXa
    rate = k[51]*y[20]*y[15];
    dy[[20,15]] += [-rate, -rate];
    
    # TM + IIa -> IIa==TM
    rate = k[52]*y[18]*y[16];
    dy[[18,16,19]] += [-rate, -rate, rate];
    
    
    # Enzyme Reactions
    # XII by XIIa
    rate = k[0]*y[31]*y[30]/(k[55]+y[31]);
    dy[[31,30]] += [-rate, rate];
    
    # PK by XIIa
    rate = k[1]*y[33]*y[30]/(k[56]+y[33]);
    dy[[33,32]] += [-rate, rate];
    
    # PK by XIIf
    rate = k[2]*y[33]*y[29]/(k[57]+y[33]);
    dy[[33,32]] += [-rate, rate];
    
    # XII by K
    rate = k[3]*y[31]*y[32]/(k[58]+y[31]);
    dy[[31,30]] += [-rate, rate];
    
    # XIIa by K
    rate = k[4]*y[30]*y[32]/(k[59]+y[30]);
    dy[[30,29]] += [-rate, rate];
    
    # XI by XIIa
    rate = k[5]*y[10]*y[30]/(k[60]+y[10]);
    dy[[10,9]] += [-rate, rate];
    
    # XII by XIa
    rate = k[6]*y[31]*y[9]/(k[61]+y[31]);
    dy[[31,30]] += [-rate, rate];
    
    # IX by XIa
    rate = k[7]*y[12]*y[9]/(k[62]+y[12]);
    dy[[12,11]] += [-rate, rate];
    
    # X by IXa
    rate = k[8]*y[5]*y[11]/(k[63]+y[5]);
    dy[[5,4]] += [-rate, rate];
    
    # X by IXa==VIIIa
    rate = k[9]*y[5]*y[15]/(k[64]+y[5]);
    dy[[5,4]] += [-rate, rate];
    
    # II by Xa
    rate = k[10]*y[17]*y[4]/(k[65]+y[17]);
    dy[[17,16]] += [-rate, rate];
    
    # II by Xa==Va
    rate = k[11]*y[17]*y[8]/(k[66]+y[17]);
    dy[[17,16]] += [-rate, rate];
    
    # V by IIa
    rate = k[12]*y[7]*y[16]/(k[67]+y[7]);
    dy[[7,6]] += [-rate, rate];
    
    # V by Xa
    rate = k[13]*y[7]*y[4]/(k[68]+y[7]);
    dy[[7,6]] += [-rate, rate];
    
    # VIII by IIa
    rate = k[14]*y[14]*y[16]/(k[69]+y[14]);
    dy[[14,13]] += [-rate, rate];
    
    # I by IIa
    rate = k[15]*y[36]*y[16]/(k[70]+y[36]);
    dy[[36,35]] += [-rate, rate];
    
    # VII==TF by Xa
    rate = k[45]*y[1]*y[4]/(k[71]+y[1]);
    dy[[1,2]] += [-rate, rate];
    
    # X by VIIa==TF
    rate = k[46]*y[5]*y[2]/(k[72]+y[5]);
    dy[[5,4]] += [-rate, rate];
    
    # IX by VIIa==TF
    rate = k[47]*y[12]*y[2]/(k[73]+y[12]);
    dy[[12,11]] += [-rate, rate];
    
    # PC by IIa==TM
    rate = k[53]*y[21]*y[19]/(k[74]+y[21]);
    dy[[21,20]] += [-rate, rate];
    
    # XI by IIa
    rate = k[54]*y[10];
    dy[[10,9]] += [-rate, rate];
    
    return dy








