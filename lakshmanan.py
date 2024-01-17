# 0  TF
# 1  VII
# 2  TF==VII
# 3  VIIa
# 4  TF==VIIa
# 5  Xa
# 6  IIa
# 7  X
# 8  TF==VIIa==X
# 9  TF==VIIa==Xa
# 10 IX
# 11 TF==VIIa==IX
# 12 IXa
# 13 II
# 14 VIII
# 15 VIIIa
# 16 IXa==VIIIa
# 17 IXa==VIIIa==X
# 18 VIIIa1-L
# 19 VIIIa2
# 20 V
# 21 Va
# 22 Xa==Va
# 23 Xa==Va==II
# 24 mIIa
# 25 TFPI
# 26 Xa==TFPI
# 27 TF==VIIa==Xa==TFPI
# 28 ATIII
# 29 Xa==ATIII
# 30 mIIa==ATIII
# 31 IXa==ATIII
# 32 IIa==ATIII
# 33 TF==VIIa==ATIII
# 34 XI
# 35 XI==IIa
# 36 XIa
# 37 XIa==IX
# 38 XIa==ATIII
# 39 C1-inh
# 40 XIa==C1-inh
# 41 IXa==X

import numpy as np # The Numpy package is used frequently for this work
from scipy.integrate import solve_ivp # The solve_ivp method from Scipy is used to simulate the model
import matplotlib.pyplot as plt # matplotlib is used for plotting the thrombin generation curves in the plotThr method

def setIC(ICvector = np.array([10e-12, 1.4e-6, 2e-8, 1e-8, 1e-8/100, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.4e-6, 2.5e-9]), includeExtras = False):
    """
    Converts a list of factor levels into an initial condition vector for the Lakshmanan model assuming 0 concentration for all remaining species. \n
    Inputs: \n
    ICvector - List of factor levels in the order TF, II, V, VII, VIIa, VIII, IX, X, XI, AT, TFPI. Leave black for default initial concentrations (10pM of TF), \n
    includeExtras - Boolean variable to determine whether or not the following inhibitors are included in the concentration vector (at their default concentrations): PCI, alpha1-AT, alpha2-M, alpha2-AP, PAI-1 and C1-inh. This variable should match the similar variable in getRates(). Default is False. \n
    Outputs: \n
    IC - Vector of initial conditions for the Lakshmanan model
    """
    IC = np.zeros(42);
    IC[[0,13,20,1,3,14,10,7,34,28,25]] = ICvector; # Indices correspond to factors in initial condition vector
    if includeExtras:
        IC[39] = 1.7e-6; #C1-inh
    return IC

def getRates(includeExtras = False):
    """
    Get the reaction rates for simulating the Lakshmanan model. \n
    Inputs: \n
    includeExtras - Only relevant for Panteleev model. Boolean variable to determine whether or not the following inhibitors are included in the rate vector (at their default concentrations): PCI, alpha1-AT, alpha2-M, alpha2-AP, PAI-1 and C1-inh. This variable should match the similar variable in setIC(). Default is False. \n    
    Outputs: \n
    List of reaction rates
    """
    k = np.array([1.7e5, 3e-3, 2.2e8, 3.1e-5, 4.4e5, 4.4e5, 2.5e9, 2.5e9, 2.3e4, 2.3e4, 7.5e6, 5.2e-1, 1.1e1, 2.2e7, 3.9e1, 2.7e7, 4.2e2, 1e1, 9.2e3, 2.5e7, 5e7, 2.8e-3, 1.3e8, 1e-3, 4.2e1, 7.5e-3, 2.1e4, 1.4e-4, 6.1e-4, 2.3e7, 4.9e8, 3e-1, 2.5e7, 7.8e1, 1.3e1, 2.3e7, 2.2e7, 3.7e-5, 1e7, 1e-5, 4.4e8, 1.2e3, 1e4, 7e2, 2.1e3, 3.3e2, 5e7, 9.9, 1.1e-4, 6.1e5, 9.9e-1, 1.1e-1, 3.2e2, 1.8e3, 1e8, 3.3e2, 3e-3]);
    return k

def getThr(k,y,maxt):
    """
    Simulate the Lakshmanan model and return the thrombin generation curve. \n
    Inputs: \n
    k - List of reaction rates, \n
    y - Initial condition vector from setIC, \n
    maxt - Simulation time specified in seconds (1200s = 20min) \n
    Outputs: \n
    Tuple containing a vector of time points and a vector of thrombin concentrations
    """
    sol = ODESolver(k,y,maxt);
    return (sol.t,sol.y[6])

def plotThr(k,y,maxt):
    """
    Simulate the Lakshmanan model and plots the thrombin generation curve. \n
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
    plt.title("Lakshmanan Model")
    plt.show()

def test():
    """
    Simulate the Lakshmanan model and plots the thrombin generation curve for default initial concentrations (Requires no inputs). \n
    Inputs: \n
    None \n
    Outputs: \n
    None, plot of thrombin generation curve is displayed.
    """
    y = setIC();
    k = getRates();
    plotThr(k,y,1200)

def ODESolver(k,y,maxt):
    """
    Simulate the Lakshmanan model and returns the ODE solution object. \n
    Inputs: \n
    k - List of reaction rates, \n
    y - Initial condition vector from setIC, \n
    maxt - Simulation time specified in seconds (1200s = 20min) \n
    Outputs: \n
    sol - ODE solution object containing timecourse data of all species
    """
    sol = solve_ivp(ODE, (0,maxt), y, args=(k, ), method='BDF', rtol=1e-4,atol=1e-16);
    return sol

def ODE(t,y,k):
    """
    Evaluate the ODE of the Lakshmanan model. \n
    Inputs: \n
    t - Current simulation time (not needed for evaluation but required by solve_ivp), \n
    y - Current concentration vector, \n
    k - List of reaction rates \n
    Outputs: \n
    dy - Gradient of the concentration vector
    """
    
    dy = np.zeros(42);

    # TF + VII <- TF==VII
    rate = k[1]*y[2];
    dy[[2,0,1]] += [-rate, rate, rate];

    # TF + VII -> TF==VII
    rate = k[0]*y[0]*y[1];
    dy[[0,1,2]] += [-rate, -rate, rate];

    # TF + VIIa <- TF==VIIa
    rate = k[3]*y[4];
    dy[[4,0,3]] += [-rate, rate, rate];

    # TF + VIIa -> TF==VIIa
    rate = k[2]*y[0]*y[3];
    dy[[0,3,4]] += [-rate, -rate, rate];

    # TF==VIIa + VII -> TF==VIIa + VIIa
    rate = k[4]*y[4]*y[1];
    dy[[1,3]] += [-rate, rate];

    # TF==VIIa + TF==VII -> TF==VIIa + TF==VIIa
    rate = k[5]*y[4]*y[2];
    dy[[2,4]] += [-rate, rate];

    # Xa + VII -> Xa + VIIa
    rate = k[6]*y[5]*y[1];
    dy[[1,3]] += [-rate, rate];

    # Xa + TF==VII -> Xa + TF==VIIa
    rate = k[7]*y[5]*y[2];
    dy[[2,4]] += [-rate, rate];

    # IIa + VII -> IIa  + VIIa
    rate = k[8]*y[6]*y[1];
    dy[[1,3]] += [-rate, rate];

    # IIa + TF==VII -> IIa  + TF==VIIa
    rate = k[9]*y[6]*y[2];
    dy[[2,4]] += [-rate, rate];

    # TF==VIIa + X <- TF==VIIa==X
    rate = k[11]*y[8];
    dy[[8,4,7]] += [-rate, rate, rate];

    # TF==VIIa + X -> TF==VIIa==X
    rate = k[10]*y[7]*y[4];
    dy[[4,7,8]] += [-rate, -rate, rate];

    # TF==VIIa==X -> TF==VIIa==Xa
    rate = k[12]*y[8];
    dy[[8,9]] += [-rate, rate];

    # TF==VIIa + Xa <- TF==VIIa==Xa
    rate = k[14]*y[9];
    dy[[9,4,5]] += [-rate, rate, rate];

    # TF==VIIa + Xa -> TF==VIIa==Xa
    rate = k[13]*y[4]*y[5];
    dy[[4,5,9]] += [-rate, -rate, rate];

    # TF==VIIa + IX <- TF==VIIa==IX
    rate = k[16]*y[11];
    dy[[11,4,10]] += [-rate, rate, rate];

    # TF==VIIa + IX -> TF==VIIa==IX
    rate = k[15]*y[4]*y[10];
    #continue from here
    dy[[4,10,11]] += [-rate, -rate, rate];

    # TF==VIIa==IX -> TF==VIIa + IXa
    rate = k[17]*y[11];
    dy[[11,4,12]] += [-rate, rate, rate];

    # Xa + II -> Xa + IIa
    rate = k[18]*y[5]*y[13];
    dy[[13,6]] += [-rate, rate];

    # IIa + VIII -> IIa + VIIIa
    rate = k[19]*y[6]*y[14];
    dy[[14,15]] += [-rate, rate];

    # VIIIa + IXa <- IXa==VIIIa
    rate = k[21]*y[16];
    dy[[16,12,15]] += [-rate, rate, rate];

    # VIIIa + IXa -> IXa==VIIIa
    rate = k[20]*y[15]*y[12];
    dy[[12,15,16]] += [-rate, -rate, rate];

    # IXa==VIIIa + X <- IXa==VIIIa==X
    rate = k[23]*y[17];
    dy[[17,7,16]] += [-rate, rate, rate];

    # IXa==VIIIa + X -> IXa==VIIIa==X
    rate = k[22]*y[7]*y[16];
    dy[[7,16,17]] += [-rate, -rate, rate];

    # IXa==VIIIa==X -> IXa==VIIIa + Xa
    rate = k[24]*y[17];
    dy[[17,5,16]] += [-rate, rate, rate];

    # VIIIa <- VIIIa1-L + VIIIa2
    rate = k[26]*y[18]*y[19];
    dy[[18,19,15]] += [-rate, -rate, rate];

    # VIIIa -> VIIIa1-L + VIIIa2
    rate = k[25]*y[15];
    dy[[15,18,19]] += [-rate, rate, rate];

    # IXa==VIIIa==X -> VIIIa1-L + VIIIa2 + X + IXa
    rate = k[27]*y[17];
    dy[[17,7,12,18,19]] += [-rate, rate, rate, rate, rate];

    # IXa==VIIIa -> VIIIa1-L + VIIIa2 + IXa
    rate = k[28]*y[16];
    dy[[16,12,18,19]] += [-rate, rate, rate, rate];

    # IIa + V -> IIa + Va
    rate = k[29]*y[6]*y[20];
    dy[[20,21]] += [-rate, rate];

    # Xa + Va <- Xa==Va
    rate = k[31]*y[22];
    dy[[22,5,21]] += [-rate, rate, rate];

    # Xa + Va -> Xa==Va
    rate = k[30]*y[5]*y[21];
    dy[[5,21,22]] += [-rate, -rate, rate];

    # Xa==Va + II <- Xa==Va==II
    rate = k[33]*y[23];
    dy[[23,13,22]] += [-rate, rate, rate];

    # Xa==Va + II -> Xa==Va==II
    rate = k[32]*y[22]*y[13];
    dy[[13,22,23]] += [-rate, -rate, rate];

    # Xa==Va==II -> Xa==Va + mIIa
    rate = k[34]*y[23];
    dy[[23,22,24]] += [-rate, rate, rate];

    # mIIa + Xa==Va -> IIa + Xa==Va
    rate = k[35]*y[24]*y[22];
    dy[[24,6]] += [-rate, rate];

    # Xa + TFPI <- Xa==TFPI
    rate = k[37]*y[26];
    dy[[26,5,25]] += [-rate, rate, rate];

    # Xa + TFPI -> Xa==TFPI
    rate = k[36]*y[5]*y[25];
    dy[[5,25,26]] += [-rate, -rate, rate];

    # TF==VIIa==Xa + TFPI <- TF==VIIa==Xa==TFPI
    rate = k[39]*y[27];
    dy[[27,9,25]] += [-rate, rate, rate];

    # TF==VIIa==Xa + TFPI -> TF==VIIa==Xa==TFPI
    rate = k[38]*y[9]*y[25];
    dy[[9,25,27]] += [-rate, -rate, rate];

    # TF==VIIa + Xa==TFPI -> TF==VIIa==Xa==TFPI
    rate = k[40]*y[4]*y[26];
    dy[[4,26,27]] += [-rate, -rate, rate];

    # Xa + ATIII -> Xa==ATIII
    rate = k[41]*y[28]*y[5];
    dy[[5,28,29]] += [-rate, -rate, rate];

    # mIIa + ATIII -> mIIa==ATIII
    rate = k[42]*y[24]*y[28];
    dy[[24,28,30]] += [-rate, -rate, rate];

    # IXa + ATIII -> IXa==ATIII
    rate = k[43]*y[12]*y[28];
    dy[[12,28,31]] += [-rate, -rate, rate];

    # IIa + ATIII -> IIa==ATIII
    rate = k[44]*y[6]*y[28];
    dy[[6,28,32]] += [-rate, -rate, rate];
    
    # TF==VIIa + ATIII -> TF==VIIa==ATIII
    rate = k[45]*y[4]*y[28];
    dy[[4,28,33]] += [-rate, -rate, rate];

    # XI + IIa -> XI==IIa
    rate = k[46]*y[34]*y[6];
    dy[[34,6,35]] += [-rate, -rate, rate];

    # XI + IIa <- XI==IIa
    rate = k[47]*y[35];
    dy[[35,34,6]] += [-rate, rate, rate];
    
    # XI==IIa -> XIa + IIa
    rate = k[48]*y[35];
    dy[[35,36,6]] += [-rate, rate, rate];
    
    # XIa + IX -> XIa==IX
    rate = k[49]*y[36]*y[10];
    dy[[36,10,37]] += [-rate, -rate, rate];
    
    # XIa + IX <- XIa==IX
    rate = k[50]*y[37];
    dy[[37,36,10]] += [-rate, rate, rate];
    
    # XIa==IX -> XIa + IXa
    rate = k[51]*y[37];
    dy[[37,36,12]] += [-rate, rate, rate];
    
    # XIa + ATIII -> XIa==ATIII
    rate = k[52]*y[36]*y[28];
    dy[[36,28,38]] += [-rate, -rate, rate];
    
    # XIa + C1-inh -> XIa==C1-inh
    rate = k[53]*y[36]*y[39];
    dy[[36,39,40]] += [-rate, -rate, rate];
    
    # IXa + X -> IXa==X
    rate = k[54]*y[12]*y[7];
    dy[[12,7,41]] += [-rate, -rate, rate];
    
    # IXa + X <- IXa==X
    rate = k[55]*y[41];
    dy[[41,12,7]] += [-rate, rate, rate];
    
    # IXa==X -> IXa + Xa
    rate = k[56]*y[41];
    dy[[41,12,5]] += [-rate, rate, rate];
    
    return dy




