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
# 34 Boc-VPR-MCA
# 35 Boc-VPR-MCA-IIa
# 36 Boc-VPR
# 37 AMC
# 38 XII
# 39 XIIa
# 40 XIIa==XII
# 41 PK
# 42 XIIa==PK
# 43 K
# 44 XII==K
# 45 Kinh
# 46 CTI
# 47 C1inh
# 48 XIIa==C1inh
# 49 XIIa==ATIII
# 50 XI
# 51 XI-IIa
# 52 XIa
# 53 XIIa==XI
# 54 XIa==ATIII
# 55 XIa==C1inh
# 56 a1AT
# 57 XIa==a1AT
# 58 a2AP
# 59 XIa==a2AP
# 60 XIa==IX
# 61 IXa==X
# 62 Xa==VIII
# 63 VIIa==IX
# 64 VIIa==X
# 65 Fbg
# 66 Fbg==IIa
# 67 Fbn1
# 68 FPA
# 69 Fbn2
# 70 Fbn12
# 71 Fbn22
# 72 FPB
# 73 Fbn2==IIa
# 74 Fbn12==IIa
# 75 Fbn12==IIa==ATIII
# 76 Fbn1==IIa
# 77 Fbn1==IIa==ATIII
# 78 Fbn2==IIa==ATIII
# 79 XIIa==CTI
# 80 eps

import numpy as np # The Numpy package is used frequently for this work
from scipy.integrate import solve_ivp # The solve_ivp method from Scipy is used to simulate the model
import matplotlib.pyplot as plt # matplotlib is used for plotting the thrombin generation curves in the plotThr method

def setIC(ICvector = np.array([10e-12, 1.4e-6, 2e-8, 1e-8, 1e-8/100, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.4e-6, 2.5e-9]), eps0=0.01, includeExtras = False):
    """
    Converts a list of factor levels into an initial condition vector for the Chatterjee model assuming 0 concentration for all remaining species. \n
    Inputs: \n
    ICvector - List of factor levels in the order TF, II, V, VII, VIIa, VIII, IX, X, XI, AT, TFPI. Leave black for default initial concentrations (10pM of TF), \n
    includeExtras - Boolean variable to determine whether or not the following inhibitors are included in the concentration vector (at their default concentrations): PCI, alpha1-AT, alpha2-M, alpha2-AP, PAI-1 and C1-inh. This variable should match the similar variable in getRates(). Default is False. \n
    Outputs: \n
    IC - Vector of initial conditions for the Chatterjee model
    """
    IC = np.zeros(81);
    IC[[0,13,20,1,3,14,10,7,50,28,25,80]] = np.append(ICvector,eps0); # np.append adds the initial value for eps onto the input. Indices correspond to factors in initial condition vector
    if includeExtras:
        IC[[56,58,47]] = [4e-5, 9.75e-7, 1.7e-6]; #a1AT, a2AP, C1-inh
    return IC

def getRates(includeExtras = False):
    """
    Get the reaction rates for simulating the Chatterjee model. \n
    Inputs: \n
    includeExtras - Only relevant for Panteleev model. Boolean variable to determine whether or not the following inhibitors are included in the rate vector (at their default concentrations): PCI, alpha1-AT, alpha2-M, alpha2-AP, PAI-1 and C1-inh. This variable should match the similar variable in setIC(). Default is False. \n    
    Outputs: \n
    List of reaction rates
    """
    k = np.array([3.1e-2, 3.2e6, 3.1e-5, 2.3e7, 4.4e5, 1.3e7, 2.3e4, 0.0105, 2.5e7, 6, 19, 2.2e7, 2.4, 1e7, 1.8, 7.5e3, 2.0e7, 1e-4, 1e7, 1e-5, 1e8, 8.2, 2.2e4, 6e-5, 1e-3, 2e7, 0.008, 4e8, 2.06, 1e8, 63.5, 1.5e7, 3.6e-4, 9e5, 0.011, 3.2e8, 5e7, 1.5e3, 7.1e3, 4.9e2, 7.1e3, 2.3e2, 1.0e8, 6100, 53.8, 5.0e-4, 1.0e8, 750, 0.033, 1e8, 3600, 40, 1e8, 45.3, 5.7, 27000, 0.011, 1e8, 2.4, 3600, 21.6, 1e8, 5, 1.3e-4, 1e8, 200, 5.7e-4, 3.19e6, 320, 1800, 100, 4300, 1e8, 41, 7.7, 1e8, 0.64, 7e-4, 1e8, 2.1, 0.023, 1e8, 0.9, 3.6e-5, 1e8, 210, 1.6e-6, 1e8, 636, 84, 1e8, 742.6, 7.4, 1e6, 0.064, 1e8, 701, 49, 1e8, 1000, 16000, 16000, 10000, 0.01, 0.005]);
    return k

def getThr(k,y,maxt):
    """
    Simulate the Chatterjee model and return the thrombin generation curve. \n
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
    Simulate the Chatterjee model and plots the thrombin generation curve. \n
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
    plt.title("Chatterjee Model")
    plt.show()

def test():
    """
    Simulate the Chatterjee model and plots the thrombin generation curve for default initial concentrations (Requires no inputs). \n
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
    Simulate the Chatterjee model and returns the ODE solution object. \n
    Inputs: \n
    k - List of reaction rates, \n
    y - Initial condition vector from setIC, \n
    maxt - Simulation time specified in seconds (1200s = 20min) \n
    Outputs: \n
    sol - ODE solution object containing timecourse data of all species
    """
    global IIaMax
    IIaMax = 0;
    sol = solve_ivp(ODE, (0,maxt), c0, args=(k, ), method='BDF', rtol=1e-4,atol=1e-16);
    return sol


def ODE(t,y,k):
    """
    Evaluate the ODE of the Chatterjee model. \n
    Inputs: \n
    t - Current simulation time (not needed for evaluation but required by solve_ivp), \n
    y - Current concentration vector, \n
    k - List of reaction rates \n
    Outputs: \n
    dy - Gradient of the concentration vector
    """
    
    global IIaMax
    IIaMax = max(IIaMax,y[7]);
    fIIaMax = IIaMax**1.6123/(IIaMax**1.6123+(2.4279e-9)**1.6123);
    epsMax = k[103]+(1-k[103])*fIIaMax;
    dy = np.zeros(81); 
    dy[80] = k[104]*(epsMax-y[80]);
    
    
    # TF + VII <- TF==VII
    rate = k[0]*y[2];
    dy[[2,0,1]] += [-rate, rate, rate];
    
    # TF + VII -> TF==VII
    rate = k[1]*y[0]*y[1];
    dy[[0,1,2]] += [-rate, -rate, rate];
    
    # TF + VIIa <- TF==VIIa
    rate = k[2]*y[4];
    dy[[4,0,3]] += [-rate, rate, rate];
    
    # TF + VIIa -> TF==VIIa
    rate = k[3]*y[0]*y[3];
    dy[[0,3,4]] += [-rate, -rate, rate];
    
    # TF==VIIa + VII -> TF==VIIa + VIIa
    rate = k[4]*y[4]*y[1];
    dy[[1,3]] += [-rate, rate];
    
    # Xa + VII -> Xa + VIIa
    rate = k[5]*y[5]*y[1];
    dy[[1,3]] += [-rate, rate];
    
    # IIa + VII -> IIa  + VIIa
    rate = k[6]*y[6]*y[1];
    dy[[1,3]] += [-rate, rate];
    
    # TF==VIIa + X <- TF==VIIa==X
    rate = k[7]*y[8];
    dy[[8,4,7]] += [-rate, rate, rate];
    
    # TF==VIIa + X -> TF==VIIa==X
    rate = k[8]*y[7]*y[4];
    dy[[4,7,8]] += [-rate, -rate, rate];
    
    # TF==VIIa==X -> TF==VIIa==Xa
    rate = k[9]*y[8];
    dy[[8,9]] += [-rate, rate];
    
    # TF==VIIa + Xa <- TF==VIIa==Xa
    rate = k[10]*y[9];
    dy[[9,4,5]] += [-rate, rate, rate];
    
    # TF==VIIa + Xa -> TF==VIIa==Xa
    rate = k[11]*y[4]*y[5]; 
    dy[[4,5,9]] += [-rate, -rate, rate];
    
    # TF==VIIa + IX <- TF==VIIa==IX
    rate = k[12]*y[11];
    dy[[11,4,10]] += [-rate, rate, rate];
    
    # TF==VIIa + IX -> TF==VIIa==IX
    rate = k[13]*y[4]*y[10];
    dy[[4,10,11]] += [-rate, -rate, rate];
    
    # TF==VIIa==IX -> TF==VIIa + IXa
    rate = k[14]*y[11];
    dy[[11,4,12]] += [-rate, rate, rate];
    
    # Xa + II -> Xa + IIa
    rate = k[15]*y[5]*y[13];
    dy[[13,6]] += [-rate, rate];
    
    # IIa + VIII -> IIa + VIIIa
    rate = k[16]*y[6]*y[14];
    dy[[14,15]] += [-rate, rate];
    
    # VIIIa + IXa <- IXa==VIIIa
    rate = k[17]*y[16];
    dy[[16,12,15]] += [-rate, rate, rate];
    
    # VIIIa + IXa -> IXa==VIIIa
    rate = k[18]*y[15]*y[12];
    dy[[12,15,16]] += [-rate, -rate, rate];
    
    # IXa==VIIIa + X <- IXa==VIIIa==X
    rate = k[19]*y[17];
    dy[[17,7,16]] += [-rate, rate, rate];
    
    # IXa==VIIIa + X -> IXa==VIIIa==X
    rate = k[20]*y[7]*y[16]; 
    dy[[7,16,17]] += [-rate, -rate, rate];
    
    # IXa==VIIIa==X -> IXa==VIIIa + Xa
    rate = k[21]*y[17];
    dy[[17,5,16]] += [-rate, rate, rate];
    
    # VIIIa <- VIIIa1-L + VIIIa
    rate = k[22]*y[18]*y[19];
    dy[[18,19,15]] += [-rate, -rate, rate];
    
    # VIIIa -> VIIIa1-L + VIIIa2
    rate = k[23]*y[15];
    dy[[15,18,19]] += [-rate, rate, rate];
    
    # IXa==VIIIa==X -> VIIIa1-L + VIIIa2 + X + IXa
    rate = k[24]*y[17];
    dy[[17,7,12,18,19]] += [-rate, rate, rate, rate, rate];
    
    # IXa==VIIIa -> VIIIa1-L + VIIIa2 + IXa
    rate = k[24]*y[16]; 
    dy[[16,12,18,19]] += [-rate, rate, rate, rate];
    
    # IIa + V -> IIa + Va
    rate = k[25]*y[6]*y[20];
    dy[[20,21]] += [-rate, rate];
    
    # Xa + Va <- Xa==Va
    rate = k[26]*y[22];
    dy[[22,5,21]] += [-rate, rate, rate];
    
    # Xa + Va -> Xa==Va
    rate = k[27]*y[5]*y[21];
    dy[[5,21,22]] += [-rate, -rate, rate];
    
    # Xa==Va + II <- Xa==Va==II
    rate = k[28]*y[23];
    dy[[23,13,22]] += [-rate, rate, rate];
    
    # Xa==Va + II -> Xa==Va==II
    rate = k[29]*y[22]*y[13];
    dy[[13,22,23]] += [-rate, -rate, rate];
    
    # Xa==Va==II -> Xa==Va + mIIa
    rate = k[30]*y[23];
    dy[[23,22,24]] += [-rate, rate, rate];
    
    # mIIa + Xa==Va -> IIa + Xa==Va
    rate = k[31]*y[24]*y[22];
    dy[[24,6]] += [-rate, rate];
    
    # Xa + TFPI <- Xa==TFPI
    rate = k[32]*y[26];
    dy[[26,5,25]] += [-rate, rate, rate];
    
    # Xa + TFPI -> Xa==TFPI
    rate = k[33]*y[5]*y[25];
    dy[[5,25,26]] += [-rate, -rate, rate];
    
    # TF==VIIa==Xa + TFPI <- TF==VIIa==Xa==TFPI
    rate = k[34]*y[27];
    dy[[27,9,25]] += [-rate, rate, rate];
    
    # TF==VIIa==Xa + TFPI -> TF==VIIa==Xa==TFPI
    rate = k[35]*y[9]*y[25];
    dy[[9,25,27]] += [-rate, -rate, rate];
    
    # TF==VIIa + Xa==TFPI -> TF==VIIa==Xa==TFPI
    rate = k[36]*y[4]*y[26];
    dy[[4,26,27]] += [-rate, -rate, rate];
    
    # Xa + ATIII -> Xa==ATIII
    rate = k[37]*y[28]*y[5];
    dy[[5,28,29]] += [-rate, -rate, rate];
    
    # mIIa + ATIII -> mIIa==ATIII
    rate = k[38]*y[24]*y[28];
    dy[[24,28,30]] += [-rate, -rate, rate];
    
    # IXa + ATIII -> IXa==ATIII
    rate = k[39]*y[12]*y[28];
    dy[[12,28,31]] += [-rate, -rate, rate];
    
    # IIa + ATIII -> IIa==ATIII
    rate = k[40]*y[6]*y[28];
    dy[[6,28,32]] += [-rate, -rate, rate];
    
    # TF==VIIa + ATIII -> TF==VIIa==ATIII
    rate = k[41]*y[4]*y[28];
    dy[[4,28,33]] += [-rate, -rate, rate];
    
    
    
    
    # Boc-VPR-MCA + IIa -> Boc-VPR-MCA==IIa
    rate = k[42]*y[34]*y[6];
    dy[[6,34,35]] += [-rate, -rate, rate];
    
    # Boc-VPR-MCA + IIa <- Boc-VPR-MCA==IIa
    rate = k[43]*y[35];
    dy[[35,6,34]] += [-rate, rate, rate];
    
    # Boc-VPR-MCA==IIa -> Boc-VPR + AMC + IIa
    rate = k[44]*y[35];
    dy[[35,36,37,6]] += [-rate, rate, rate, rate];
    
    # XII -> XIIa
    rate = k[45]*y[38];
    dy[[38,39]] += [-rate, rate];
    
    # XIIa + XII -> XIIa==XII
    rate = k[46]*y[38]*y[39];
    dy[[38,39,40]] += [-rate, -rate, rate];
    
    # XIIa + XII <- XIIa==XII
    rate = k[47]*y[40]/y[80];
    dy[[40,38,39]] += [-rate, rate, rate];
    
    # XIIa==XII -> 2 XIIa
    rate = k[48]*y[40];
    dy[[40,39]] += [-rate, 2*rate];
    
    #XIIa + PK -> XIIa==PK
    rate = k[49]*y[39]*y[41];
    dy[[39,41,42]] += [-rate, -rate, rate];
    
    #XIIa + PK <- XIIa==PK
    rate = k[50]*y[42]/y[80];
    dy[[42,39,41]] += [-rate, rate, rate];
    
    # XIIa==PK -> XIIa + K
    rate = k[51]*y[42];
    dy[[42,39,43]] += [-rate, rate, rate];
    
    # XII + K -> XII==K
    rate = k[52]*y[38]*y[43];
    dy[[38,43,44]] += [-rate, -rate, rate];
    
    #XII + K <- XII==K
    rate = k[53]*y[44]/y[80];
    dy[[44,38,43]] += [-rate, rate, rate];
    
    #XII==K -> XIIa + K
    rate = k[54]*y[44];
    dy[[44,39,43]] += [-rate, rate, rate];
    
    # PK + K -> K + K
    rate = k[55]*y[41]*y[43];
    dy[[41,43]] += [-rate, rate];
    
    # K -> Kinh
    rate = k[56]*y[43];
    dy[[43,45]] += [-rate, rate];
    
    # XIIa + CTI -> XIIa==CTI
    rate = k[57]*y[39]*y[46];
    dy[[39,46,79]] += [-rate, -rate, rate];
    
    # XIIa + CTI <- XIIa==CTI
    rate = k[58]*y[79];
    dy[[79,39,46]] += [-rate, rate, rate];
    
    # XIIa + C1inh -> XIIa==C1inh
    rate = k[59]*y[39]*y[47];
    dy[[39,47,48]] += [-rate, -rate, rate];
    
    # XIIa + ATIII -> XIIa==ATIII
    rate = k[60]*y[39]*y[28];
    dy[[39,28,49]] += [-rate, -rate, rate];
    
    # XI + IIa -> XI==IIa
    rate = k[61]*y[50]*y[6];
    dy[[50,6,51]] += [-rate, -rate, rate];
    
    # XI + IIa <- XI==IIa
    rate = k[62]*y[51];
    dy[[51,6,50]] += [-rate, rate, rate];
    
    # XI==IIa -> XIa + IIa
    rate = k[63]*y[51];
    dy[[51,6,52]] += [-rate, rate, rate];
    
    # XIIa + XI -> XIIa==XI
    rate = k[64]*y[39]*y[50];
    dy[[39,50,53]] += [-rate, -rate, rate];
    
    # XIIa + XI <- XIIa==XI
    rate = k[65]*y[53]/y[80];
    dy[[53,39,50]] += [-rate, rate, rate];
    
    # XIIa==XI -> XIIa + XIa
    rate = k[66]*y[53];
    dy[[53,39,52]] += [-rate, rate, rate];
    
    # XIa + XI -> 2 XIa
    rate = k[67]*y[52]*y[50];
    dy[[50,52]] += [-rate, rate];
    
    #XIa + ATIII -> XIa==ATIII
    rate = k[68]*y[52]*y[28];
    dy[[52,28,54]] += [-rate, -rate, rate];
    
    # XIa + C1inh -> XIa==C1inh
    rate = k[69]*y[52]*y[47];
    dy[[52,47,55]] += [-rate, -rate, rate];
    
    # XIa + a1AT -> XIa==a1AT
    rate = k[70]*y[52]*y[56];
    dy[[52,56,57]] += [-rate, -rate, rate];
    
    #XIa + a2AP -> XIa==a2AP
    rate = k[71]*y[52]*y[58];
    dy[[52,58,59]] += [-rate, -rate, rate];
    
    # XIa + IX -> XIa==IX
    rate = k[72]*y[52]*y[10];
    dy[[52,10,60]] += [-rate, -rate, rate];
    
    # XIa + IX <- XIa==IX
    rate = k[73]*y[60]/y[80];
    dy[[60,52,10]] += [-rate, rate, rate];
    
    # XIa==IX -> XIa + IXa
    rate = k[74]*y[60];
    dy[[60,52,12]] += [-rate, rate, rate];
    
    # IXa + X -> IXa==X
    rate = k[75]*y[12]*y[7];
    dy[[12,7,61]] += [-rate, -rate, rate];
    
    # IXa + X <- IXa==X
    rate = k[76]*y[61]/y[80];
    dy[[61,7,12]] += [-rate, rate, rate];
    
    # IXa==X -> IXa + Xa
    rate = k[77]*y[61];
    dy[[61,12,5]] += [-rate, rate, rate];
    
    # Xa + VIII -> Xa==VIII
    rate = k[78]*y[5]*y[14];
    dy[[5,14,62]] += [-rate, -rate, rate];
    
    # Xa + VIII <- Xa==VIII
    rate = k[79]*y[62]/y[80];
    dy[[62,5,14]] += [-rate, rate, rate];
    
    # Xa==VIII -> Xa + VIIIa
    rate = k[80]*y[62];
    dy[[62,5,15]] += [-rate, rate, rate];
    
    # VIIa + IX -> VIIa==IX
    rate = k[81]*y[3]*y[10];
    dy[[3,10,63]] += [-rate, -rate, rate];
    
    # VIIa + IX <- VIIa==IX
    rate = k[82]*y[63];
    dy[[63,10,3]] += [-rate, rate, rate];
    
    # VIIa==IX -> VIIa + IXa
    rate = k[83]*y[63];
    dy[[63,3,12]] += [-rate, rate, rate];
    
    # VIIa + X -> VIIa==X
    rate = k[84]*y[3]*y[7];
    dy[[3,7,64]] += [-rate, -rate, rate];
    
    # VIIa + X <- VIIa==X
    rate = k[85]*y[64]/y[80];
    dy[[64,3,7]] += [-rate, rate, rate];
    
    # VIIa==X -> VIIa + Xa
    rate = k[86]*y[64];
    dy[[64,3,5]] += [-rate, rate, rate];
    
    # Fbg + IIa -> Fbg==IIa
    rate = k[87]*y[65]*y[6];
    dy[[65,6,66]] += [-rate, -rate, rate];
    
    # Fbg + IIa <- Fbg==IIa
    rate = k[88]*y[66];
    dy[[66,6,65]] += [-rate, rate, rate];
    
    # Fbg==IIa -> Fbn1 + IIa + FPA
    rate = k[89]*y[66];
    dy[[66,6,67,68]] += [-rate, rate, rate, rate];
    
    # Fbn1 + IIa -> Fbn1==IIa
    rate = k[90]*y[6]*y[67];
    dy[[67,6,76]] += [-rate, -rate, rate];
    
    # Fbn1 + IIa <- Fbn1==IIa
    rate = k[91]*y[76];
    dy[[76,6,67]] += [-rate, rate, rate];
    
    # Fbn1==IIa -> Fbn2 + IIa + FPB
    rate = k[92]*y[76];
    dy[[76,6,69,72]] += [-rate, rate, rate, rate];
    
    # 2 Fbn1 -> Fbn12
    rate = k[93]*y[67]**2;
    dy[[67,70]] += [-2*rate, rate];
    
    # 2 Fbn1 <- Fbn12
    rate = k[94]*y[70];
    dy[[70,67]] += [-rate, 2*rate];
    
    # Fbn12 + IIa -> Fbn12==IIa
    rate = k[95]*y[70]*y[6];
    dy[[70,6,74]] += [-rate, -rate, rate];
    
    # Fbn12 + IIa <- Fbn12==IIa
    rate = k[96]*y[74];
    dy[[74,6,70]] += [-rate, rate, rate];
    
    # Fbn12==IIa -> Fbn22 + IIa + FPB
    rate = k[97]*y[74];
    dy[[74,6,71,72]] += [-rate, rate, rate, rate];
    
    # Fbn2 + IIa -> Fbn2==IIa
    rate = k[98]*y[6]*y[69];
    dy[[69,6,73]] += [-rate, -rate, rate];
    
    # Fbn2 + IIa <- Fbn2==IIa
    rate = k[99]*y[73];
    dy[[73,6,69]] += [-rate, rate, rate];
    
    # Fbn12==IIa + ATIII -> Fbn12==IIa==ATIII
    rate = k[100]*y[74]*y[28];
    dy[[74,28,75]] += [-rate, -rate, rate];
    
    # Fbn1==IIa + ATIII -> Fbn1==IIa==ATIII
    rate = k[101]*y[76]*y[28];
    dy[[76,28,77]] += [-rate, -rate, rate];
    
    # Fbn2==IIa + ATIII -> Fbn2==IIa==ATIII
    rate = k[102]*y[73]*y[28];
    dy[[73,28,78]] += [-rate, -rate, rate];
    
    return dy





