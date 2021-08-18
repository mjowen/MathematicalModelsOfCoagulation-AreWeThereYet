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
# 34 TM
# 35 TM==IIa
# 36 PC
# 37 TM==IIa==PC
# 38 APC
# 39 APC==Va
# 40 Va5
# 41 Va3
# 42 APC==Va5
# 43 APC==Va3
# 44 Va53
# 45 HCF
# 46 LCA1
# 47 APC==LCA1
# 48 TM==IIa==APC
# 49 Xa==Va5
# 50 Xa==Va3
# 51 Xa==Va5==II
# 52 Xa==Va3==II
# 53 TM==mIIa
# 54 TM==mIIa==PC
# 55 Xa==Va53
# 56 Xa==Va53==II
# 57 II==Va

import numpy as np # The Numpy package is used frequently for this work
from scipy.integrate import solve_ivp # The solve_ivp method from Scipy is used to simulate the model
import matplotlib.pyplot as plt # matplotlib is used for plotting the thrombin generation curves in the plotThr method

def setIC(ICvector = np.array([10e-12, 1.45e-6, 2e-8, 1e-8, 1e-8/100, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.45e-6, 2.5e-9])):
    """
    Converts a list of factor levels into an initial condition vector for the Brummel model assuming 0 concnetration for all remaining species. \n
    Inputs: \n
    ICvector - List of factor levels in the order TF, II, V, VII, VIIa, VIII, IX, X, XI, AT, TFPI. Leave black for default initial concentrations (10pM of TF) \n
    Outputs: \n
    IC - Vector of initial conditions for the Brummel model
    """
    IC = np.zeros(58);
    IC[[0,13,20,1,3,14,10,7,28,25]] = np.delete(ICvector,8); # np.delete removes XI from the input. Indices correspond to factors in initial condition vector
    return IC

def getRates():
    """
    Get the reaction rates for simulating the Brummel model. \n
    Inputs: \n
    None \n
    Outputs: \n
    List of reaction rates
    """
    k = np.array([3.1e-3, 3.2e6, 3.1e-3, 2.3e7, 4.4e5, 1.3e7, 2.3e4, 1.05, 2.5e7, 6, 19, 2.2e7, 2.4, 1e7, 1.8, 7.5e3, 2.0e7, 5e-3, 1e7, 1e-3, 1e8, 8.2, 2.2e4, 6e-3, 1e-3, 2e7, 0.075, 1.5e8, 103, 1e8, 63.5, 2.3e8, 3.6e-4, 9e5, 1.1e-4, 3.2e8, 5e7, 4.2e3, 7.1e3, 4.9e2, 7.1e3, 2.3e2, 0.33, 1e8, 100, 1e8, 0.41, 7.1e3, 0.7, 1e8, 1, 0.192, 0.028, 0.15, 10.3, 10.3, 4.6e7, 4.6e7, 0.0035, 5.7e3, 3e6, 70, 1e8, 4.05e6]);
    return k

def getThr(k,y,maxt):
    """
    Simulate the Brummel model and return the thrombin generation curve. \n
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
    Simulate the Brummel model and plots the thrombin generation curve. \n
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
    plt.title("Brummel Model")
    plt.show()

def test():
    """
    Simulate the Brummel model and plots the thrombin generation curve for default initial concentrations (Requires no inputs). \n
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
    Simulate the Brummel model and returns the ODE solution object. \n
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
    Evaluate the ODE of the Brummel model. \n
    Inputs: \n
    t - Current simulation time (not needed for evaluation but required by solve_ivp), \n
    y - Current concentration vector, \n
    k - List of reaction rates \n
    Outputs: \n
    dy - Gradient of the concentration vector
    """
    
    dy = np.zeros(58); 
    
    
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





    # TM + IIa <- TM==IIa
    rate = k[42]*y[35];
    dy[[35,34,6]] += [-rate, rate, rate];
    
    # TM + IIa -> TM==IIa
    rate = k[43]*y[34]*y[6];
    dy[[34,6,35]] += [-rate, -rate, rate];
    
    # TM==IIa + PC <- TM==IIa==PC
    rate = k[44]*y[37];
    dy[[37,35,36]] += [-rate, rate, rate];
    
    # TM==IIa + PC -> TM==IIa==PC
    rate = k[45]*y[35]*y[36];
    dy[[35,36,37]] += [-rate, -rate, rate];
    
    # TM==IIa==PC -> TM==IIa + APC
    rate = k[46]*y[37];
    dy[[37,35,38]] += [-rate, rate, rate];
    
    # TM==IIa + ATIII -> IIa==ATIII + TM
    rate = k[47]*y[35]*y[28];
    dy[[28,35,32,34]] += [-rate, -rate, rate, rate];
    
    # APC + Va <- APC==Va
    rate = k[48]*y[39];
    dy[[39,38,21]] += [-rate, rate, rate];
    
    # APC + Va -> APC==Va
    rate = k[49]*y[21]*y[38];
    dy[[21,38,39]] += [-rate, -rate, rate];
    
    # APC==Va -> APC + Va5
    rate = k[50]*y[39];
    dy[[39,38,40]] += [-rate, rate, rate];
    
    # APC==Va -> APC + Va3
    rate = k[51]*y[39];
    dy[[39,38,41]] += [-rate, rate, rate];
    
    # APC + Va5 <- APC==Va5
    rate = k[48]*y[42];
    dy[[42,38,40]] += [-rate, rate, rate];
    
    # APC + Va5 -> APC==Va5
    rate = k[49]*y[38]*y[40];
    dy[[38,40,42]] += [-rate, -rate, rate];
    
    # APC + Va3 <- APC==Va3
    rate = k[48]*y[43];
    dy[[43,38,41]] += [-rate, rate, rate];
    
    # APC + Va3 -> APC==Va3
    rate = k[49]*y[38]*y[41];
    dy[[38,41,43]] += [-rate, -rate, rate];
    
    # APC==Va3 -> APC + Va53
    rate = k[50]*y[43];
    dy[[43,38,44]] += [-rate, rate, rate];
    
    # APC==Va5 -> APC + Va53
    rate = k[51]*y[42];
    dy[[42,38,44]] += [-rate, rate, rate];
    
    # Va3 -> HCF + LCA1
    rate = k[52]*y[41];
    dy[[41,45,46]] += [-rate, rate, rate];
    
    # Va53 -> HCF + LCA1
    rate = k[52]*y[44];
    dy[[44,45,46]] += [-rate, rate, rate];
    
    # APC + LCA1 <- APC==LCA1
    rate = k[48]*y[47];
    dy[[47,38,46]] += [-rate, rate, rate];
    
    # APC + LCA1 -> APC==LCA1
    rate = k[49]*y[38]*y[46];
    dy[[38,46,47]] += [-rate, -rate, rate];
    
    # APC + TM==IIa <- TM==IIa==APC
    rate = k[44]*y[48];
    dy[[48,38,35]] += [-rate, rate, rate];
    
    # APC + TM==IIa -> TM==IIa==APC
    rate = k[45]*y[38]*y[35];
    dy[[35,38,48]] += [-rate, -rate, rate];
    
    # Xa + Va5 <- Xa==Va5
    rate = k[53]*y[49];
    dy[[49,5,40]] += [-rate, rate, rate];
    
    # Xa + Va5 -> Xa==Va5
    rate = k[27]*y[5]*y[40];
    dy[[5,40,49]] += [-rate, -rate, rate];
    
    # Xa + Va3 <- Xa==Va3
    rate = k[53]*y[50];
    dy[[50,5,41]] += [-rate, rate, rate];
    
    # Xa + Va3 -> Xa==Va3
    rate = k[27]*y[5]*y[41];
    dy[[5,41,50]] += [-rate, -rate, rate];
    
    # Xa==Va5 + II <- Xa==Va5==II
    rate = k[28]*y[51];
    dy[[51,49,13]] += [-rate, rate, rate];
    
    # Xa==Va5 + II -> Xa==Va5==II
    rate = k[29]*y[49]*y[13];
    dy[[49,13,51]] += [-rate, -rate, rate];
    
    # Xa==Va5==II -> Xa==Va5 + mIIa
    rate = k[54]*y[51];
    dy[[51,49,24]] += [-rate, rate, rate];
    
    # Xa==Va3 + II <- Xa==Va3==II
    rate = k[28]*y[52];
    dy[[52,13,50]] += [-rate, rate, rate];
    
    # Xa==Va3 + II -> Xa==Va3==II
    rate = k[29]*y[13]*y[50];
    dy[[13,50,52]] += [-rate, -rate, rate];
    
    # Xa==Va3==II -> Xa==Va3 + mIIa
    rate = k[55]*y[52];
    dy[[52,50,24]] += [-rate, rate, rate];
    
    # Xa==Va5 + mIIa -> IIa + Xa==Va5
    rate = k[56]*y[49]*y[24];
    dy[[24,6]] += [-rate, rate];
    
    # Xa==Va3 + mIIa -> IIa + Xa==Va3
    rate = k[57]*y[50]*y[24];
    dy[[24,6]] += [-rate, rate];
    
    # Xa==Va3 -> HCF + LCA1 + Xa
    rate = k[58]*y[50];
    dy[[50,45,46,5]] += [-rate, rate, rate, rate];
    
    # Xa==Va3==II -> HCF + LCA1 + Xa + II
    rate = k[58]*y[52];
    dy[[52,45,46,5,13]] += [-rate, rate, rate, rate, rate];
    
    # IXa + X -> IXa + Xa
    rate = k[59]*y[12]*y[7];
    dy[[7,5]] += [-rate, rate];
    
    # mIIa + V -> mIIa + Va
    rate = k[60]*y[24]*y[20];
    dy[[20,21]] += [-rate, rate];
    
    # TM + mIIa <- TM==mIIa
    rate = k[42]*y[53];
    dy[[53,34,24]] += [-rate, rate, rate];
    
    # TM + mIIa -> TM==mIIa
    rate = k[43]*y[34]*y[24];
    dy[[24,34,53]] += [-rate, -rate, rate];
    
    # TM==mIIa + PC <- TM==mIIa==PC
    rate = k[44]*y[54];
    dy[[54,53,36]] += [-rate, rate, rate];
    
    # TM==mIIa + PC -> TM==mIIa==PC
    rate = k[45]*y[53]*y[36];
    dy[[53,36,54]] += [-rate, -rate, rate];
    
    # TM==mIIa==PC -> TM==mIIa + APC
    rate = k[46]*y[54];
    dy[[54,53,38]] += [-rate, rate, rate];
    
    # TM==mIIa + ATIII -> mIIa==ATIII + TM
    rate = k[47]*y[53]*y[28];
    dy[[53,28,30,34]] += [-rate, -rate, rate, rate];
    
    # Xa + Va53 <- Xa==Va53
    rate = k[53]*y[55];
    dy[[55,5,44]] += [-rate, rate, rate];
    
    # Xa + Va53 -> Xa==Va53
    rate = k[27]*y[5]*y[44];
    dy[[5,44,55]] += [-rate, -rate, rate];
    
    # Xa==Va53 + II <- Xa==Va53==II
    rate = k[28]*y[56];
    dy[[56,55,13]] += [-rate, rate, rate];
    
    # Xa==Va53 + II -> Xa==Va53==II
    rate = k[29]*y[13]*y[55];
    dy[[13,55,56]] += [-rate, -rate, rate];
    
    # Xa==Va53==II -> Xa==Va53 + mIIa
    rate = k[55]*y[56];
    dy[[56,55,24]] += [-rate, rate, rate];
    
    # Xa==Va53 + mIIa -> Xa==Va53 + IIa
    rate = k[57]*y[55]*y[24];
    dy[[24,6]] += [-rate, rate];
    
    # Xa==Va53 -> HCF + LCA1 + Xa
    rate = k[58]*y[55];
    dy[[55,45,46,5]] += [-rate, rate, rate, rate];
    
    # Xa==Va53==II -> HCF + LCA1 + Xa + II
    rate = k[58]*y[56];
    dy[[56,45,46,5,13]] += [-rate, rate, rate, rate, rate];
    
    # II + Va <- II==Va
    rate = k[61]*y[57];
    dy[[57,13,21]] += [-rate, rate, rate];
    
    # II + Va -> II==Va
    rate = k[62]*y[13]*y[21];
    dy[[13,21,57]] += [-rate, -rate, rate];
    
    # Xa==Va5 + APC -> Xa==Va53 + APC
    rate = k[63]*y[49]*y[38];
    dy[[49,55]] += [-rate, rate];
    
    return dy


