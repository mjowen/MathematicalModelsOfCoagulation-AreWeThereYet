#### IMPORT ####
import numpy as np # The Numpy package is used frequently for this work
from scipy.integrate import trapz # The trapezium fuction from the Scipy package is used for the sumStat function
import matplotlib.pyplot as plt # The matplotlib package is used for the ETP scatter plots in the ETP Correlation method

# Import the 8 models as individual module
import hockin
import danforth
import chatterjee
import brummel
import bungay
import panteleev
import tyurin
import lakshmanan

#### Function Definitions ####
def sumStat(t,thr):
    """
    Calulate the 6 summary statistics used for the sensitivity analysis. \n
    Inputs: \n
    t - List of time points, \n
    thr - List of thrombin concentrations (in moles) at each time point. \n
    Outputs: \n
    List of summary statistics: \n
    ETP (Endogenous Thrombin Potential), Peak Concentration, Time to Peak Concentration,
    Lagtime (Time to 1% of peak concentration), Maximum Increasing Rate and 
    Minimum Decreasing Rate.
    """
    # ETP is found using scipy's trapezium method to integrate the thrombin generation curve
    etp = trapz(thr,t);
    peak = np.max(thr);
    # The list index for the peak concentration is found using numpy's argmax and the corresponding timepoint is extracted
    peakIdx = np.argmax(thr);
    ttp = t[peakIdx];
    # Lagtime is found as the first timepoint such that the concentration exceeds 1% of the peak concentration
    for i in range(thr.size):
        if thr[i]>peak*0.01:
            break
    lagtime = t[i];
    # Numpy's gradient method is used for finding the maximum and minimum rates
    grad = np.gradient(thr, t);
    maxIncRate = np.max(grad);
    minDecRate = np.min(grad);
    stats = np.array([etp, peak, ttp, lagtime, maxIncRate, minDecRate]);
    return stats

#### Options ####
# Toggles for running separate parts of the analysis
plotThr = True;
runETPCorr = False;
runRateSA = False;
runICSA = False;


#### Thrombin Generation Plots ####
# Plot a thrombin generation curve for each model
if plotThr:
    maxt=1200; # Simulation time in seconds (1200s = 20 min)
    
    hockinK = hockin.getRates();
    danforthK = danforth.getRates();
    chatterjeeK = chatterjee.getRates();
    brummelK = brummel.getRates();
    bungayK = bungay.getRates();
    panteleevK = panteleev.getRates();
    tyurinK = tyurin.getRates();
    lakshmananK = lakshmanan.getRates();
    
    hockin.plotThr(hockinK,hockin.setIC(),maxt);
    danforth.plotThr(danforthK,danforth.setIC(),maxt);
    chatterjee.plotThr(chatterjeeK,chatterjee.setIC(),maxt);
    brummel.plotThr(brummelK,brummel.setIC(),maxt);
    bungay.plotThr(bungayK,bungay.setIC(),maxt);
    panteleev.plotThr(panteleevK,panteleev.setIC(),maxt);
    tyurin.plotThr(tyurinK,tyurin.setIC(),maxt);
    lakshmanan.plotThr(lakshmananK,lakshmanan.setIC(),maxt);


#### ETP Correlation ####
# Run the ETP correlation for 15 individuals and plot scatter plots
if runETPCorr:
    #### Load Rates ####
    # Preload the reaction rates for each model
    hockinK = hockin.getRates(includeExtras = True);
    danforthK = danforth.getRates(includeExtras = True);
    chatterjeeK = chatterjee.getRates(includeExtras = True);
    brummelK = brummel.getRates(includeExtras = True);
    bungayK = bungay.getRates(includeExtras = True);
    panteleevK = panteleev.getRates(includeExtras = True);
    tyurinK = tyurin.getRates(includeExtras = True);
    lakshmananK = lakshmanan.getRates(includeExtras = True);
    
    maxt=1200; # Simulation time in seconds (1200s = 20 min)
    baseIC = np.array([1.4e-6, 2e-8, 1e-8, 1e-8/100, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.4e-6, 2.5e-9]); #Baseline initial conditions/ factor levels. Order is II, V, VII, VIIa, VIII, IX, X, XI, AT, TFPI
    data = np.genfromtxt('paperData.csv', delimiter=',', skip_header=1); #Extract the patient data from the CSV file
    # Patient data is extracted as data[patient row index, factor/ETP column index]   row then column (person then data type)
    # Columns are ordered as ETP % (+TF), ETP % (-TF), TF (pM), II (%), V (%), VII (%), VIII (%), IX (%), X (%), XI (%), AT (%), TFPI (%)
    etp = np.zeros((15,8,2)); #15 patients, 8 models, 2 simulations per patient (with or without 5pM added TF)
    
    for i in range(15):
        # Find donor specific factor levels using baseIC. VIIa is used as 1% of VII
        donorIC = np.array([data[i,2]*1e-12, baseIC[0]*data[i,3]/100, baseIC[1]*data[i,4]/100, baseIC[2]*data[i,5]/100, baseIC[3]*data[i,5]/100, baseIC[4]*data[i,6]/100, baseIC[5]*data[i,7]/100, baseIC[6]*data[i,8]/100, baseIC[7]*data[i,9]/100, baseIC[8]*data[i,10]/100, baseIC[9]*data[i,11]/100]); #TF, II, V, VII, VIIa, VIII, IX, X, XI
        # Convert the donor specific factor levels into the form used for the individual models with 0 for all other concentrtions
        hockinIC = hockin.setIC(donorIC,includeExtras = True);
        danforthIC = danforth.setIC(donorIC,includeExtras = True);
        chatterjeeIC = chatterjee.setIC(donorIC,includeExtras = True);
        brummelIC = brummel.setIC(donorIC,includeExtras = True);
        bungayIC = bungay.setIC(donorIC,includeExtras = True);
        panteleevIC = panteleev.setIC(donorIC,includeExtras = True);
        tyurinIC = tyurin.setIC(donorIC,includeExtras = True);
        lakshmananIC = lakshmanan.setIC(donorIC,includeExtras = True);
        
        # Simulate each model using the getThr function. Extract ETP and store in the etp array
        (t,thr) = hockin.getThr(hockinK, hockinIC, maxt);
        sums = sumStat(t, thr);
        etp[i,0,0] = sums[0];
        
        (t,thr) = danforth.getThr(danforthK, danforthIC, maxt);
        sums = sumStat(t, thr);
        etp[i,1,0] = sums[0];
        
        (t,thr) = chatterjee.getThr(chatterjeeK, chatterjeeIC, maxt);
        sums = sumStat(t, thr);
        etp[i,2,0] = sums[0];
        
        (t,thr) = brummel.getThr(brummelK, brummelIC, maxt);
        sums = sumStat(t, thr);
        etp[i,3,0] = sums[0];
        
        (t,thr) = bungay.getThr(bungayK, bungayIC, maxt);
        sums = sumStat(t, thr);
        etp[i,4,0] = sums[0];
        
        (t,thr) = panteleev.getThr(panteleevK, panteleevIC, maxt);
        sums = sumStat(t, thr);
        etp[i,5,0] = sums[0];
        
        (t,thr) = tyurin.getThr(tyurinK, tyurinIC, maxt);
        sums = sumStat(t, thr);
        etp[i,6,0] = sums[0];
        
        (t,thr) = lakshmanan.getThr(lakshmananK, lakshmananIC, maxt);
        sums = sumStat(t, thr);
        etp[i,7,0] = sums[0];
        
        # Repeat the same process with 5pM of added TF
        donorIC[0] += 5e-12;
        hockinIC = hockin.setIC(donorIC,includeExtras = True);
        danforthIC = danforth.setIC(donorIC,includeExtras = True);
        chatterjeeIC = chatterjee.setIC(donorIC,includeExtras = True);
        brummelIC = brummel.setIC(donorIC,includeExtras = True);
        bungayIC = bungay.setIC(donorIC,includeExtras = True);
        panteleevIC = panteleev.setIC(donorIC,includeExtras = True);
        tyurinIC = tyurin.setIC(donorIC,includeExtras = True);
        lakshmananIC = lakshmanan.setIC(donorIC,includeExtras = True);
        
        (t,thr) = hockin.getThr(hockinK, hockinIC, maxt);
        sums = sumStat(t, thr);
        etp[i,0,1] = sums[0];
        
        (t,thr) = danforth.getThr(danforthK, danforthIC, maxt);
        sums = sumStat(t, thr);
        etp[i,1,1] = sums[0];
        
        (t,thr) = chatterjee.getThr(chatterjeeK, chatterjeeIC, maxt);
        sums = sumStat(t, thr);
        etp[i,2,1] = sums[0];
        
        (t,thr) = brummel.getThr(brummelK, brummelIC, maxt);
        sums = sumStat(t, thr);
        etp[i,3,1] = sums[0];
        
        (t,thr) = bungay.getThr(bungayK, bungayIC, maxt);
        sums = sumStat(t, thr);
        etp[i,4,1] = sums[0];
        
        (t,thr) = panteleev.getThr(panteleevK, panteleevIC, maxt);
        sums = sumStat(t, thr);
        etp[i,5,1] = sums[0];
        
        (t,thr) = tyurin.getThr(tyurinK, tyurinIC, maxt);
        sums = sumStat(t, thr);
        etp[i,6,1] = sums[0];
        
        (t,thr) = lakshmanan.getThr(lakshmananK, lakshmananIC, maxt);
        sums = sumStat(t, thr);
        etp[i,7,1] = sums[0];
        
        
        print("Completed donor "+str(i+1)+" of 15")
    
    # Plot ETP Correlation scatter plots
    modelNames = ["Hockin","Danforth","Chatterjee","Brummel","Bungay","Panteleev","Tyurin","Lakshmanan"];
    for i in range(len(modelNames)):
        plt.scatter(etp[:,i,:], data[:,[1,0]]);
        plt.title(modelNames[i]);
        plt.xlabel("Predicted ETP (M$\cdot$sec)");
        plt.ylabel("Experimental ETP (%)");
        plt.show();



#### Sensitivity Analysis ####
##  Reaction Rate Sensitivity  ##
if runRateSA:
    # Specify which model to perform the sensitivity analysis on
    import hockin as modelSA 
    
    maxt=1200; # Simulation time in seconds (1200s = 20 min)
    y = modelSA.setIC(); # Get model default initial condition vector
    k = modelSA.getRates(); # Get reaction rate vector
    parMult = np.linspace(0.5,1.5,11); # 50% to 150% scale factors
    (t,thr) = modelSA.getThr(k, y, maxt); # Simulate model with default initial concentrations to get the default summary statistics
    baseSums = sumStat(t,thr); 
    sums = np.zeros((baseSums.size,parMult.size)); # Array for storing summary statistics
    Sk = np.zeros(k.size); # Vector for storing sensitivities
    for i in range(k.size):
        for j in range(parMult.size):
            newK = modelSA.getRates(); # Get reaction rate vector and ...
            newK[i] = newK[i]*parMult[j]; # Scale ith rate
            (t,thr) = modelSA.getThr(newK, y, maxt); #Simulate model and ...
            sums[:,j] = sumStat(t,thr); # Store summary statistics
        Sk[i] = np.linalg.norm(np.std(sums,1)/baseSums); #Calculate sensitivity from summary statistics
        print("Completed parameter "+str(i+1)+" of "+str(k.size))
    Sk = Sk/np.sum(Sk); #Normalise sensitivities to sum to 1


##  Initial Concentration Sensitivity  ##
if runICSA:
    import hockin as modelSA;
    
    maxt=1200; # Simulation time in seconds (1200s = 20 min)
    factorNames = ["TF","FII","FV","FVII","FVIIa","FVIII","FIX","FX","FXI","AT","TFPI"];
    y = modelSA.setIC(); # Get model default initial condition vector
    k = modelSA.getRates(); # Get reaction rate vector
    parMult = np.linspace(0.5,1.5,11); # 50% to 150% scale factors
    (t,thr) = modelSA.getThr(k, y, maxt); # Simulate model with default initial concentrations to get the default summary statistics
    baseSums = sumStat(t,thr);
    sums = np.zeros((baseSums.size,parMult.size)); # Array for storing summary statistics
    Sy = np.zeros(len(factorNames)); # Vector for storing sensitivities
    for i in range(len(factorNames)):
        for j in range(parMult.size):
            newY = np.array([10e-12, 1.4e-6, 2e-8, 1e-8, 1e-8/100, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.4e-6, 2.5e-9]); # Set baseline concentrations
            newY[i] = newY[i]*parMult[j]; # Scale ith initial concentration
            newY = modelSA.setIC(newY); # Convert scaled initial concentration vector to correct form for model simulation
            (t,thr) = modelSA.getThr(k, newY, maxt); #Simulate model and ...
            sums[:,j] = sumStat(t,thr); # Store summary statistics
        Sy[i] = np.linalg.norm(np.std(sums,1)/baseSums); #Calculate sensitivity from summary statistics
        print("Completed "+str(factorNames[i]));
    Sy = Sy/np.sum(Sy); #Normalise sensitivities to sum to 1


