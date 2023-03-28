# Mathematical Models Of Coagulation - Are We There Yet?
 Supplemental Python code for "Mathematical models of coagulation - are we there yet?"


<!-- ABOUT THE PROJECT -->
## About The Project
This project contains the supplemental Python code for "Mathematical models of coagulation - are we there yet?". This Python code (written in Python ver. 3.7.9) is a translation from the original Matlab code used for the paper. It has the ability to both run the 8 models used for the paper and to run the correlational and both sensitivity analysis methods introduced and used in the paper.



<!-- GETTING STARTED -->
## Getting Started
To use this code follow the installation steps below.

### Prerequisites and Version Numbers
This code is written and tested in Python ver. 3.7.9 however any Python version (not including Python 2 or earlier) should be suitable. In addition to this Numpy (ver. 1.19.3), Scipy (ver. 1.7.0) and Matplotlib (ver 3.4.2) are also used throughout this project.

### Installation
There are 2 methods for downloading the project.
1. Clone the repo
   ```sh
   git clone https://github.com/mjowen/MathematicalModelsOfCoagulation-AreWeThereYet.git
   ```
2. Download as ZIP
    - Visit [https://github.com/mjowen/MathematicalModelsOfCoagulation-AreWeThereYet](https://github.com/mjowen/MathematicalModelsOfCoagulation-AreWeThereYet)
    - Click Code -> Download ZIP
    - Unzip and set the folder "MathematicalModelsOfCoagulation-AreWeThereYet-main" as the working directory in Python


<!-- USAGE EXAMPLES -->
## Usage
### Using *main.py*
*main.py* runs the analysis presented in the paper. By default it will plot thrombin generation curves for each model but this can be changed in the *Options* section. By toggling the variables *plotThr*, *runETPCorr*, *runRateSA* and *runICSA* to *True* or *False* you can control whether running *main.py* plots the thrombin generation curves, runs the ETP correlation analysis, runs the reaction rate sensitivity analysis or the initial condition (coagulation factor) sensitivity analysis. 

The ETP correlation analysis will read in the data of 15 donors, given in *paperData.csv*, and will plot the scatterplots for the simulations.

Both sensitivity analyses methods return a list of sensitivities (*Sk* for reaction rate sensitivities and *Sy* for coagulation factor sensitivities). The coagulation factor sensitivities fall in the standard order of: TF, II, V, VII, VIIa, VIII, IX, X, XI, AT and TFPI. The reaction rate sensitivities return the sensitivity of each reaction rate (in the order used in the code) and will need to be cross-referenced against the model's ODE function to find which reaction is corresponds to.

### How to use the individual models
Each model can be also used individually. The models are stored in there own files and can be freely imported in Python, for example:
```python
import hockin
```
#### setIC
To use a model an initial condition vector has to be made which can be done using the setIC fucntion. This will translate a vector of coagulation factors (in the standard order of: TF, II, V, VII, VIIa, VIII, IX, X, XI, AT and TFPI) into an initial condtion vector to be used by the model. If no coagulation factors are specified, the default coagulation factor concentration will be used. 

#### getRates
Calling the models getRates function will return the reaction rates. Should you want to make any changes to the reaction rates, call the getRates function, make any changes to the returned variable and pass this into the later functions. All reaction rates use the default units of moles and seconds (regardless of the units used in the original description of the model).

#### getThr
The getThr function will solve the model and return a tuple of a time vector and a thrombin concentration vector. It requires inputs of an initial condition vector, a reaction rate vector and a maximum time for the simulation (simulated time rather than real time) specified in seconds (all analysis from the paper uses 20 min = 1200 sec which it the length of time the assay was run for in our data).

#### plotThr
The plotThr function uses all the same inputs as the getThr function but instead of returning the thrombin concentration a thrombin generation curve is plotted.

#### test
The test function requires no inputs and will go through the steps of forming an initial condition vector (of default factor levels), forming the reaction rate vector and calling the plotThr function to plot a default thrombin generation curve.

#### Example
Below is an example for solving the Hockin model for 15pM of TF (rather than the default of 10pM) as well as a change in the reaction rate for **TF + VII -> TF==VII** (reaction rate index 1) and plotting the resulting thrombin generation curve.
```python
import numpy as np
import hockin
y = hockin.setIC(np.array([15e-12, 1.45e-6, 2e-8, 1e-8, 1e-8/100, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.45e-6, 2.5e-9]));
k = hockin.getRates();
k[1] = k[1]*2;
hockin.plotThr(k,y,1200)
```

#### Solving for factors other than thrombin
Should you wish to find factor concentration curves for factors other than thrombin you can change the index in the getThr function to the index of the factor you want to measure (indicies for the factors are given at the start of each model file). Note: The factor concentrations are not inclusive, ie if you wish to measure levels of Xa in the Hockin model it may not be sufficient to only use the index 5 (the index of Xa), if you wish to include levels of TF:VIIa:Xa or levels of Xa:Va then this needs to be specified explicitly. An example line for calculating total Xa (not included TFPI or AT inhibited) is given below.
```python
return (sol.t,sol.y[5]+sol.y[9]+sol.y[22]+sol.y[23])
```
should replace the line
```python
return (sol.t,sol.y[6])
```

<!-- LICENSE -->
## License
Licensed under the [GNU GENERAL PUBLIC LICENSE](LICENSE)
    
    Copyright (C) 2021  Matthew John Owen
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    


<!-- CONTACT -->
## Contact

Matt J. Owen - matthew.owen@nottingham.ac.uk

Project Link: [https://github.com/mjowen/MathematicalModelsOfCoagulation-AreWeThereYet](https://github.com/mjowen/MathematicalModelsOfCoagulation-AreWeThereYet)

