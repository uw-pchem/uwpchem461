# Chem 461 - Physical Chemistry Laboratory
This class gives practical, hands-on experience in applying physical chemistry concepts in
the laboratory. It provides practice in professional scientific writing and editing. It 
lays out more in-depth (and less "cookbook") approaches to analyzing data. This repository
focuses more on the data analysis as it provides an introduction to statistically analyzing
data through scientific computing.

## Reading the data

The first step towards data analysis is to read the data into a variable name 
(or data structure). The following modules can be used:
```
get_OO_Data.py
```

## Experiments

### Exp 42: Infrared Spectra of HCl and DCl
This experiment estimates the vibrational and rotational transition energy of
HCl, using a harmonic oscillator model (with anharmonic corrections) for the
energy spectrum. Following, the thermodynamic and statistical mechanical quantities 
of the HCl system are estimated. This estimated quantities of HCl are compared
to that of DCl to determine the isotope effect between the systems. The
analysis is demonstrated in:
```
exp42_analyze_data.ipynb
```
The following code is also used for baseline corrections of the IR spectra:
```
getbase.py
```
This code is a module within the class object for this experiment.

### Exp 1: Spectrometry and Chemical Kinetics
This experiment estimates the rate constant for the product and intermediate
species of penicillin hydrolysis. The reaction rate is first-order, thus
fitting an exponential curve to the absorbance versus time data estimates
certain parameters of the initial concentration, final concentration and the
rate constant. Analysis for the fitting procedure is demonstrated in tutorials
1 and 2. The following code is used for block averaging the large dataset:
```
blockavgM.py
```
This code is a module within the class object for this experiment.

### Exp 3: Heats of Combustion
This experiment determines the thermodynamic quantities of the combustion of
cyclohexane, cyclohexene, benzene, hexane, to study the effects of resonance
stability. The estimated heat capacity of a calorimetry and the temperature
difference is used. The change in temperature is estimated based on the fit
model described in the following code:
```
caloexp.py
```
This code is a module within the class object for this experiment.

### Exp 9: Liquid-Vapor Equilibrium in Binary Systems
This experiment estimates the Flory-Huggins parameter for the Regular
Solution Theory that is used as a correction to Raoult's law, which fails to
describe the phase-diagram for azeotropes. The fit model procedure for the 
Flory-Huggins solution theory is described in the follwing code:
```
rst.py
```
This code is a module within the class object for this experiment.

## Tutorials
The tutorial (jupyter) notebooks demonstrate the techniques of the modules
used for analysing the dataset of the experiments. For demonstrations of how to
read a file see:
```
tutorial1.ipynb
```
This tutorial also introduces the least square fitting procedure, based on the
module:
```
zlstsq.ipynb
```
For demonstrations on how to write functions (e.g the exponential functions)
used for fitting data (e.g. absorbance vs time data) see
```
tutorial2.ipynb
```
For all the estimated quanties derived from the fit parameters, there must be
an error analysis. In this class, guassian error propagation is used, which is
introduced in
```
tutorial3.ipynb
```

## Outline for Chem461 Winter Quarter
### Tasks
    1. Set up github page
    2. Set up and test-install environment: miniconda, jupyter notebook
    3. Translate the following matlab codes:
        a. codes for readind data file
            i. get_OO_Data.m
            ii. get_Putty_Data.m
            iii. get_Text_Data.m
        b. code for Exp. 3
            i. CaloExp.m
        c. code for Exp. 9
            i. RST_Exp9.m
            ii. Find_T.m
        d. code for Exp. 42
            i. GetBase.m
        e. code for Exp. 1
            i. BlockAvgM.m
        f. code for Tutorials
            i. Zlstsq.m
    4. Translate the following documents for python:
        a. Experiment 42 Analyze Data.pdf
        b. Week 1 Tutorial.pdf
        c. Weak 2 Tutorial.pdf
        d. Week 3 Tutorial.pdf
### Time Allocated for Tasks
    Week 1: Orientation, set hours
    Week 2: Make outline, do Task 1
    Week 3: do Task 2 // Task 4a - with Sarah's help
    Week 4: do Task 3a // Task 4a (with Sarah's help)
    Week 5: do Task 3b // Task 4b (with Sarah's help)
    Week 6: do Task 3c // Task 4b (with Sarah's help) 
    Week 7: do Task 3d // Task 4c (with Sarah's help)
    Week 8: do Task 3e // Task 4c (with Sarah's help)
    Week 9: do Task 3f // Task 4d (with Sarah's help)
