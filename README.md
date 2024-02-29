# Chem 461 - Physical Chemistry Laboratory
This class gives practical, hands-on experience in applying physical chemistry
concepts in the laboratory. It provides practice in professional scientific
writing and editing. It lays out more in-depth (and less "cookbook") approaches
to analyzing data. This repository focuses more on the data analysis as it
provides an introduction to statistically analyzing data through scientific
computing.

## Reading the data

The first step towards data analysis is to read the data into a variable name 
(or data structure). Use the method called `uwpchem.Opener.getdata()`.

## Experimental Analysis

### Exp 42: Infrared Spectra of HCl and DCl
This experiment determines the vibrational and rotational transition energies
of HCl, using a harmonic oscillator model (with anharmonic corrections) for the
IR spectra. Following, the thermodynamic and statistical mechanical quantities 
of the HCl system are estimated. This estimated quantities of HCl are compared
to that of DCl to determine the isotope effect between the systems. After
reading the dataset, perform a baseline correction to clean and smoothen the
baseline of the dataset. For such baseline correstion use the method called 
`uwpchem.Analyse.getbase()`. This procedure is demonstrated in the jupter
notebook called "Demo.ipynb".

### Exp 1: Spectrometry and Chemical Kinetics
This experiment estimates the rate constant for the product and intermediate
species of penicillin hydrolysis. The reaction rate is first-order, thus
fitting an exponential curve to the absorbance versus time data estimates
the parameters of the initial concentration, final concentration and the rate
constant. A fitting procedure is demonstrated in the jupyter notebook called 
"tutorial.ipynb".

### Exp 3: Heats of Combustion
This experiment determines the thermodynamic quantities of the combustion of
cyclohexane, cyclohexene, benzene, hexane, to study the effects of resonance
stability. The estimated heat capacity of a calorimetry and the temperature
difference is used. The change in temperature is estimated based on the model
that is used to fit the temperature vs time dataset. This model as a fit
function given by the method called `uwpchem.Analyse.caloexp()`. An example of 
this fitting procedure is demonstrated in the jupyter notebook called
"Demo.ipynb".

### Exp 9: Liquid-Vapor Equilibrium in Binary Systems
This experiment estimates the Flory-Huggins parameter, for the Regular
Solution Theory, that is used as a correction to Raoult's law, which fails to
describe the phase-diagram for azeotropes. The model of the Regular Solution
theory is used to fit the phase-diagram. This model as a fit function is given
by the method called `uwpchem.Analyse.rst()`. An example of this fitting
procedure is demonstrated in the jupyter notebook called "Demo.ipynb".

## Tutorials
The jupyter notebook called "tutorial.ipynb" provides instructions on importing
a package, manipulating a data, plotting a data, fitting a data, calculating
and plotting associated error bars based on the estimated errors of a fit,
and importing the `uwpchem` package and its modules `Opener` and `Analyse`.

Further instructions on importing the `uwpchem` package and its modules
`Opener` and `Analyse` are given in the jupyter notebook called "Demo.ipynb".
This notebook also demonstrate how to use the associated methods in the
modules: `getdata()` for reading the data; `getbase()` for interpolating a
baseline for baseline correction on the IR spectrum; `caloexp()` for fitting 
the data from the calorimetry; and `rst()` for fitting the data from the
phase-diagram of azeotropes.

<!--- 
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
-->
