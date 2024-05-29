import math
import numpy as np
import pandas as pd
import scipy as sp


class Analyse():
    """analyse dataset from chem 461"""

    def __init__(self):
        """Initializing"""

    def getbase(self, ds, hibar, nfitpts=30, adjust=25):
        """
        determine a baseline from experiment 42 data
        Input:
        ds - n by 2 array, data of absorbance vs wavenumber
        hibar - scalar, height-value below which a baseline is interpolated
        nfitpts - scalar, number of points used to interpolate a baseline
        Output:
        baseline = n by 1 array
        """
        ## parse dataset
        dsarr = np.array(ds)
        x = dsarr[:, 0]
        y = dsarr[:, 1]
        npts = len(x)
        ## define a uniform increment, dN, across the length of the dataset
        dN = npts/nfitpts
        K = np.linspace(0.5, (nfitpts - 0.5), nfitpts)*dN
        K = K.round().astype(int)
        ## determine the sign of the derivative across the dataset
        KI = K[1:-1]
        dK = np.full(len(KI), dN/adjust)
        dK = dK.round().astype(int)
        MD = np.sign(y[KI + dK] - y[KI])
        MD = np.append([-1], MD)
        MD = np.append(MD, [1])
        MD = MD.astype(int)
        ## nudge the y-values downards until they are below at hibar
        nudgefac = np.rint(0.0004*npts)
        if (nudgefac < 1): nudgefac = 1
        idx = np.where(y[K] > hibar)[0]
        while(idx.size):
            ## if y is increasing then move in the reverse direction
            ## so use MD to determine the slope of y
            K[idx] = K[idx] - MD[idx]*nudgefac
            idx = np.where(y[K] > hibar)[0]

        ## now keep nudging the y-values until a little upwards
        idx = np.linspace(0, nfitpts)[0]
        while(idx.size):
            KN = K - MD*nudgefac
            KN = KN.round().astype(int)
            idx = np.where((y[KN] - y[K]) < 0)[0];
            K[idx] = KN[idx]

        ## now use the good indexes, which represents the fit points
        ## that are spread out uniformly across the dataset and below hibar, 
        ## to interpolate a baseline
        xk = np.flip(x[K])  # arrange values in increasing sequence
        yk = np.flip(y[K])
        yb = sp.interpolate.CubicSpline(xk, yk)
        baseline = yb(x)

        return baseline

    @staticmethod
    def caloexp(time, rate_heatgain, rate_heatloss, rate_reaction, DeltaT,
        tstart=145, Tstart=291, Troom=293):
        """
        The function is defined for scipy's fitting procedure in order to
        optimized the temperature rise T as a function of time, t
        This function models the temperature rise of the Paar Calorimeter in
        experiment 3:
        dT/dt = rate_heatgain - rate_heatloss*(T - Troom) + rate_reaction*T_rxn
            = rate_heatgain - rate_heatloss*(T - Troom) 
            + rate_reaction*DeltaT*exp(-rate_reaction*(t - tstart))
        Input:
        time - n by 1 array [sec]
        parameters to optimize:
        rate_heatgain - rate of heat gain from the stirrer [K/sec]
        rate_heatloss - the rate of heat loss [1/sec]
        rate_reaction - rate of heat gain due to reaction [1/sec]
            heat is on from t=0  to  t= 1/(rate_reaction/2)
        DeltaT - net rise of temperature due to chemical reaction [K]
        e.g. of guess parameters to optimize: pars = [0.003 0.001 0.02 5]
        default paremeters as karg
        tstart=30 - begining of reaction (heat of combution) [sec]
        Tstart=292 - the first temperature in the data [K]
        Troom=290 - temperature of the reservoir (or room) [K]
        Output:
        Tcurve - n by 1 array, the temperature range [K]
        """
        # Set up the heating region that goes from tstart 
        idheat = [i for i, j in enumerate(time) if j >= tstart]
        Zt = [time[i] - tstart for i in idheat]
        # The temperature in the Zt range from start to end
        TZ = [DeltaT*rate_reaction*(math.exp(-i*rate_heatloss) -
           math.exp(-i*rate_reaction))/(rate_reaction - rate_heatloss) 
           for i in Zt]
        Tmod = [0] * len(time)
        for i, j in zip(idheat, TZ):
            Tmod[i] = j
    
        # the heat in and heat loss term is added to the Tmod 
        # (due to reaction heating)
        Eloss = [math.exp(-i*rate_heatloss) for i in time]
        Zloss = []
        if(rate_heatloss < 1e-7):  # take the limit if user sets rate_loss=0
            Zloss = time
        else:
            Zloss = [(1 - i)/rate_heatloss for i in Eloss]
    
        TI = [i*rate_heatgain + j*Tstart + (1 - j)*Troom
            for i, j in zip(Zloss, Eloss)]
        Tcurve = [i + j for i, j in zip(TI, Tmod)]
    
        return Tcurve

    def rst_exp9(self, var_pars, FHP):
        """
        A function to simulate the liquid-vapor phase-diagram for two liquids
        Using Regular Solution Theory (RST) for the data from Exp 9, Chem 461
        Needed: the boiling temperature (K) and the Heat of Vaporization
        (KJ/mol) of two liquids A and B
        Input:
        FHP - The Flory Huggins (attraction) parameter must be set
          FHP = 0 give ideal solution (Raoult's Law) phase diagram
          FHP ranges from about -1 to 1.
        vap_pars - the parameters of vaporization, 
            pars = [TboilA, enthalpy_vapA, TboilB, enthalpy_vapB]
            TboilA - ref. boiling point of pure A (preferably at 1 atm) [Kelvin]
            TboilB - ref. boiling point of B
            enthalpy_vapA - heat of vaporization of pure A [KJ/mol]
            enthalpy_vapB - heat of vaporization of pure B
            eg:  var_pars  =  [383.8, 33.18, 390.6, 43.29];
        Output:
          xB - the mole fraction of B in liquid's phase
          yB - the mole fraction of B in vapor's phase and 
          Tbest - the equilibrium boiling temperature at which 
            xB and yB are determined
        NB: The total pressure is set to 1 atm,  but can be changed.
        """
        # constants
        npts = 200;
        xB = np.linspace(0, 1, npts)
        xA = 1 - xB
        R = 8.314  # gas constant in J/(mol K)
        TboilA = var_pars[0]; enthalpy_vapA = var_pars[1]*1e3/R/TboilA;
        TboilB = var_pars[2]; enthalpy_vapB = var_pars[3]*1e3/R/TboilB;
        PboilA = 1; PboilB = 1;  # [atm]
        P0 = 1; T0 = 360; Ptot = 1;

        # The Clausius Clapeyron Relation
        PA_pure = lambda Tx : np.exp(enthalpy_vapA*(Tx - TboilA)/Tx)*P0
        PB_pure = lambda Tx : np.exp(enthalpy_vapB*(Tx - TboilB)/Tx)*P0
        # xA = (Ptot - PB_pure(T0))/(PA_pure(T0) - PB_pure(T0))
        # yA = (xA*PA_pure(T0))/Ptot

        # Using regulary solution theory, we calculate the partial pressures
        # method 1: calculate at fixed T0
        # lnPPB = np.log(xB) + FHP*(xA**2);
        # PB = np.exp(lnPPB)*PB_pure(T0)
        # yB = xB*np.exp(FHP*(xB**2))*PB_pure(T0)/Ptot
        # method 2:
        # Find the equilibrium temperature for the liquid-vapor mole fractions
        # First, extend the temperature range to find the answer
        # the extended range should be ten times the magnitude of FHP
        pad = 10*np.rint(FHP);
        range_T = np.array(
            np.minimum(TboilA, TboilB) - pad, np.maximum(TboilA, TboilB) + pad
            )
        Tbest = np.zeros(npts)
        for k in range(0, npts):
            PA = lambda Tx : xA[k]*np.exp(FHP*(xB[k]**2))*PA_pure(Tx)
            PB = lambda Tx : xB[k]*np.exp(FHP*(xA[k]**2))*PB_pure(Tx)
            find_T = lambda Tx : PA(Tx) + PB(Tx) - Ptot
            Tbest[k] = sp.optimize.newton(find_T, range_T)

        yB = xB*np.exp(FHP*(xA**2))*PB_pure(Tbest)/Ptot

        return xB, yB, Tbest


class Opener():
    """read dataset from chem461"""

    def __init__(self):
       """Initializing""" 

    def getdata(self, dnfn):
        """
        method to read text/csv data using pandas' package
        Input:
        dnfn - string, path to file
        Output:
        data - array, data array
        """
        # the dataframe of pandas
        df = pd.read_csv(dnfn)
        # convert dataframe into numpy array
        ds = df.to_numpy()

        return ds

    def get_txt_data(self, dnfn, headerlines=2):
        """
        method to read the cvs data for exp 42
        """
        ds = []
 
        with open(dnfn, 'rt') as f:
 
            # skip header lines
            for i in range(headerlines):
                next(f)
 
            # readline and convert string to number array
            for line in f:
                linetrunc_x = line[0:6].strip()
                linetrunc_y = line[8:].strip()
                # convert string array to number array
                data = [float(linetrunc_x), float(linetrunc_y)]
                data_abs = [data[0], -np.log10(data[1]/100)] 
                ds.append(data_abs)
 
        return ds, dnfn
