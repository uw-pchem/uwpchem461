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
        determine a baseline from experiment 42 dataset
        Input:
        ds - n by 2 array, dataset
        hibar - scalar, height value below which baseline is established
        nfitpts - scalar, number of fit points used to interpolate a baseline
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

    def caloexp(self, time, pars):
        """
        This models the temperature rise of the Paar Calorimeter
        dT/dt = RIn - RLoss*(Yloop - Tres) + beta*DeltaT*exp(-beta*(time -
            tstart))
        It is defined for scipy's fitting procedure in order to optimized
        that temperature rise as a function of time
        Input:    
        time - n by 1 array [sec]
        pars - parameters to optimize in pars = [RIn RLoss beta DeltaT];
            rate_in - the stirrer temperature rate [K/sec]
            rate_loss - the rate of loss [1/sec]
            beta - rate of heat injection due to reaction [1/sec]
                heat is on from t=0  to  t= 1/beta/2
            DeltaT - net rise of temperature due to chemical reaction [K]
            e.g. pars = [0.003 0.001 0.02 5]
        default paremeters as karg
        tstart=30 - begining of reaction (heat of combution) [sec]
        Tstart=292 - the first temperature in the data [K]
        Tres=290 - temperature of the reservoir (or room) [K]
        Output:
        temperature_curve - n by 1 array, the temperature range [K]
        """
        # Parse the input
        rate_heatgain = pars[0]
        rate_heatloss = pars[1]
        beta = pars[2]
        DeltaT = pars[3]
        var1 = input("What is the starting time (seconds) of the combustion? ")
        var2 = input("What is the temperature (Kelvin) at the starting time? ")
        var3 = input("What is the room temperature (Kelvin)? ")
        tstart = float(var1)
        Tstart = float(var2)
        Tres = float(var3)
        # Set up the heating part that goes from tstart 
        # this is the heating region.
        idheat = [i for i, j in enumerate(time) if j >= tstart]
        # the set of indicies for the time that cover the heating range.
        # Nheat = length(idheat);  % number of points in that range.
        Zt = [time[i] - tstart for i in idheat]
        # The temperature in the Zt range from start end
        TZ = [beta*DeltaT*(math.exp(-i*rate_heatloss) -
           math.exp(-i*beta))/(beta - rate_heatloss) for i in Zt]
        Tmod = [0] * len(time)
        for i, j in zip(idheat, TZ):
            Tmod[i] = j
    
        # the heat in and heat loss term is added to the Tmod 
        # (due to reaction heating)
        Eloss = [math.exp(-i*rate_heatloss) for i in time]
        Zloss = []
        if( rate_heatloss < 1e-7):  # take the limit if user sets rate_loss=0
            Zloss = time
        else:
            Zloss = [(1 - i)/rate_heatloss for i in Eloss]
    
        TI = [i*rate_heatgain + j*Tres + (1 - j)*Tres
            for i, j in zip(Zloss, Eloss)]
        temperature_curve = [i + j for i, j in zip(Tmod, TI)]
    
        return temperature_curve


class Opener():
    """read dataset from chem461"""

    def __init__(self):
       """Initializing""" 

    def get_data(self, dnfn):
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


#     def get_oo_data(self, dnfn, start=13, headerlines=0, d=","):
#         """
#         method to read ocean optics data
#         Input:
#         dnfn - path to file
#         start - skips first thirteen chars of each line in ocean optics data
#         Output:
#         ds - dataset
#         dnfn - path to file returned
#         """
#         ds = np.array([[0, 0]])
# 
#         with open(dnfn, 'rt') as f:
# 
#             # skip header lines
#             for i in range(headerlines):
#                 next(f)
# 
#             # readline and convert string to number array
#             for line in f:
#                 linetrunc = line[start:].strip()
#                 # convert string array to number array
#                 linedata = np.fromstring(linetrunc, dtype=float, sep=d)
#                 np.append(ds, [linedata])
# 
#         return ds[1:,1:], dnfn


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


# class Fitpolynomial():
#     """
#     Use the method of least-squares to fit a polynomial to a dataset
#     """
# 
#     def __init__(self):
#         """
#         par -  k by 1 array, the elements are the fit parameters,
#             where n is the degree of the polynomial
#         yf - n by 1 array, the estimated y-values based on the fit 
#             parameters:
#             yf = M*par
#             yf = par[0]*x^k + par[1]*x^(k - 1) ... + par[k]*x^0
#         M - n by k array, model matrix:
#             M = [x^k + x^(k - 1) ... + x^0]
#         V - k by k array, the reduced "variance-covariance" matrix:
#             V = M^(-1)*M
#         stefit - standard errors of the fit's estimated y-values
#         stepar - standard errors of the parameters
#         steyf - standard deviation of the fit's estimated y-values
#         """
#         self.par = []
#         self.stepar = []
#         self.stefit = []
#         self.steyf = []
# 
#     def zlstsq(self, xdata, yata, degree):
#         """
#         the least-squares fitting procedure to the dataset:
#         Input:
#         xdata - n by 1 array, x-values of dataset
#         ydata - n by 1 array, y-values of dataset
#         degree - scalar, the degree of the polynomial
#         """
#         x = np.array(xdata); y = np.array(ydata)
#         # the model matrix: M = [x^k + x^(k - 1) ... + x^0]
#         M = np.array([x**k for k in reversed(range((degree + 1)))])
#         M = M.transpose()
#         # this is the heart of this method, which depends on the fact that
#         # there exists an inverse matrix of the model matrix
#         inverseM = linalg.inv(M)


