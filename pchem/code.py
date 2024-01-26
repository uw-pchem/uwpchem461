import numpy as np
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
        # print(idx.size)
        # print(K.size)
        # print(MD.size)
        # print(idx); print(K); print(MD)
        # print(y[K])
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

    def caloexp(pars, time, tstart=30, Tstart=292, Tres=290):
        """
        This models the temperature rise of the Paar Calorimeter
        dT/dt = RIn - RLoss*(Yloop - Tres) + beta*DeltaT*exp(-beta*(time -
            tstart))
        It is defined for scipy's fitting procedure in order to optimized
        that temperature rise as a function of time
        Input:    
        time - n by 1 array [sec]
        pars - parameters to optimize in pars = [RIn RLoss beta DeltaT];
            RIn - the stirrer temperature rate [K/sec]
            RLoss - the rate of loss [1/sec]
            beta - rate of heat injection due to reaction [1/sec]
                heat is on from t=0  to  t= 1/beta/2
            DeltaT - net rise of temperature due to chemical reaction [K]
            e.g. pars = [0.003 0.001 0.02 5]
        default paremeters as karg
        tstart=30 - begining of reaction (heat of combution) [sec]
        Tstart=292 - the first temperature in the data [K]
        Tres=290 - temperature of the reservoir (or room) [K]
        """
    
        # Set up the heating part that goes from tstart  this is the heating region.
        idheat = [i for i, j in enumerate(time) if j > tstart]
        # the set of indicies for the time that cover the heating range.
        # Nheat = length(idheat);  % number of points in that range.
        Zt = [time[i] - tstart for i in idheat]
        # The temperature in the Zt range from start end
        TZ = [beta*DeltaT*(math.exp(-RLoss*i) - math.exp(-beta*i))/(beta-RLoss) for i in Zt]
        Tmod = [0] * len(time)
        for i, j in zip(idheat, TZ):
            Tmod[i] = j
    
        # the heat in and heat loss term is added to the Tmod (due to reaction heating)
        Eloss = [math.exp(-RLoss*i) for i in time]
        Zloss = []
        if( RLoss < 1e-7):  # take the limit if user sets RLoss=0
            Zloss = time
        else:
            Zloss = [(1-i)/RLoss for i in Eloss]
    
        TI = [RIn*i + Tres*j + Tres*(1-j) for i, j in zip(Zloss, Eloss)]
        return [i + j for i, j in zip(Tmod, TI)]


class Opener():
    """read dataset from chem461"""

    def __init__(self):
       """Initializing""" 

    def get_oo_data(self, dnfn, start=13, headerlines=0, d=","):
        """
        method to read ocean optics data
        Input:
        dnfn - path to file
        start - skips first thirteen chars of each line in ocean optics data
        Output:
        ds - dataset
        dnfn - path to file returned
        """
        ds = np.array([[0, 0]])

        with open(dnfn, 'rt') as f:

            # skip header lines
            for i in range(headerlines):
                next(f)

            # readline and convert string to number array
            for line in f:
                linetrunc = line[start:].strip()
                # convert string array to number array
                linedata = np.fromstring(linetrunc, dtype=float, sep=d)
                np.append(ds, [linedata])

        return ds[1:,1:], dnfn


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
