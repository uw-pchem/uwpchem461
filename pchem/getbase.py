import numpy as np
import scipy as sp

Class Analyse():
    """analyse dataset from chem 461"""

    def __init__(self):
        """Initializing"""

    def getbase(self, ds, hibar, nfitpts=17):
        """
        determine a baseline from experiment 42 dataset
        Input:
        ds - n by 2 array, dataset
        hibar - scalar, height value below which baseline is established
        nfitpts - scalar, number of points to be used for fitting a baseline
        Output:
        baseline = n by 1 array
        """
        ## parse dataset
        x = ds[:,0]
        y = ds[:,1]
        npts = len(x)
        ## define a uniform increment, dN, across the length of the dataset
        dN = npts/nfitpts
        K = np.linspace(0.5, (nfitpts - 0.5), nfitpts)*dN
        K = K.round().astype(int)
        ## determine the sign of the derivative across the dataset
        dK = np.full(nfitpts, 2, dtype=int)
        KI = K[1:-1]
        MD = np.sign(y[KI + dk] - y[KI])
        MD = np.append([-1], MD)
        MD = np.append(MD, [1])
        ## nudge the y-values downards until they are below at hibar
        nudgefac = np.rint(0.0004*npts)
        if (nudgefac < 1): nudgefac = 1
        idx = np.where(y[K] > hibar)
        while(idx.size):
            ## if y is increasing then move in the reverse direction
            ## so use MD to determine the slope of y
            K[idx] = K[idx] - MD[idx]*nudgefac
            idx = np.where(y[K] > hibar)

        ## now keep nudging the y-values until a little upwards
        idx = np.linspace(0, nfitpts)
        while(idx.size):
            KN = K - MD*nudgefac
            idx = np.where((y[KN] - y[K]) < 0);
            K[idx] = KN[idx]

        ## now use the good index points, which represents the fit points
        ## that is spread out uniformly across the dataset and below hibar, 
        ## to interpolate a baseline
        yb = sp.interpolate.CubicSpline(xk, yk)
        baseline = yb(x)

        return baseline
