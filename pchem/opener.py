import numpy as np

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
                ds.append(data)

        return ds, dnfn
