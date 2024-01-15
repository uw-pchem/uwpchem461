import numpy as np

class Opener():
    """read dataset from chem461"""
    def __init__(self):

    return


    def get_oo_data(self, dnfn, start=13, headerlines=0):
        """
        method to read ocean optics data
        Input:
        dnfn - path to file
        start - skips first thirteen chars of each line in ocean optics data
        Output:
        ds - dataset
        dnfn - path to file returned
        """
        ds = np.arr([])

        with open(dnfn, 'rt') as f:

            # skip header lines
            for i in range(headerlines)
                next(f)

            # readline and convert string to number array
            for line in f:
                linetrunc = line[start:]
                # convert string array to number array
                linedata = np.asarray(linetrunc, dtype=float)
                ds.append(linedata)

        return ds, dnfn


    def get_txt_data(self, dnfn):
        """
        method to read the cvs data for exp 42
        """
        ds, _ = self.get_oo_data(dnfn, start=0, headerlines=2)

        return ds, dnfn
