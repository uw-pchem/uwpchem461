# unit test for the code on regular solution theory
import numpy as np
## import os.path
## from pathlib import Path

from chem461.pchem import Opener, Analyse


def test_rst():
    """
    check the datatype of the liquid and vapor mole fractions of a molecules B
    as well as the equilibrium boiling temperature at which the mole fractions
    are determined
    """
    ## abspath = os.path.abspath("data/exp3_dataset.csv")
    ## dnfn = Path(abspath)
    analyse = Analyse()
    var_pars = [383.8, 33.18, 390.6, 43.29]; FHP = 0;
    xB, yB, Tvap = analyse.rst_exp9(var_pars, FHP)
    # check the data-type
    try:
        truthvalue = np.issubdtype(Tvap.dtype, np.floating)
        assert(truthvalue)
    except AssertionError as a:
        raise a("The data-type of the temperature's array should be floating")
    else:
        print("The boiling temperature's array is floating.")
    ## check shape of file
    try:    
        assert(yB.shape[0] == yB.size)
    except AssertionError as a:
        raise a("The vapor's mole fraction should be an array with shape n by 1")
    else:
        print("The vapor's mole fraction is an array with shape ", yB.shape)

    return
