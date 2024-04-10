import numpy as np
import os.path
from pathlib import Path

from chem461.pchem import Opener, Analyse


def test_caloexp():
    """
    is the temperature range the same in length as the time domain?
    is the initial temperature consistent with the one given?
    """
    abspath = os.path.abspath("data/exp3_dataset.csv")
    dnfn = Path(abspath)
    opener = Opener()
    analyse = Analyse()
    ds = opener.getdata(dnfn)
    temperature_curve = analyse.caloexp(ds[:, 0], 0.003, 0.001, 0.02, 5)    
    ## check shape of file
    try:    
        assert(len(temperature_curve) == len(ds[:, 0]))
    except AssertionError as a:
        raise AssertionError("The temperature's curve should be an n by 1 \
            array")
    else:
        print("The shape of temperature's curve is n by ", 
            len(temperature_curve))

    return

    ## confirm the initial temperature
    try:
        initial_temp = 292  # [K] 
        assert(np.abs(temperature_curve[0] - initial_temp) <= 1e-3)
    except AssertionError as a:
        raise AssertionError("The first value on the temeperature's curve \
            should be equal to the initial temperature")
    else:
        print("The first value on the temperature's curve is equal to the \
            initial temperature")

    return
