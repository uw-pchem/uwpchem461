import numpy as np
import os.path
from pathlib import Path

from pchem import Opener, Analyse


def test_getbase():
    """
    is the baseline array the same in length as the dataset?
    """
    abspath = os.path.abspath("data/exp42_dataset.csv")
    dnfn = Path(abspath)
    opener = Opener()
    analyse = Analyse()
    ## check shape of file
    try:
        ds, _ = opener.get_txt_data(dnfn)
        dat = np.array(ds)
        hibar = 0.14
        baseline = analyse.getbase(dat, hibar)
        assert(len(baseline) == len(dat[:, 1]))
    except AssertionError as a:
        raise AssertionError("The baseline should be an n by 1 array")
    else:
        print("The shape of baseline is ", baseline.shape)

    return
