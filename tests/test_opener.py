from pathlib import Path
import os.path

from chem461.pchem import Opener


def test_get_txt_data():
    """is the dataset an n by 2 array?"""
    abspath = os.path.abspath("data/exp42_dataset.csv")
    dnfn = Path(abspath)
    # dnfn = "../data/exp42_dataset.csv"
    testopener = Opener()

    # check that file exists
    try:
        file_path = dnfn.resolve(strict=True)
    except FileNotFoundError:
        print("The file path does not exists")
    else:
        print("The file path exists")
    
    # check shape of file
    try:
        ds, _ = testopener.get_txt_data(dnfn)
        assert(len(ds[0]) == 2)
    except AssertionError as a:
        raise AssertionError("The dataset should be an n by 2 array")
    else:
        print("The shape of ds is ", len(ds[:]), " by ds.shape ", len(ds[0]))

    try:
        ds = testopener.getdata(dnfn)
        # ignore the third empty column that is read by pandas
        ds_trunc = ds[:, 0:2]
        assert(len(ds_trunc[1]) == 2)
    except AssertionError as a:
        raise AssertionError("The dataset should be an n by 2 array")
    else:
        print("The shape of ds is ", ds_trunc.shape) 

    return
