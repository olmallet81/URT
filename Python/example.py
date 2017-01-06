#==================================================================================================
#                     Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

import numpy as np
import pandas as pd

import CyURT as urt
 
if __name__ == "__main__":

    y = pd.read_csv('../data/y.csv', sep=',', header=None)
    x = pd.read_csv('../data/x.csv', sep=',', header=None)

    yd = np.asarray(y).reshape(y.size)
    yf = yd.astype(np.float32)
    xd = np.asarray(x, order='F')
    xf = xd.astype(np.float32)

    # running OLS regression as in ./examples/example1.cpp using double precision type
    fit = urt.OLS_d(yd, xd, True)
    fit.show()

    # running OLS regression as in ./examples/example1.cpp using single precision type
    fit = urt.OLS_f(yf, xf, True)
    fit.show()

    # running first ADF test as in ./examples/example2.cpp using double precision type
    test = urt.ADF_d(yd, lags=10, trend='ct')
    test.show()

    # running second ADF test as in ./examples/example2.cpp using double precision type
    test.method = 'AIC'
    test.bootstrap = True
    test.niter = 10000
    test.show()
