#==================================================================================================
#                     Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

import numpy as np
import CyURT as urt
from timeit import default_timer as timer 

if __name__ == "__main__":

    sizes = [100,150,200,250,300,350,400,450,500,1000,1500,2000,2500,3000,3500,4000,4500,5000]

    for i in range(len(sizes)):

        # generating Wiener process
        data = np.cumsum(np.random.normal(size=sizes[i]))
        # uncomment this line and comment the one above to switch to single precision
        #data = np.cumsum(np.random.normal(size=sizes[i])).astype(np.float32)

        if sizes[i] < 1000: niter = 10000
        else: niter = 1000
        
        start = timer()
        for k in range(niter):
            test = urt.ADF_d(data, method='AIC')
            # uncomment this line and comment the one above to switch to single precision
            #test = urt.ADF_f(data, method='AIC')
            test.statistic()
        end = timer()

        print '{:8d}'.format(sizes[i]), '{:8.1f}'.format(end - start)

