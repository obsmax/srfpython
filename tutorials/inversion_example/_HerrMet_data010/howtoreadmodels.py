import matplotlib.pyplot as plt
from srfpython.HerrMet.files import RunFile


fig = plt.figure()
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)


with RunFile('_HerrMet.run') as db:
    generator = db.get(  # < see also getzip and getpack
        llkmin=None,  # if you want to select based on the quality of the solution
        limit=1000,    # maximum number of "best" models you want
        step=None,    # if you want to take one model every "step" 
        algo=None)    # not important, leave to None
        
    for modelid, chainid, weight, llk, nlayer, (Z, VP, VS, RH), (W, T, M, F, DV) in generator:
    
        # Z is the array with the top depth of each layer
        # VP is the vp value in each layer in km/s
        # VS is the vs value in each layer in km/s        
        # RH is the density in g/cc
        
        #W,T,M,F,DV are the corresponding dispersion curves


        ax1.plot(VS, RH, 'o')
        ax2.plot(VP, RH, 'o')        
        ax3.plot(VS, VP, 'o')        
        
plt.show()
