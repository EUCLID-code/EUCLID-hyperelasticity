#=====================================================================
# INITIALIZATIONS:
#=====================================================================
#sys and core
import sys
sys.path.insert(0, "../")
from core import *
#config
import config as c
#CUDA
initCUDA(c.cuda)

datasets = []
counter_load = -1
for loadstep in c.loadsteps:
    counter_load +=1
#=====================================================================
# DATA:
#=====================================================================
    data = loadFemData(c.femDataPath+'/'+str(loadstep), AD=True, noiseLevel = c.noiseLevel, noiseType = 'displacement')
    datasets.append(data)
    data.convertToNumpy()
#=====================================================================
# LHS & RHS OF THE WEAK FORM:
#=====================================================================
    numFeatures = getNumberOfFeatures()
    if counter_load == 0:
        LHS = np.zeros([numFeatures,numFeatures])
        RHS = np.zeros([numFeatures])
    LHS , RHS = extendLHSRHS(data,c,LHS,RHS)
#=====================================================================
# Lp-NORM PENALTY:
#=====================================================================
    if loadstep == c.loadsteps[-1]:
        theta = applyPenaltyLpIteration(datasets,LHS,RHS,c)
#=====================================================================
# PRINT RESULTS:
#=====================================================================
        saveResultsLp(theta,c,counter_load=counter_load)
    c.penaltyLp = c.penaltyLp_init

