# If you need GPU-CUDA acceleration | -1: CPU, 0: CUDA0, 1: CUDA1
cuda = -1

# some essential definitions
dim = 2 
numNodesPerElement = 3
voigtMap = [[0,1],[2,3]]

#=====================================================================
# MODEL:
#=====================================================================
str_model = 'NeoHookeanJ2'
#str_model = 'NeoHookeanJ4'
#str_model = 'GentThomas'
#str_model = 'Isihara'
#str_model = 'HainesWilson'

#=====================================================================
# MESH:
#=====================================================================
str_mesh = 'plate_hole_1k_'
# str_mesh = 'plate_hole_60k_'

#=====================================================================
# NOISE LEVEL:
#=====================================================================
noiseLevelData = '0'
# noiseLevelData = '1e-4'
# noiseLevelData = '1e-3'

if noiseLevelData == '0':
    str_noise = 'FEM_data/' # no noise
elif noiseLevelData == '1e-4':
    str_noise = 'kr_denoisedData/denoised_60k_to_60k_noise=0.0001/'
elif noiseLevelData == '1e-3':
    str_noise = 'kr_denoisedData/denoised_60k_to_60k_noise=0.001/'

# additional noise
noiseLevel = 0.0001

#=====================================================================
# DATA PATH CLUSTER:
#=====================================================================
# femDataPath = '/cluster/scratch/mflaschel/data_08_25/' + str_noise + str_mesh + str_model

#=====================================================================
# DATA PATH LOCAL:
#=====================================================================
femDataPath = '../' + str_noise + str_mesh + str_model

#=====================================================================
# LOAD STEPS:
#=====================================================================
if (str_model[-1] == '2') or (str_model[-1] == '4'):
    loadsteps = [10,20,30,40]
else:
    loadsteps = [10,20,30,40,50,60,70,80]

#=====================================================================
# (HYPER-)PARAMETERS:
#=====================================================================
# hyperparameters for cost function
balance = 100

# Lp-regularization
penaltyLp = 1e-4 # initial lp-norm penalty parameter
penaltyLp_init = penaltyLp
p = 1.0/4.0 # type of lp-norm penalty
numIncrements = 5 # max. number of penalty parameter increments
factorIncrements = 5 # factor of penalty parameter increments

# fixed-point iteration solver
numGuesses = 1 # number of initial guesses for fixed-point iteration
numIterations = 200 # max. number of iterations for fixed-point iteration
lowestCost = -1
lowestCostGuessID = -1
threshold_iter = 1e-6 # absolute threshold parameter during fixed-point iteration

# final threshold
threshold = 1e-2 # absolute threshold parameter

# results
saveResultsName = str_model


