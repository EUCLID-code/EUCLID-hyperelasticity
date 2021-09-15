#=====================================================================
# INITIALIZATIONS:
#=====================================================================
from .__importList__ import *
from .utilities import *
from .features import *
#=====================================================================
# MODELS LIST:
#=====================================================================
from .models_weak import *
#=====================================================================
# ENERGY REQUIREMENTS:
#=====================================================================
def checkEnergyRequirementsRigorous(datasets,theta):
    """
    Check if physical requirements are fulfilled.
    
    _Input Arguments_
    
    - `datasets`

    - `theta` - material parameters
        
    _Output Arguments_
    
    - `fulfillEnergyRequirements` - specify if the physical requirements are fulfilled
    
    ---
    
    """
    print('\n-----------------------------------------------------')
    print('Check if the energy requirements are fulfilled.')
    fulfillEnergyRequirements = True
    #In order to obtain zero strain energy for undeformed materials, the strain
    #energy is normalized, i.e. a constant value is added to the strain energy.
    F = torch.zeros((1,4))
    F[0,0] = 1.0
    F[0,3] = 1.0
    C = computeCauchyGreenStrain(F)
    I1, I2, I3 = computeStrainInvariants(C)
    Q = computeFeatures(I1, I2, I3).numpy()[0,:]
    Psi = Q.dot(theta)
    print("Constant that needs to be added to the energy: ", -Psi)
    #The energy must be positive at every point in the data set.
    counter = 0
    for data in datasets:
        counter += 1
        Psi , Psi_total = computeStrainEnergyTheta(data,theta)
        Psi_min = np.min(Psi)
        if Psi_min < 0.0:
            print("Energy is negative at at least one point in the " + str(counter) + ". data set:")
            print(Psi_min)
            print('The energy requirements are not fulfilled.')
            print('-----------------------------------------------------\n')
            fulfillEnergyRequirements = False
            return fulfillEnergyRequirements
        else:
            print("Energy is positive in the whole " + str(counter) + ". data set.")
    #The energy must go to infinite for large deformations.
    deformationPathList = ['tension','bitension','compression','bicompression','simpleShear','pureShear']
    pathx = getPathX()
    print('Number of large deformation checks: ', len(deformationPathList)*len(pathx))
    for deformationPath in deformationPathList:
        Psi_temp = 0
        for x in pathx:
            Q = computeFeaturesGivenDeformation(deformationPath,x)
            Psi = Q.dot(theta)
#            print("Deformation (" + deformationPath + "): The energy is: ", Psi)
            if np.isnan(Psi) or np.isinf(Psi):
                break
            if Psi < (Psi_temp - 1e-10):
                print('The energy requirements are not fulfilled.')
                print('-----------------------------------------------------\n')
                fulfillEnergyRequirements = False
                return fulfillEnergyRequirements
            Psi_temp = np.copy(Psi)
    if fulfillEnergyRequirements:
        print('The energy requirements are fulfilled.')
    print('-----------------------------------------------------\n')
    return fulfillEnergyRequirements

#=====================================================================
# ARBITRARILY GENERATED DEFORMATIONS:
#=====================================================================
def computeFeaturesGivenDeformation(deformationPath,x):
    """
    Compute features for given deformation.
    
    _Input Arguments_
    
    - `deformationPath` - type of deformation

    - `x` - scalar deformation measure
        
    _Output Arguments_
    
    - `Q` - features
    
    ---
    
    """
    F = torch.zeros((1,4))
    if deformationPath == None:
        F[:,0] = 1.0
        F[:,3] = 1.0        
    elif deformationPath == 'tension':
        F[:,0] = x
        F[:,3] = 1.0
    elif deformationPath == 'bitension':
        F[:,0] = x
        F[:,3] = x
    elif deformationPath == 'compression':
        F[:,0] = 1/x
        F[:,3] = 1.0
    elif deformationPath == 'bicompression':
        F[:,0] = 1/x
        F[:,3] = 1/x
    elif deformationPath == 'simpleShear':
        F[:,0] = 1.0
        F[:,1] = x - 1.0
        F[:,3] = 1.0
    elif deformationPath == 'pureShear':
        F[:,0] = x
        F[:,3] = 1/x 
    C = computeCauchyGreenStrain(F)
    I1, I2, I3 = computeStrainInvariants(C)
    Q = computeFeatures(I1, I2, I3).numpy()[0,:]
    return Q

def getPathX():
    """
    Generate a scalar deformation path.
    
    _Input Arguments_
    
    - _none_
    
    _Output Arguments_
    
    - `pathx` - scalar deformation path
        
    ---
    
    """
    # a value of 2 results in NAN for shear
    pathx = [1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,\
             1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,\
             2.1,3,4,5,6,7,8,9,\
             10,20,30,40,50,60,70,80,90,\
             1e2,2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2,\
             1e3,2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,\
             1e4,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,\
             1e5,2e5,5e5,1e6,2e6,5e6,1e7,2e7,5e7,1e8,2e8,5e8,1e9]
    return pathx


