from .__importList__ import *

def computeFeatures(I1, I2, I3):
    """
    Compute the features dependent on the right Cauchy-Green strain invariants.
    Note that the features only depend on I1 and I3.
    
    _Input Arguments_
    
    - `I1` - 1st invariant
    
    - `I2` - 2nd invariant

    - `I3` - 3rd invariant
    
    _Output Arguments_
    
    - `x` - features
    
    ---
    
    """
    #Generalized Mooney-Rivlin.
    #The Gent-Thomas model cannot be represented by the generalized
    #Mooney-Rivlin model. An additional feature has to be added.
    considerGentThomas = True
    #Polynomial terms degree.
    N = 7
    #Volumetric terms degree.
    M = 7
#    print("Generate a generalized Mooney-Rivlin feature library with..." )
#    print("N: " + str(N))
#    print("M: " + str(M))
    K1 = I1 * torch.pow(I3,-1/3) - 3.0
    K2 = (I1 + I3 - 1) * torch.pow(I3,-2/3) - 3.0
    J = torch.sqrt(I3)
    #Calculate the number of features.
    numFeatures = 0
    #Polynomial terms (dependent on K1 and K2).
    for n in range(N):
        numFeatures += n + 2
#        numFeatures += 1 # REMOVE K2
    #Volumetric terms (dependent on J).
    for m in range(M):
        numFeatures += 1
    #Additional Gent-Thomas feature.
    if considerGentThomas:
        numFeatures += 1
    #Calculate the features.
    x = torch.zeros(I1.shape[0],numFeatures)
    i=-1;
    #Polynomial terms (dependent on K1 and K2).
    for p in range(1,N+1):
        for q in range(p+1):
            i+=1; x[:,i:(i+1)] = K1**(p-q) * K2**q
#        for q in range(1): # REMOVE K2
#            i+=1; x[:,i:(i+1)] = K1**(p-q) # REMOVE K2
    #Volumetric terms (dependent on J):
    for m in range(1,M+1):
        i+=1; x[:,i:(i+1)] = (J-1)**(2*m)
    #Additional Gent-Thomas feature.
    if considerGentThomas:
        i+=1; x[:,i:(i+1)] = torch.log((K2+3.0)/3.0)

    return x

def getNumberOfFeatures():
    """
    Compute number of features.
    
    _Input Arguments_
    
    - _none_
    
    _Output Arguments_
    
    - `features.shape[1]` - number of features
    
    ---
    
    """
    features = computeFeatures(torch.zeros(1,1),torch.zeros(1,1),torch.zeros(1,1))
    return features.shape[1]








