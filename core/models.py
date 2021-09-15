#=====================================================================
# INITIALIZATIONS:
#=====================================================================
from .__importList__ import *
from .utilities import *
from .features import *
#=====================================================================
# THRESHOLD:
#=====================================================================
def applyThreshold(LHS,RHS,theta,c):
    """
    Apply the threshold algorithm to the given linear system to enforce sparsity on the solution vector theta.
    
    _Input Arguments_
    
    - `LHS` - left hand side of symmetric linear system
    
    - `RHS` - right hand side of symmetric linear system
    
    - `theta` - material parameters
    
    - `c` - see `config`
    
    _Output Arguments_
    
    - `theta` - material parameters
    
    ---
    
    """
    print('\n-----------------------------------------------------')
    print('Apply threshold: ', c.threshold)
    print()
    for i in range(len(theta)):
        #Find parameters less than threshold.
        theta_less = np.abs(theta) < c.threshold
        if all(theta[theta_less] == 0.0):
            print("Converged after iteration:", i)
            #Convergence of this method is guaranteed.
            break
        #Find parameters greater than threshold.
        theta_greater = np.abs(theta) >= c.threshold
        theta[theta_less] = 0.0
        LHS_temp = LHS[theta_greater,:]
        LHS_temp = LHS_temp[:,theta_greater]
        theta[theta_greater] = np.linalg.solve(LHS_temp,RHS[theta_greater])
    print("Solution:")
    print(theta)
    print('-----------------------------------------------------\n')
    return theta

#=====================================================================
# MISCELLANEOUS:
#=====================================================================
def computeLpNorm(vector,p):
    """
    Compute the Lp-regularization term for a given vector.
    
    _Input Arguments_
    
    - `vector`
    
    - `p` - type of the Lp-regularization term
    
    _Output Arguments_
    
    - `norm`
    
    ---
    
    """
    norm = np.power(np.sum(np.power(np.abs(vector),p)),1/p)
    return norm

def computeDoubleContraction(A,B):
    """
    Double contraction in Voigt notation.
    
    _Input Arguments_
    
    - `A` - vector containing components of a 2x2 matrix ([A_11 A_12 A_21 A_22])
    
    - `B` - vector containing components of a 2x2 matrix ([B_11 B_12 B_21 B_22])
    
    _Output Arguments_
    
    - `C` - scalar value of the double contraction
    
    ---
    
    """
    C = A[0]*B[0] + A[1]*B[2] + A[2]*B[1] + A[3]*B[3]
    return C
