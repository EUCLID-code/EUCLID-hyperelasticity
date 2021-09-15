#=====================================================================
# INITIALIZATIONS:
#=====================================================================
from .__importList__ import *
from .utilities import *
from .features import *
#=====================================================================
# SPARSE REGRESSION, WEAK FORM, LOCAL FE INTERPOLATION:
#=====================================================================
def extendLHSRHS(data,c,LHS,RHS):
    """
    Given symmetric linear system, add information from additional data.
    
    _Input Arguments_
    
    - `data`

    - `c` - see `config`

    - `LHS` - left hand side of symmetric linear system
    
    - `RHS` - right hand side of symmetric linear system
    
    _Output Arguments_
    
    - `LHS` - updated left hand side of symmetric linear system

    - `RHS` - updated right hand side of symmetric linear system
    
    ---
    
    """
    print('\n-----------------------------------------------------')
    print('Consider new load step in LHS and RHS.')
    print('Balance factor: ', c.balance)
    LHS_weak = computeWeakLHS(data,c)
    #Assumption: force is measured globally at the boundary.
    LHS_step , RHS_step = considerReactionGlobal(data,LHS_weak,c)
    LHS += LHS_step
    RHS += RHS_step
    print('-----------------------------------------------------\n')
    return LHS , RHS

def computeInternalForceTheta(data,theta,c):
    """
    Compute nodal internal forces for given material parameters.
    
    _Input Arguments_
    
    - `data`

    - `theta` - material parameters
    
    - `c` - see `config`

    _Output Arguments_
    
    - `f_int_nodes` - internal force at each node
        
    ---
    
    """
    #Computes the internal force by integrating the product of the 1st Piola Kirchhoff stress tensor
    #and the B-matrix over each element.
    f_int_nodes = np.zeros(c.dim*data.numNodes) #allocate memory
    for ele in range(data.numElements):                
        d_features_dF_element = assembleFeatureDerivative(data,ele)
        P_element = np.transpose(d_features_dF_element).dot(theta)
        B_element = assembleB(data,ele,c)
        f_int_element = np.transpose(B_element).dot(P_element) * data.qpWeights[ele]
        f_int_nodes = assembleGlobalVector(f_int_nodes,f_int_element,data.connectivity,ele,c)
    return f_int_nodes

def computeStrainEnergyTheta(data,theta):
    """
    Compute strain energy density for given material parameters.
    
    _Input Arguments_
    
    - `data`

    - `theta` - material parameters
    
    _Output Arguments_
    
    - `Psi` - strain energy density at each element
    
    - `Psi_total` - total strain energy density
    
    ---
    
    """
    features = data.featureSet.features
    Psi = np.zeros(data.numElements)
    Psi_total = 0.0
    for ele in range(data.numElements):
        Psi[ele] = features[ele,:].dot(theta)
        Psi_total += Psi[ele] * data.qpWeights[ele]
    return Psi , Psi_total

def computeWeakOverdetermined(data,c):
    """
    Compute left and right hand side of the internal-external force
    balance (weak formulation) at all free and fixed degrees of freedom.
    
    _Input Arguments_
    
    - `data`

    - `c` - see `config`
    
    _Output Arguments_
    

    - `LHS_free` - left hand side of non-symmetric linear system

    - `RHS_free` - right hand side of non-symmetric linear system

    - `LHS_fix` - left hand side of non-symmetric linear system

    - `RHS_fix` - right hand side of non-symmetric linear system
    
    ---
    
    """
    print('\n-----------------------------------------------------')
    print('Assemble overdetermined system of equations.')
    LHS = computeWeakLHS(data,c)
    non_dirichlet_nodes = ~(zipper(data.dirichlet_nodes)) # ~ can cause errors -> use np.logical_not
    LHS_free = LHS[non_dirichlet_nodes,:]
    RHS_free = np.zeros(LHS_free.shape[0])
    LHS_fix = np.zeros([len(data.reactions),LHS_free.shape[1]])
    RHS_fix = np.zeros(len(data.reactions))
    reaction_counter = -1
    for reaction in data.reactions:
        reaction_counter += 1
        force = reaction.force
        dofs = zipper(reaction.dofs)
        one_vector = np.ones(LHS[dofs,:].shape[0])
        LHS_fix[reaction_counter] = np.transpose(one_vector).dot(LHS[dofs,:])
        RHS_fix[reaction_counter] = force
    print('-----------------------------------------------------\n')
    return LHS_free , RHS_free , LHS_fix , RHS_fix

def computeWeakLHS(data,c):
    """
    Compute left hand side of overdetermined linear system for given data.
    Reaction forces are not consider yet.
    
    _Input Arguments_
    
    - `data`

    - `c` - see `config`
    
    _Output Arguments_
    
    - `LHS` - left hand side of overdetermined linear system
    
    ---
    
    """
    numFeatures = data.featureSet.features.shape[1]
    LHS = np.zeros((c.dim*data.numNodes,numFeatures)) #allocate memory
    for ele in range(data.numElements):                
        d_features_dF_element = assembleFeatureDerivative(data,ele)
        B_element = assembleB(data,ele,c)
        LHS_element = np.transpose(B_element).dot(np.transpose(d_features_dF_element)) * data.qpWeights[ele]
        LHS = assembleGlobalMatrix(LHS,LHS_element,data.connectivity,ele,c)
    return LHS

def considerReactionGlobal(data,LHS,c):
    """
    Given overdetermined linear system, add information from global reaction forces.
    
    _Input Arguments_
    
    - `data`

    - `LHS` - left hand side of symmetric linear system
    
    - `c` - see `config`
        
    _Output Arguments_
    
    - `LHS` - updated left hand side of symmetric linear system

    - `RHS` - right hand side of symmetric linear system
    
    ---
    
    """
    #Considers Reaction forces to be known globally.
    non_dirichlet_nodes = ~(zipper(data.dirichlet_nodes))
    LHS_bulk = LHS[non_dirichlet_nodes,:]
    #The overdetermined linear system can be solved in a least-squares sense.
    LHS_bulk = 2*np.transpose(LHS_bulk).dot(LHS_bulk)
    LHS_reaction = np.zeros_like(LHS_bulk)
    RHS_reaction = np.zeros(LHS_bulk.shape[0])
    for reaction in data.reactions:
        force = reaction.force
        dofs = zipper(reaction.dofs)
        one_vector = np.ones(LHS[dofs,:].shape[0])
        LHS_reaction += 2*np.outer( np.transpose(LHS[dofs,:]).dot(one_vector) , np.transpose(one_vector).dot(LHS[dofs,:]) )
        RHS_reaction += 2*np.transpose(LHS[dofs,:]).dot(one_vector) * force    
    LHS = LHS_bulk+c.balance*LHS_reaction
    RHS = c.balance*RHS_reaction
    return LHS , RHS

#=====================================================================
# HELP MATRICES:
#=====================================================================
def assembleFeatureDerivative(data,ele):
    """
    Apply chain rule to calculate feature derivatives at element level.
        
    _Input Arguments_
    
    - `data`

    - `ele` - element number
            
    _Output Arguments_
    
    - `d_features_dF_element` - derivative at element level
    
    ---
    
    """
    d_features_dI1_element = data.featureSet.d_features_dI1[ele,:]
    d_features_dI3_element = data.featureSet.d_features_dI3[ele,:]
    dI1dF_element = data.dI1dF[ele,:]
    dI3dF_element = data.dI3dF[ele,:]
    d_features_dF_element = np.outer(d_features_dI1_element,dI1dF_element)+np.outer(d_features_dI3_element,dI3dF_element)
    return d_features_dF_element

def assembleB(data,ele,c):
    """
    Define a matrix (B-matrix) such that the gradient of the test functions can be expressed as a product of the B-matrix and the nodal values of the test functions.
        
    _Input Arguments_
    
    - `data`

    - `ele` - element number
    
    - `c` - see `config`
        
    _Output Arguments_
    
    - `B_element` - B-matrix
    
    ---
    
    """
    B_element = np.zeros((c.dim*c.dim,c.dim*c.numNodesPerElement)) #allocate memory
    for a in range(c.numNodesPerElement):
        dN1 = data.gradNa[a][ele,0]
        dN2 = data.gradNa[a][ele,1]
#        B_element[:,2*a:2*a+2] = np.array([(dN1,0.0),(0.0,dN1),(dN2,0.0),(0.0,dN2)])
        #Note that the second and third row are swopped in order to ensure that the double contraction
        #of the 1st Piola-Kirchhoff stress tensor and the gradient of the test functions can be expressed
        #as a scalar product of their Voigt notation counterparts.
        B_element[:,2*a:2*a+2] = np.array([(dN1,0.0),(dN2,0.0),(0.0,dN1),(0.0,dN2)])
    return B_element

#=====================================================================
# ARRAY REARRANGEMENT:
#=====================================================================
def assembleGlobalVector(vector_global,vector_element,connectivity,ele,c):
    """
    Add a vector defined at an element to the corresponding position in the global vector.
    
    _Input Arguments_
    
    - `vector_global` - global vector

    - `vector_element` - vector at element level

    - `connectivity` - contains connectivity information of the elements

    - `ele` - element number

    - `c` - see `config`
    
    _Output Arguments_
    
    - `vector_global` - updated global vector
    
    ---
    
    """
    #Adds a vector defined at an element to the corresponding position in the global vector.
    for a in range(c.numNodesPerElement):
        vector_global[2*connectivity[a][ele]] += vector_element[2*a]
        vector_global[2*connectivity[a][ele]+1] += vector_element[2*a+1]
    return vector_global
    
def assembleGlobalMatrix(matrix_global,matrix_element,connectivity,ele,c):
    """
    Add a matrix defined at an element to the corresponding position in the global matrix.
    
    _Input Arguments_
    
    - `matrix_global` - global matrix

    - `matrix_element` - matrix at element level

    - `connectivity` - contains connectivity information of the elements

    - `ele` - element number

    - `c` - see `config`
    
    _Output Arguments_
    
    - `matrix_global` - updated global matrix
    
    ---
    
    """
    for a in range(c.numNodesPerElement):
        matrix_global[2*connectivity[a][ele],:] += matrix_element[2*a,:]
        matrix_global[2*connectivity[a][ele]+1,:] += matrix_element[2*a+1,:]
    return matrix_global

def zipper(matrix):
    """
    Take a Nx2 matrix and assemble a 2Nx1 vector of the same type.
    
    _Input Arguments_
    
    - `matrix`
    
    _Output Arguments_
    
    - `vector`
    
    ---
    
    """
    vector = np.zeros(2*matrix.shape[0], dtype=type(matrix[0,0]))
    for i in range(matrix.shape[0]):
        vector[2*i] = matrix[i,0:1]
        vector[2*i+1] = matrix[i,1:2]
    return vector