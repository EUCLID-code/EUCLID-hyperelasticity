from .__importList__ import *
from .utilities import *

class Reference:
    def __init__(self):
        """
        Generate `Reference` object.
    
        """
        self.F = torch.tensor([[1.,0.,0.,1.]])
        self.J = computeJacobian(self.F)
        self.C = computeCauchyGreenStrain(self.F)
        self.I1, self.I2, self.I3 = computeStrainInvariants(self.C)
        self.I1.requires_grad = True
        self.I2.requires_grad = True
        self.I3.requires_grad = True
        self.dI1dF = computeStrainInvariantDerivatives(self.F,1)
        self.dI2dF = computeStrainInvariantDerivatives(self.F,2)
        self.dI3dF = computeStrainInvariantDerivatives(self.F,3)
        self.W = torch.tensor([[0.]])
        self.P = torch.tensor([[0.,0.,0.,0.]])

class Reaction:
    def __init__(self, dofs, force):
        """
        Generate `Reaction` object.
        
        """
        #dofs is a boolean matrix of size: numNodes x 2
        #dofs[a,i] == True: means node 'a' is constrained in 'i'-th direction 
        #and contributes to this particular reaction force measurement
        self.force = force #scalar force measurement
        self.dofs = dofs #boolean: numNodes x 2

class FeatureSet:
    def __init__(self, features=None, d_features_dI1=None, d_features_dI2=None, d_features_dI3=None, dd_features_dI1dI1=None, dd_features_dI1dI3=None, dd_features_dI3dI1=None, dd_features_dI3dI3=None):
        """
        Generate `FeatureSet` object.
        
        """
        #All matrices have rows = numElements, cols=numFeatures
        self.features = features 
        self.d_features_dI1 = d_features_dI1 
        self.d_features_dI2 = d_features_dI2 
        self.d_features_dI3 = d_features_dI3 
#        self.dd_features_dI1dI1 = dd_features_dI1dI1 
#        self.dd_features_dI1dI3 = dd_features_dI1dI3 
#        self.dd_features_dI3dI1 = dd_features_dI3dI1 
#        self.dd_features_dI3dI3 = dd_features_dI3dI3 
        
    def convertTensorToNumpy(self,T):
        """
        Convert to numpy.
        
        """
        T = T.cpu().detach().numpy()
        return T

    def convertToNumpy(self):
        """
        Convert to numpy.
        
        """
        self.features = self.convertTensorToNumpy(self.features)
        self.d_features_dI1 = self.convertTensorToNumpy(self.d_features_dI1)
        self.d_features_dI2 = self.convertTensorToNumpy(self.d_features_dI2)
        self.d_features_dI3 = self.convertTensorToNumpy(self.d_features_dI3)
#        self.dd_features_dI1dI1 = self.convertTensorToNumpy(self.dd_features_dI1dI1)
#        self.dd_features_dI1dI3 = self.convertTensorToNumpy(self.dd_features_dI1dI3)
#        self.dd_features_dI3dI1 = self.convertTensorToNumpy(self.dd_features_dI3dI1)
#        self.dd_features_dI3dI3 = self.convertTensorToNumpy(self.dd_features_dI3dI3)

class FemDataset:
    def __init__(self, path,
        x_nodes, u_nodes, dirichlet_nodes,
        reactions,
        connectivity, gradNa, qpWeights,
        F, J, C, 
        I1, I2, I3, 
        dI1dF, dI2dF, dI3dF,
        featureSet):
        """
        Generate `FemDataset` object.
        """

        self.path = path

        #---------------------------------------------------------------------------------
        #nodal data
        #---------------------------------------------------------------------------------
        self.numNodes = x_nodes.shape[0]

        self.x_nodes = x_nodes #node positions: numNodes x 2
        self.u_nodes = u_nodes #node displacements: numNodes x 2
        
        #dirichlet_nodes is a boolean matrix of size: numNodes x 2
        #dirichlet_nodes[a,i] == True: means node 'a' is constrained in 'i'-th direction 
        self.dirichlet_nodes = dirichlet_nodes #boolean: numNodes x 2
        
        #---------------------------------------------------------------------------------
        #reactions
        #---------------------------------------------------------------------------------
        self.reactions = reactions #all reactions as described in 'Reaction' class
        
        #---------------------------------------------------------------------------------
        #element data
        #---------------------------------------------------------------------------------
        self.numElements = qpWeights.shape[0]

        #qpWeights is the area (2D) or volume (3D) of each quadrature point.
        self.qpWeights = qpWeights # 1D matrix of length numElements

        #connectivity is a list of 3 1D-matrix. Each matrix contains integers.
        #Length of each matrix is equal to numElements. 
        #For 2D-triangular elements with a={0,1,2}:
        #connectivity[a][e]: is the index of the 'a'-th node of the 'e'-th element.
        #Note that node indices range from 0 to (numNodes-1).
        self.connectivity = connectivity 
        
        #gradNa is a list of 3 matrices, each of size: numElements x 2
        #For 2D-triangular elements with a={0,1,2} and i={0,1}:
        #gradNa[a][e,i]: is the 'i'-th component of the shape function gradient 
        #of the 'a'-th node of the 'e'-th element; evaluated at that element's 
        #single quadrature point
        self.gradNa = gradNa
        
        #---------------------------------------------------------------------------------
        #strain data
        #---------------------------------------------------------------------------------
        #The following are all 2D matrices with number of rows = numElements.
        #Each matrix represents certain strain-related quantity 
        #evaulated at the respecive element's quadrature points
        self.F = F #deformation gradient: numElements x 4
        self.J = J #determinant of F: numElements x 1
        self.C = C #right Cauchy-Green tensor C=F^T*F : numElements x 4
        self.I1 = I1 #first invariant of C : numElements x 1
        self.I2 = I2 #second invariant of C : numElements x 1
        self.I3 = I3 #third invariant of C : numElements x 1
        self.dI1dF = dI1dF #derivative of I1 with F : numElements x 4
        self.dI2dF = dI2dF #derivative of I2 with F : numElements x 4
        self.dI3dF = dI3dF #derivative of I3 with F : numElements x 4
        self.featureSet = featureSet #features defined in NN feature-library : numElements x (?)

    def convertTensorToNumpy(self,T):
        """
        Convert to numpy.
        
        """
        T = T.cpu().detach().numpy()
        return T

    def convertToNumpy(self):
        """
        Convert to numpy.
        
        """
        #nodal data
        self.x_nodes = self.convertTensorToNumpy(self.x_nodes)
        self.u_nodes = self.convertTensorToNumpy(self.u_nodes)
        self.dirichlet_nodes = self.convertTensorToNumpy(self.dirichlet_nodes)
        
        #reactions
        for i in range(len(self.reactions)):
            self.reactions[i].dofs = self.convertTensorToNumpy(self.reactions[i].dofs)
        
        #element data
        for i in range(len(self.connectivity)):
            self.connectivity[i] = self.convertTensorToNumpy(self.connectivity[i])
        for i in range(len(self.gradNa)):
            self.gradNa[i] = self.convertTensorToNumpy(self.gradNa[i])
        self.qpWeights = self.convertTensorToNumpy(self.qpWeights)
        
        #strain data
        self.F = self.convertTensorToNumpy(self.F)
        self.J = self.convertTensorToNumpy(self.J)
        self.C = self.convertTensorToNumpy(self.C)
        self.I1 = self.convertTensorToNumpy(self.I1)
        self.I2 = self.convertTensorToNumpy(self.I2)
        self.I3 = self.convertTensorToNumpy(self.I3)
        self.dI1dF = self.convertTensorToNumpy(self.dI1dF)
        self.dI2dF = self.convertTensorToNumpy(self.dI2dF)
        self.dI3dF = self.convertTensorToNumpy(self.dI3dF)
        self.featureSet.convertToNumpy()