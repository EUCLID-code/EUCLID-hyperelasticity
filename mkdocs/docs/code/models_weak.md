#### `assembleB(data,ele,c):`

Define a matrix (B-matrix) such that the gradient of the test functions can be expressed as a product of the B-matrix and the nodal values of the test functions.

_Input Arguments_

- `data`
   
- `ele` - element number
   
- `c` - see `config`
   
_Output Arguments_

- `B_element` - B-matrix
   
---
   

#### `assembleFeatureDerivative(data,ele):`

Apply chain rule to calculate feature derivatives at element level.

_Input Arguments_

- `data`
   
- `ele` - element number
       
_Output Arguments_

- `d_features_dF_element` - derivative at element level
   
---
   

#### `assembleGlobalMatrix(matrix_global,matrix_element,connectivity,ele,c):`

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
   

#### `assembleGlobalVector(vector_global,vector_element,connectivity,ele,c):`

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
   

#### `computeInternalForceTheta(data,theta,c):`

Compute nodal internal forces for given material parameters.

_Input Arguments_

- `data`
   
- `theta` - material parameters
   
- `c` - see `config`
   
_Output Arguments_

- `f_int_nodes` - internal force at each node
   
---
   

#### `computeStrainEnergyTheta(data,theta):`

Compute strain energy density for given material parameters.

_Input Arguments_

- `data`
   
- `theta` - material parameters
   
_Output Arguments_

- `Psi` - strain energy density at each element
   
- `Psi_total` - total strain energy density
   
---
   

#### `computeWeakLHS(data,c):`

Compute left hand side of overdetermined linear system for given data.
Reaction forces are not consider yet.

_Input Arguments_

- `data`
   
- `c` - see `config`
   
_Output Arguments_

- `LHS` - left hand side of overdetermined linear system
   
---
   

#### `computeWeakOverdetermined(data,c):`

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
   

#### `considerReactionGlobal(data,LHS,c):`

Given overdetermined linear system, add information from global reaction forces.

_Input Arguments_

- `data`
   
- `LHS` - left hand side of symmetric linear system
   
- `c` - see `config`
   
_Output Arguments_

- `LHS` - updated left hand side of symmetric linear system
   
- `RHS` - right hand side of symmetric linear system
   
---
   

#### `extendLHSRHS(data,c,LHS,RHS):`

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
   

#### `zipper(matrix):`

Take a Nx2 matrix and assemble a 2Nx1 vector of the same type.

_Input Arguments_

- `matrix`
   
_Output Arguments_

- `vector`
   
---
   