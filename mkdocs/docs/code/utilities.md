#### `computeCauchyGreenStrain(F):`

Compute right Cauchy-Green strain tensor from deformation gradient.

_Input Arguments_

- `F` - deformation gradient in Voigt notation
   
_Output Arguments_

- `C` - Cauchy-Green strain tensor in Voigt notation
   
---
   

#### `computeJacobian(F):`

Compute Jacobian from deformation gradient.

_Input Arguments_

- `F` - deformation gradient in Voigt notation
   
_Output Arguments_

- `J` - Jacobian
   
---
   

#### `computeStrainInvariantDerivatives(F,i,secondDerivative=False):`

Compute derivatives of the invariants of the Cauchy-Green strain tensor with respect to the deformation gradient.
Plane strain is assumed.

_Input Arguments_

- `F` - deformation gradient in Voigt notation
   
- `i` - specify the invariant that should be differentiated 
   
- `secondDerivative` - specify if second derivative should be computed 
   
_Output Arguments_

- `dIdF` - derivative (note that the size of `dIdF` depends on the choice of `secondDerivative`)
   
---
   

#### `computeStrainInvariants(C):`

Compute invariants of the Cauchy-Green strain tensor.
Plane strain is assumed.

_Input Arguments_

- `C` - Cauchy-Green strain tensor in Voigt notation
   
_Output Arguments_

- `I1` - 1st invariant
   
- `I2` - 2nd invariant
   
- `I3` - 3rd invariant
   
---
   