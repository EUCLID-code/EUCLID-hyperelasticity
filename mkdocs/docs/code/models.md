#### `applyThreshold(LHS,RHS,theta,c):`

Apply the threshold algorithm to the given linear system to enforce sparsity on the solution vector theta.

_Input Arguments_

- `LHS` - left hand side of symmetric linear system
   
- `RHS` - right hand side of symmetric linear system
   
- `theta` - material parameters
   
- `c` - see `config`
   
_Output Arguments_

- `theta` - material parameters
   
---
   

#### `computeDoubleContraction(A,B):`

Double contraction in Voigt notation.

_Input Arguments_

- `A` - vector containing components of a 2x2 matrix ([A_11 A_12 A_21 A_22])
   
- `B` - vector containing components of a 2x2 matrix ([B_11 B_12 B_21 B_22])
   
_Output Arguments_

- `C` - scalar value of the double contraction
   
---
   

#### `computeLpNorm(vector,p):`

Compute the Lp-regularization term for a given vector.

_Input Arguments_

- `vector`
   
- `p` - type of the Lp-regularization term
   
_Output Arguments_

- `norm`
   
---
   