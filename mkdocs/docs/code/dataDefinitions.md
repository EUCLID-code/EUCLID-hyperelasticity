# `class Reaction:`

_Attributes_

- `force` - scalar force value

- `dofs` - degrees of freedom corresponding to the reaction force

_Methods_

- `__init__(...):` - Generate `Reaction` object.

---
   

# `class FeatureSet:`

_Attributes_

- `features`

- `d_features_dI1` - derivative

- `d_features_dI2` - derivative

- `d_features_dI3` - derivative


_Methods_

- `__init__(...):` - Generate `FeatureSet` object.

- `convertTensorToNumpy(...):` - Convert to numpy.

- `convertToNumpy(...):` - Convert to numpy.


---
   

# `class FemDataset:`

_Attributes_

- `path`

- `numNodes` - number of nodes

- `x_nodes` - node reference positions

- `u_nodes` - node displacements

- `dirichlet_nodes` - Dirichlet nodes

- `reactions`

- `numElements` - number of finite elements

- `qpWeights` - weights of the quadrature points

- `connectivity` - contains connectivity information of the elements

- `gradNa` - shape function gradient

- `F` - deformation gradient in Voigt notation

- `J` - Jacobian

- `C` - Cauchy-Green strain tensor in Voigt notation

- `I1` - 1st invariant

- `I2` - 2nd invariant

- `I3` - 3rd invariant

- `dI1dF` - derivative

- `dI2dF` - derivative

- `dI3dF` - derivative

- `featureSet`

_Methods_

- `__init__(...):` - Generate `FemDataset` object.

- `convertTensorToNumpy(...):` - Convert to numpy.

- `convertToNumpy(...):` - Convert to numpy.

---
   