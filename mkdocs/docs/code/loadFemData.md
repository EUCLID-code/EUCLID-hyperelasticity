#### `loadFemData(path, AD = True, noiseLevel = 0., noiseType = 'displacement', denoisedDisplacements = None):`

Load finite element data and add noise (optional).
Note that the loaded finite element data might already be perturbed by noise.
In that case, adding additional noise is not necessary.

_Input Arguments_

- `path` - path to finite element data
   
- `AD` - specify if automatic differention is needed
   
- `noiseLevel` - noise level
   
- `noiseType` - specify whether noise should be added to displacements or strains
   
- `denoisedDisplacements` - pass denoised displacement data if available
   
_Output Arguments_

- `dataset` - finite element dataset
   
---
   