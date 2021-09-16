# Example
It follows a step-by-step description of the implemented algorithm for an exemplary dataset.
In particular, EUCLID (Efficient Unsupervised Constitutive Law Identification & Discovery) is applied to displacement and reaction force data that were generated based on the Neo-Hookean material model NH2.
The goal is to use the data to discover the material model, without knowing its specific mathematical form a priori.

## Data and Parameters
The input data for EUCLID and the parameters for the optimization process are defined in `config.py`.
First, the material model based on which the data were generated, the mesh, the load steps and the noise level are specified.
```Python
str_model = 'NeoHookeanJ2'
str_mesh = 'plate_hole_1k_'
loadsteps = [10,20,30,40]
noiseLevelData = '0'
noiseLevel = 0.0001
```
Here, `noiseLevelData` is set to zero meaning that we want to load data that has not yet been corrupted by noise
and `noiseLevel` is set to be greater than zero as artificial noise should be added to that data.
We have to make sure that `femDataPath` describes the correct path pointing to the finite element data.
Higher resolution data can be downloaded from the <a href="https://www.research-collection.ethz.ch/handle/20.500.11850/505693" target="_blank">ETH Research Collection</a>.

After specifying which data to use, the hyperparameters for EUCLID are chosen.
```Python
balance = 100
penaltyLp = 1e-4
p = 1.0/4.0
factorIncrements = 5
numGuesses = 1
threshold = 1e-2
```
Setting `balance = 100` increases the influence of boundary data on the objective function
and `penaltyLp` defines the initial value of the penalty factor in front of the Lp-regularization term added to the objective function.
The latter parameter will be successively increased by the factor `factorIncrements` until a physically admissible material model is found.
For this example, we choose only one initial guess for the minimization probem
and the threshold below which material parameters should be discarded is defined as `threshold = 1e-2`. 

## Discovery
Executing `main_Lp.py` starts the material model discovery process.
The parameters defined in `config.py` are loaded as `import config as c`.
A loop over all load steps is defined and for each load step, the finite element data is loaded and perturbed by the specified noise.
```Python
for loadstep in c.loadsteps:
	data = loadFemData(c.femDataPath+'/'+str(loadstep), AD=True, noiseLevel = c.noiseLevel, noiseType = 'displacement')
	datasets.append(data)
```
The left hand side `LHS` and right hand side `RHS` of the system of linear equations, which are obtained from differentiating the minimization function (without Lp-regularization), are assembled.
Finally, the material model is discovered from the data by executing:
```Python
theta = applyPenaltyLpIteration(datasets,LHS,RHS,c)
```

## Results
In `applyPenaltyLpIteration` the fixed-point iteration is applied for the initial penalty factor in front of the sparsity promoting Lp-regularization term.
As the initial value has been choosen negligibly small, the obtained solution vector for the material parameters `theta` is dense, i.e., it contains a large number of nonzero entries.
```Matlab
-----------------------------------------------------
Apply Lp-norm penalty.
Lp-norm penalty factor:  0.0001
Lp-norm: p= 0.25
Number of initial guesses: 1

Solution with the lowest cost:
[ 6.08591003e-01 -1.26556172e+01 -1.50304242e+00  3.47039792e+00
  0.00000000e+00 -8.40926870e-01  4.26303110e+00 -6.44342608e+00
  2.68558116e+00  1.91942291e-01 -5.80249192e-01  0.00000000e+00
  1.13330467e+00 -7.01748020e-01 -8.85663806e-03  0.00000000e+00
  7.39918129e-02 -4.10729435e-02 -1.23457739e-01  9.60345278e-02
  0.00000000e+00  9.97033677e-04  0.00000000e+00 -1.09769271e-02
  1.66016891e-02 -3.06059710e-03 -3.42983131e-03 -1.06609625e-05
  7.32603777e-05 -1.98700714e-04  0.00000000e+00  1.03872040e-03
 -1.85894059e-03  1.16813894e-03 -2.13678851e-04  1.50478155e+00
  0.00000000e+00 -1.44416532e-01  7.52157815e-01 -1.03935826e+00
  5.20008688e-01 -8.75636549e-02  3.75995295e+01]
-----------------------------------------------------
```
Such a solution is not desired as the resulting material model would have many material parameters and would be physically uninterpretable.
Further, the resulting model shown above is not physically admissible,
which can be proven by calculating the strain energy density along specific deformation paths (see `model_requirements.py`).
Hence, in `applyPenaltyLpIteration` the penalty factor `penaltyLp` is successively increased until the physical requirements are fulfilled.
This procedure finally results in:
```Matlab
-----------------------------------------------------
Apply Lp-norm penalty.
Lp-norm penalty factor:  0.0125
Lp-norm: p= 0.25
Number of initial guesses: 1

Solution with the lowest cost:
[ 4.97871855e-01 -5.74695560e-05  0.00000000e+00  0.00000000e+00
  0.00000000e+00  7.35505936e-05  0.00000000e+00  0.00000000e+00
  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
 -7.84065373e-06  0.00000000e+00  0.00000000e+00  0.00000000e+00
  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.49898747e+00
  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
  0.00000000e+00  0.00000000e+00  0.00000000e+00]
-----------------------------------------------------
```
It can be seen that the material parameters have been shrunk, and the physical requirements are now fulfilled.
After a final postprocessing step, in which parameters with absolute value below a `threshold` are set to zero, a sparse material parameter vector is found.
```Matlab
-----------------------------------------------------
Apply threshold:  0.01

Converged after iteration: 1
Solution:
[0.50052861 0.         0.         0.         0.         0.
 0.         0.         0.         0.         0.         0.
 0.         0.         0.         0.         0.         0.
 0.         0.         0.         0.         0.         0.
 0.         0.         0.         0.         0.         0.
 0.         0.         0.         0.         0.         1.49892849
 0.         0.         0.         0.         0.         0.
 0.        ]
-----------------------------------------------------
```
The resulting material model has hence a small amount of material parameters and is physically interpretable.
It is observed that EUCLID identifies the correct features of the Neo-Hookean material model.
The material parameter values vary only slightly compared to the true parameters as a consequence of the noise added to the displacement data.



