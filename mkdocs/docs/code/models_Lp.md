#### `applyPenaltyLpIteration(datasets,LHS,RHS,c):`

Start with a small penalty parameter and solve the Lp-regularized optimization problem, increase the penalty parameter until physical requirements are fulfilled.

_Input Arguments_

- `datasets`
   
- `LHS` - left hand side of symmetric linear system
   
- `RHS` - right hand side of symmetric linear system
   
- `c` - see `config`
   
_Output Arguments_

- `theta` - material parameters
   
---
   

#### `applyPenaltyLpRandomStart(datasets,LHS,RHS,c):`

Solve the Lp-regularized optimization problem with random initial guesses.

_Input Arguments_

- `datasets`
   
- `LHS` - left hand side of symmetric linear system
   
- `RHS` - right hand side of symmetric linear system
   
- `c` - see `config`
   
_Output Arguments_

- `theta` - material parameters
   
- `at_least_one_converged` - information about the convergence
   
---
   

#### `applyPenaltyLpThreshold(datasets,LHS,RHS,theta,c):`

Solve the Lp-regularized optimization problem for a given initial guess using fixed-point iteration.
During the iterative process, the threshold algorithm is applied.

_Input Arguments_

- `datasets`
   
- `LHS` - left hand side of symmetric linear system
   
- `RHS` - right hand side of symmetric linear system
   
- `theta` - initial guess of material parameters
   
- `c` - see `config`
   
_Output Arguments_

- `theta` - material parameters
   
- `converged` - information about the convergence
   
---
   

#### `checkLocalMinimumLp(datasets,theta,c):`

Check if the provided solution is a local minimum of the cost function.

_Input Arguments_

- `datasets`
   
- `theta` - material parameters
   
- `c` - see `config`
   
_Output Arguments_

- _none_
   
---
   

#### `computeCostLp(datasets,theta,c):`

Compute cost function for L2 minimization with Lp-regularization.

_Input Arguments_

- `datasets`
   
- `theta` - material parameters
   
- `c` - see `config`
   
_Output Arguments_

- `Cost_weak` - cost corresponding to the weak linear momentum balance
   
- `Cost_penaltyLp` - cost corresponding to the Lp-regularization
   
- `Cost_total` - total cost
   
---
   

#### `saveResultsLp(theta,c,counter_load=None):`

Save the results and the chosen parameters in a text file.

_Input Arguments_

- `theta` - material parameters
   
- `c` - see `config`
   
- `counter_load` - current load step
   
_Output Arguments_

- _none_
   
---
   