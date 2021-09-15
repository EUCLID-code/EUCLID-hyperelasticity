#=====================================================================
# INITIALIZATIONS:
#=====================================================================
from .__importList__ import *
from .utilities import *
from .features import *
#=====================================================================
# MODELS LIST:
#=====================================================================
from .models_requirements import *
from .models_weak import *
from .models import *
#=====================================================================
# Lp-NORM PENALTY:
#=====================================================================
def applyPenaltyLpThreshold(datasets,LHS,RHS,theta,c):
    """
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
    
    """
    converged = False
    temp = np.zeros_like(theta)
    LHS_Lp = np.zeros_like(LHS)
    for i in range(c.numIterations):
        #Save parameters.
        temp = np.copy(theta)
        #Find parameters greater than threshold.
        theta_greater = np.abs(theta) >= c.threshold_iter
        #Compute the LHS of the fixed point mapping.
        penalty_weighted = np.zeros_like(theta)
        for idx in range(len(penalty_weighted)):
            if theta_greater[idx]:
                    penalty_weighted[idx] = np.power(np.abs(theta[idx]),c.p-2)
        LHS_Lp = LHS+c.p*c.penaltyLp*np.diag(penalty_weighted)
        #Apply the fixed-point iteration only to parameters greater than the
        #threshold.
        LHS_Lp = LHS_Lp[theta_greater,:]
        LHS_Lp = LHS_Lp[:,theta_greater]
        theta[theta_greater] = np.linalg.solve(LHS_Lp,RHS[theta_greater])
        for idx in range(len(theta)):
            if not(theta_greater[idx]):
                theta[idx] = 0
        if all( np.absolute(temp-theta) < 1e-3):
            converged = True
            # print("Converged after fixed-point iteration:", i)
            break
    return theta , converged

def applyPenaltyLpRandomStart(datasets,LHS,RHS,c):
    """
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
    
    """
    #Random initial guesses.
    print('\n-----------------------------------------------------')
    print('Apply Lp-norm penalty.')
    print('Lp-norm penalty factor: ', c.penaltyLp)
    print('Lp-norm: p=', c.p)
    print('Number of initial guesses:', c.numGuesses)
    print()
    theta = np.linalg.solve(LHS,RHS)
    at_least_one_converged = False
    for restarts in range(c.numGuesses):
        if restarts == 0:
            theta_rand = theta
        elif restarts < 50:
            theta_rand = np.zeros_like(RHS)
            for idx in range(len(theta_rand)):
                theta_rand[idx] = (2*np.random.rand()-1)
        else:
            theta_rand = np.zeros_like(RHS)
            for idx in range(len(theta_rand)):
                theta_rand[idx] = 10*(2*np.random.rand()-1)
        theta_rand , converged = applyPenaltyLpThreshold(datasets,LHS,RHS,theta_rand,c)
        if converged:
            Cost_weak , Cost_penaltyLp , Cost_total = computeCostLp(datasets,theta_rand,c)
            if not(at_least_one_converged):
                Cost_total_min = Cost_total
                theta = np.copy(theta_rand)
                c.lowestCost = Cost_total_min
                c.lowestCostGuessID = restarts
                at_least_one_converged = True
            if Cost_total < Cost_total_min:
                Cost_total_min = Cost_total
                theta = np.copy(theta_rand)
                c.lowestCost = Cost_total_min
                c.lowestCostGuessID = restarts
#            print(theta_rand)
#            print("Cost:", Cost_weak,"+" , Cost_penaltyLp,"+" , Cost_penaltyPsi,"=" , Cost_total)
        else:
            print("Fixed-point iteration not converged for guess number: ", restarts+1)
    print("Solution with the lowest cost:")
    print(theta)
    print('-----------------------------------------------------\n')
    return theta , at_least_one_converged

def applyPenaltyLpIteration(datasets,LHS,RHS,c):
    """
    Start with a small penalty parameter and solve the Lp-regularized optimization problem, increase the penalty parameter until physical requirements are fulfilled.
    
    _Input Arguments_
    
    - `datasets`

    - `LHS` - left hand side of symmetric linear system
    
    - `RHS` - right hand side of symmetric linear system

    - `c` - see `config`
    
    _Output Arguments_
    
    - `theta` - material parameters
    
    ---
    
    """
    for i in range(c.numIncrements):
        if i > 0:
            c.penaltyLp = c.penaltyLp * c.factorIncrements
        theta , at_least_one_converged = applyPenaltyLpRandomStart(datasets,LHS,RHS,c)
        checkLocalMinimumLp(datasets,theta,c)
        theta = applyThreshold(LHS,RHS,theta,c)
        fulfillEnergyRequirements = checkEnergyRequirementsRigorous(datasets,theta)
        if fulfillEnergyRequirements and at_least_one_converged:
            break
    return theta

#=====================================================================
# LOCAL MINIMUM:
#=====================================================================
def checkLocalMinimumLp(datasets,theta,c):
    """
    Check if the provided solution is a local minimum of the cost function.
    
    _Input Arguments_
    
    - `datasets`

    - `theta` - material parameters

    - `c` - see `config`
    
    _Output Arguments_
    
    - _none_
    
    ---
        
    """
    numChecks = 10
    perturbation_magnitude = 1e-3
    print('\n-----------------------------------------------------')
    print('Check if the solution is a local minimum.')
    print('Number of checks: ', numChecks)
    local_minimum = True
    _ , _ , Cost_theta = computeCostLp(datasets,theta,c)
    for i in range(numChecks):
        theta_perturbed = np.copy(theta)
        for idx in range(len(theta_perturbed)):
            if np.abs(theta[idx]) > c.threshold:
                theta_perturbed[idx] += perturbation_magnitude*(2*np.random.rand()-1)
#        print(theta_perturbed)
        _ , _ , Cost_theta_perturbed = computeCostLp(datasets,theta_perturbed,c)
        if Cost_theta > Cost_theta_perturbed:
            local_minimum = False
            print("Solution is not a local minimum.")
            break
    if local_minimum:
        print("Solution is likely to be a local minimum.")
    print('-----------------------------------------------------\n')

#=====================================================================
# COST:
#=====================================================================
def computeCostLp(datasets,theta,c):
    """
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
    
    """
    Cost_weak = 0.0
    Cost_penaltyLp = 0.0
    Cost_total = 0.0
    for d in range(len(datasets)):
        data = datasets[d]
        #Internal forces.
        non_dirichlet_nodes = ~(zipper(data.dirichlet_nodes))
        f_int_nodes = computeInternalForceTheta(data,theta,c)
        Cost_weak += np.sum(np.power(f_int_nodes[non_dirichlet_nodes],2))
    #    LHS_weak = computeWeakLHS(data,dim,numNodesPerElement)
    #    LHS_weak = LHS_weak[non_dirichlet_nodes,:]
    #    Cost_weak_test = np.transpose(theta).dot(np.transpose(LHS_weak)).dot(LHS_weak).dot(theta)
        #External forces.
        for reaction in data.reactions:
            force = reaction.force
            dofs = zipper(reaction.dofs)
            Cost_weak += c.balance*np.power(np.sum(f_int_nodes[dofs])-force,2)
        #Penalty terms.
        Cost_penaltyLp += c.penaltyLp*np.sum(np.power(np.abs(theta),c.p))
    Cost_total = Cost_weak + Cost_penaltyLp
    return Cost_weak , Cost_penaltyLp , Cost_total

#=====================================================================
# RESULTS:
#=====================================================================
def saveResultsLp(theta,c,counter_load=None):
    """
    Save the results and the chosen parameters in a text file.
    
    _Input Arguments_
    
    - `theta` - material parameters

    - `c` - see `config`

    - `counter_load` - current load step
    
    _Output Arguments_
    
    - _none_
    
    ---
    
    """
    if not os.path.exists("../EUCLID_results"):
        os.makedirs('../EUCLID_results')
    print('\n-----------------------------------------------------')
    print('Save results.')
    resultfile = "../EUCLID_results/results_" + c.saveResultsName + ".txt"
    f=open(resultfile,"a+") # append
    f.write("-------------------- Data --------------------\n")
    f.write(c.femDataPath + "\n")
    if counter_load == None:
        f.write("Load Steps: " + str(c.loadsteps) + "\n")
    else:
        f.write("Load Steps: " + str(c.loadsteps[0:counter_load+1]) + "\n")
    f.write("Additional Noise Level: " + str(c.noiseLevel) + "\n")
    f.write("-------------------- Hyperparameters --------------------\n")
    f.write("Balance Factor: " + str(c.balance) + "\n")
    f.write("Lp-Norm Penalty Factor (initial): " + str(c.penaltyLp_init) + "\n")
    f.write("Lp-Norm Penalty Factor: " + str(c.penaltyLp) + "\n")
    f.write("p: " + str(c.p) + "\n")
    f.write("Max. Number Penalty Increments: " + str(c.numIncrements) + "\n")
    f.write("Factor Penalty Increments: " + str(c.factorIncrements) + "\n")
    f.write("Number Initial Guesses: " + str(c.numGuesses) + "\n")
    f.write("Max. Number Iterations: " + str(c.numIterations) + "\n")
    f.write("Lowest Cost: " + str(c.lowestCost) + "\n")
    f.write("Lowest Cost Guess Identifier: " + str(c.lowestCostGuessID) + "\n")
    f.write("Threshold (Fixed-Point Iteration): " + str(c.threshold_iter) + "\n")
    f.write("Threshold: " + str(c.threshold) + "\n")
    f.write("-------------------- Results --------------------\n")
    f.write("Theta (Lp-Norm Penalty + Threshold): \n" + str(theta) + "\n")
    f.write("\n\n\n")
    f.close()
    print('-----------------------------------------------------\n')


    
    
    
