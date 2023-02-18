from joblib import Parallel, delayed
import pickle

from run_am_simulation import run_am_simulation

#############################################################
# User Input Parameters: 
#############################################################

# Simulation Parameters:
population_size = 1000
num_variants    = 100
num_generations = 20
num_iterations  = 20

seed = 1234567
max_cores = 20 # Use multiprocessing.cpu_count() to check # available

avoid_inbreeding = True # Avoid mating between relatives
save_covariances = True # Save full covariance matrix
save_history = True     # Save output from each generation

# Genetic Variables:
vg1 = .65 # Genetic variance of Trait 1
vg2 = .35 # Genetic variance of Trait 2
rg  = .20 # Genetic correlation between Traits 1 and 2

min_maf = .1 # Minimum minor allele frequency 
max_maf = .5 # Maximum minor allele frequency

prop_vg1_latent = 0.2 # Proportion of heritability not captured by Trait 1 PGS
prop_vg2_latent = 0.3 # Proportion of heritability not captured by Trait 2 PGS

# Mate Correlations: 
am11 = 0.00 # Mate1-Trait1, Mate2-Trait1
am12 = 0.05 # Mate1-Trait1, Mate2-Trait2
am21 = 0.10 # Mate1-Trait2, Mate2-Trait1
am22 = 0.15 # Mate1-Trait2, Mate2-Trait2

# Vertival Transmission: 
f11 = 0.0 # Offspring Trait 1 regressed on Parental Trait 1
f12 = 0.1 # Offspring Trait 1 regressed on Parental Trait 2
f21 = 0.2 # Offspring Trait 2 regressed on Parental Trait 1
f22 = 0.3 # Offspring Trait 1 regressed on Parental Trait 2

# Environmental Components
re = .1 # Environmental correlation between Trait 1 and Trait 2

#############################################################
# End User Input Variables
#############################################################

# Run the following code to complete the simulation:

# All the variables created above:
input_variables = [population_size, num_variants, num_generations, seed, max_cores, avoid_inbreeding, save_covariances, save_history, vg1, vg2, rg, min_maf, max_maf, prop_vg1_latent, prop_vg2_latent, am11, am12, am21, am22, f11, f12, f21, f22, re]

# Run the simulations in parallel:
results = Parallel(n_jobs=max_cores)(delayed(run_am_simulation)(*input_variables) for _ in range(num_iterations))

# Collect the results:
results_summary     = results[:][0]
covariance_results  = results[:][1]
history_results     = results[:][2]
phenotype_results   = results[:][3]

# Save_Output
with open('Bivariate_AM_Simulation_Results.pkl', 'wb') as f:
    pickle.dump((results_summary, covariance_results, history_results, phenotype_results), f)