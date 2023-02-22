## Multivariate Simulation of Human Population Growth

----

### Broad Summary:

These scripts collectively simulate the genetic and phenotypic development of a human population over time. Importantly, they are able to model *assortative mating*—the phenomenon in which individuals that are high in a given trait are more likely to mate with other individuals who are high on that same trait— both within-trait and between traits. Similarly, they are also able to model within- and cross-trait *vertical transmission*— the phenomenon when parental traits directly influence an offspring trait through non-genetic means. By operating bivariately, this simulation is able to produce results that are much closer to reality than would be possible with a univariate simulation.

The mathematical expectations in these scripts are largely taken from [previously published work](https://link.springer.com/article/10.1007/s10519-020-10032-w) by our research group in which we used structural equation modeling to estimate familial effects. However, this simulation offers a vectorized implementation of the math from our prevoius publication, allowing us to model many many more relationships without becoming computationally intractable. 

For further information on the math involved— including some principles multivariate path tracing, descriptions of the matrices used, and explanations for why matrices are diagonal vs. symmetric vs. full, please see [our guide on Overleaf](https://www.overleaf.com/read/vjvshhnmfcdq).

Please feel free to reach out to me at jaba5258@colorado.edu with any quesitons— Thanks!

----

### Specific Steps:

 1. #### main.py:  ####
   - This script is pretty straightforward. Here, the user is able to supply the parameters for their simulation, including details of the simulation itself (e.g., the seed and the number of iterations), the trait's components (e.g., how heritable it is), and the desired correlation between mates. Once the parameters have all been specified, the simulation is run in parallel across a user-specified number of cores, the results are properly formatted, and the output is saved.

 2. #### run_am_simulation.py:  ####
   - This script takes the user-supplied parameters and sets the stage fo the simulation. This includes building all of the necessary matrices and adding some noise to the phenotypic (co)variances due to pleiotropy. Additionally, it is here that we're creating the "grand truth" for our genetic variables by setting the allele frequencies (drawn from a uniform distribution) and effect sizes (drawn from a multivariate normal distribution) at Time 0. 

 3. #### am_simulate.py:  ####
   - 'am_simulate' begins by creating genotypes and phenotypes from scratch for our population's 0th generation. In doing so, we are also creating empty results dataframes that will be populated as the simultaion goes on. Once Gen0 has been created, the simulation then iteratively calls 'assort_mate' and 'reproduce' to create the remaining generations before outputting the results.  
 
  4. #### assort_mate.py:  ####
   - 
