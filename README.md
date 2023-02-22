## Multivariate Simulation of Human Population Growth

----

### Broad Summary:

These scripts collectively simulate the genetic and phenotypic development of a human population over time. The most notable feature of this simulation is that it is able to *bivariately* model parenting effects and mating patterns— two things which can significantly alter the genotypic and phenotypic architectures of traits in a given population, and which are often ignored by genetic studies/ simulations. 

Imagine, for example, that male-female partners tend to be highly correlated on height, but only slightly correlated on IQ. Now imagine that taller males tend to couple up higher IQ females, but that female height has no relation to their male partner's IQ level. This phenomenon is known *assortative mating*, and there is actually widespread evidence for it occurring [within](https://www.biorxiv.org/content/10.1101/2022.03.19.484997v2.full) and [between](https://www.science.org/doi/abs/10.1126/science.abo2059) traits (including with male height and female IQ!). Similarly, you can imagine a situation where having a parent high in one trait (e.g., education level) might lead to the development of the same trait in the offspring, or of a different trait (e.g., anxiety). While this parent-offspring similarity could be due to shared genetic factors, it could also be due to shared environmental factors as well (a phenomenon known as *vertical transmission*). 

This simulation allows for the user to model assortative mating and vertical transmission both within- and between- traits. In doing so, they produce results that are much closer to real life, and allow for the user to test other models for bias due to these confounding effects. 

The mathematical expectations in these scripts are largely taken from [previously published work](https://link.springer.com/article/10.1007/s10519-020-10032-w) by our research group in which we used structural equation modeling to estimate familial effects. However, this simulation offers a vectorized implementation of the math from our prevoius publication, allowing us to model many many more relationships without becoming computationally intractable. For further information on the math involved— including some principles multivariate path tracing, descriptions of the matrices used, and explanations for why matrices are diagonal vs. symmetric vs. full, please see [our guide on Overleaf](https://www.overleaf.com/read/vjvshhnmfcdq).

Please feel free to reach out to me at jaba5258@colorado.edu with any quesitons— Thanks!

----

### Specific Steps:

 1. #### main.py:  ####
   - This script is pretty straightforward. Here, the user is able to supply the parameters for their simulation, including details of the simulation itself (e.g., the seed and the number of iterations), the trait's components (e.g., how heritable it is), and the desired correlation between mates. Once the parameters have all been specified, the simulation is run in parallel across a user-specified number of cores, the results are properly formatted, and the output is saved.

 2. #### run_am_simulation.py:  ####
   - This script takes the user-supplied parameters and sets the stage fo the simulation. This includes building all of the necessary matrices and adding some noise to the phenotypic (co)variances due to pleiotropy. Additionally, it is here that I'm creating the "grand truth" for our genetic variables by setting the allele frequencies (drawn from a uniform distribution) and effect sizes (drawn from a multivariate normal distribution) at Time 0. 

 3. #### am_simulate.py:  ####
   - 'am_simulate' begins by creating genotypes and phenotypes from scratch for our population's 0th generation. In doing so, I am also creating empty results dataframes that will be populated as the simultaion goes on. Once Gen0 has been created, the simulation then iteratively calls 'assort_mate' and 'reproduce' to create the remaining generations before outputting the results.  
 
  4. #### assort_mate.py:  ####
   - Here, I'm pairing up males and females in our sample based on their respective values for Trait 1 and Trait2. To do this, I will calculate the pairwise [Mahalanobis distance](https://en.wikipedia.org/wiki/Mahalanobis_distance) between all opposite sex pairs in our sample, taking the desired Trait 1 - Trait 2 covariance matrix into account. Once the distance between all pairs has been calculated, I'm using the [Hungarian Algorithm](https://en.wikipedia.org/wiki/Hungarian_algorithm) to find the optimal assignment for each individual. Once paired, I am randomly choosing a number of offspring for each couple such that the next generation will be the exact size desired.
     - Note: While this simulation (and the study it is based on) are explicitly modeling opposite sex pairs, this is in no way to suggest that same sex relationships are not important/ worth examining (they are!). It is simply because we are interested here in the genetic consequences that result when an assorted mate pair procreates. 

  5. #### reproduce.py:  ####
   - Once all of the individuals have been paired up and their number of offspring has been assigned, it is time to create the genotypes and phenotypes of each of the offspring. This process is very similar to what was done in am_simulate to create generation 0—  the primary difference being that both are now a function of their parents' genotypes/ phenotypes. Once this has been completed, the information for each individual is added the dataframe created in am_simulate.
