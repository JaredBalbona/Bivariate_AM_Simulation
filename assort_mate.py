import numpy as np
import pandas as pd

from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment

def assort_mate(phen_summary, pheno_mate_cor, population_size, avoid_inbreeding):
    # Divide sample by sex, and rowbind the results:
    females_phen = phen_summary[phen_summary['Sex'] == 0]
    males_phen   = phen_summary[phen_summary['Sex'] == 1]

    fm_concat = np.concatenate((females_phen[['Y1', 'Y2']].values, males_phen[['Y1', 'Y2']].values), axis=0)

    # Calculate the inverse of the mates' phenotypic covariance matrix:
    inv_pheno_mate_cor = np.linalg.inv(pheno_mate_cor)

    # Calculate the Mahalanobis distances between each pair:
    fm_dists = cdist(fm_concat, fm_concat, 'mahalanobis', VI=inv_pheno_mate_cor)
    fm_dists = fm_dists[:len(females_phen), len(females_phen):] # The Male-Female quadrant of the distance matrix
    fm_dists[np.isnan(fm_dists)] = np.inf # Change NaN's to np.inf so that they'll be ignored

    # Use the Hungarian Algorithm to find the optimal assignment:
    row_ind, col_ind = linear_sum_assignment(fm_dists)

    # Extract the matched rows from females_phen[['Y1', 'Y2']] and males_phen[['Y1', 'Y2']]
    females_phen = pd.DataFrame(females_phen.iloc[row_ind])
    males_phen = pd.DataFrame(males_phen.iloc[col_ind])

    # Make the current indices the new index:
    males_phen = males_phen.reset_index(drop=True)
    females_phen = females_phen.reset_index(drop=True)

    # Check: 
    pheno_mate_cor; pd.concat([females_phen[['Y1','Y2']], males_phen[['Y1','Y2']]], axis=1).cov() # Looks good!

    # Make sure none of the mates have common ancestors:
    if avoid_inbreeding:
        no_sib_inbreeding = np.not_equal(males_phen[['Father_ID']], females_phen[['Father_ID']])
        no_cousin_inbreeding = pd.DataFrame(np.logical_and(
            np.not_equal(males_phen["Fathers_Father_ID"], females_phen["Fathers_Father_ID"]),
            np.not_equal(males_phen["Fathers_Father_ID"], females_phen["Mothers_Father_ID"])) &
            np.logical_and(np.not_equal(males_phen["Mothers_Father_ID"], females_phen["Fathers_Father_ID"]),
            np.not_equal(males_phen["Mothers_Father_ID"], females_phen["Mothers_Father_ID"])) &
            np.logical_and(np.not_equal(males_phen["Fathers_Mother_ID"], females_phen["Fathers_Mother_ID"]),
            np.not_equal(males_phen["Fathers_Mother_ID"], females_phen["Mothers_Mother_ID"])) &
            np.logical_and(np.not_equal(males_phen["Mothers_Mother_ID"], females_phen["Fathers_Mother_ID"]),
            np.not_equal(males_phen["Mothers_Mother_ID"], females_phen["Mothers_Mother_ID"])),
            columns=["no_cousin_inbreeding"])
        avoid_inbreeding = np.logical_and(no_sib_inbreeding.values, no_cousin_inbreeding.values)
        males_phen   = males_phen[avoid_inbreeding]
        females_phen = females_phen[avoid_inbreeding]

    # Create exactly the number of offspring desired
    n_offspring = np.random.poisson(population_size / len(males_phen), size=len(males_phen))

    # If this amount != population_size, add or subtract as needed:
    off_by_amount = (np.sum(n_offspring) - population_size)
    nonzero_indices = np.nonzero(n_offspring)
    indices_to_update = np.random.choice(nonzero_indices[0], abs(off_by_amount), replace=False)
    n_offspring[indices_to_update] -= (2 * (off_by_amount > 0) * 1) - 1

    # Add in the mate information:
    females_mate_ID = males_phen[['ID']].rename(columns={'ID':'females_mate_IDs'})
    males_mate_ID   = females_phen[['ID']].rename(columns={'ID':'males_mate_IDs'})

    females_phen = pd.concat([females_phen, females_mate_ID], axis =1).assign(n_offspring=n_offspring)
    males_phen = pd.concat([males_phen, males_mate_ID], axis =1).assign(n_offspring=n_offspring)

    return females_phen, males_phen

