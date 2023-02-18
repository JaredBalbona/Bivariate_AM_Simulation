import numpy as np
import pandas as pd
from am_simulate import am_simulate

def run_am_simulation(population_size, num_variants, num_generations, seed, max_cores,  avoid_inbreeding, save_covariances, save_history, vg1, vg2, rg, min_maf, max_maf, prop_vg1_latent, prop_vg2_latent, am11, am12, am21, am22, f11, f12, f21, f22, re):  
    np.random.seed(seed)
        
    # Create VG Matrix:
    vg_mat = np.matrix([[1, rg, rg, 1]], dtype = float).reshape(2,2)

    # Observed Genetic Variables #
    # Proportion of genetic variance that is observed (i.e., captured in the PGS):
    vg_obs_1 = vg1 * (1 - prop_vg1_latent)
    vg_obs_2 = vg2 * (1 - prop_vg2_latent)

    delta11 = np.sqrt(vg_obs_1)
    delta22 = np.sqrt(vg_obs_2)
    delta_mat = np.matrix([[delta11, 0, 0, delta22]], dtype = float).reshape(2,2)

    gen_cov_obs = rg * delta11 * delta22
    vg_obs_mat = delta_mat @ vg_mat @ delta_mat.T

    # Latent Genetic Variables # 
    # Proportion of genetic variance that is latent (i.e., not captured in the PGS):
    vg_lat_1 = vg1 * prop_vg1_latent
    vg_lat_2 = vg2 * prop_vg2_latent

    a11 = np.sqrt(vg_lat_1)
    a22 = np.sqrt(vg_lat_2)
    a_mat = np.matrix([[a11, 0, 0, a22]], dtype = float).reshape(2,2)

    gen_cov_lat = rg * a11 * a22
    vg_lat_mat = a_mat @ vg_mat @ a_mat.T

    # Total genetic covariance:
    gen_cov = gen_cov_obs + gen_cov_lat

    # Add in background Y correlations due to pleiotropy, environmental correlation, and xVT (i.e., noise):
    ycov = ecov = 0
    for i in range(8):
        for j in range(10):
            ycov =  gen_cov + ecov + 2 * ((f11 * f21) + (f12 * f22) + (f11 * ycov * f22) + (f12 * ycov * f21))
        vf1 = 2 * (f11**2 + f12**2) + 2 * (f11 * ycov * f12)
        vf2 = 2 * (f22**2 + f21**2) + 2 * (f22 * ycov * f21)
        ve1 = 1 - (vg1 + vf1)
        ve2 = 1 - (vg2 + vf2)
        ecov = re * np.sqrt(ve1 * ve2)

    # Create Environmental Covariance Matrix:
    e_covmat = np.matrix([[ve1, ecov, ecov, ve2]], dtype = float).reshape(2,2)

    # Create Mate Correlation Matrix:
    am_mat = np.matrix([[am11, am12, am21, am22]], dtype = float).reshape(2,2)
    am_mat_list = [am_mat] * (num_generations + 1)

    # Create Phenotypic Covariance Matrix:
    y_covmat = np.matrix([[1, ycov, ycov, 1]], dtype = float).reshape(2,2)

    # Create Parental Effects Matrix:
    f_mat = np.matrix([[f11, f12, f21, f22]], dtype = float).reshape(2,2)
    f_covmat = 2 * (f_mat @ y_covmat @ f_mat.T)

    # Check to make sure everything is correct:
    y_covmat; a_mat @ vg_mat @ a_mat.T + delta_mat @ vg_mat @ delta_mat.T + e_covmat + f_covmat # Looks good! 
    vg_obs_mat; np.matrix([[vg_obs_1, gen_cov_obs, gen_cov_obs, vg_obs_2]], dtype = float).reshape(2,2) # Looks good! 
    vg_lat_mat; np.matrix([[vg_lat_1, gen_cov_lat, gen_cov_lat, vg_lat_2]], dtype = float).reshape(2,2) # Looks good! 

    # Create allele frequencies and effect sizes for the causal variants:
    maf_vector = np.random.uniform(min_maf, max_maf, num_variants)  
    gentp_var = maf_vector * (2 * (1 - maf_vector))

    var_effects = (np.array(np.sqrt(1 / (num_variants * gentp_var))))
    var_effects.shape = (num_variants, 1) # Make into a column vector

    alphas_pre = np.random.multivariate_normal([0,0], vg_mat, num_variants)
    alphas = alphas_pre * np.concatenate((var_effects, var_effects), axis=1)

    variant_info = pd.DataFrame({"MAF": maf_vector, "Alpha1": alphas[:, 0], "Alpha2": alphas[:, 1]})

    # Run Simulation: 
    results_summary, covariances, history, phenotype_summary = am_simulate(variant_info=variant_info, num_generations=num_generations, population_size=population_size, avoid_inbreeding=avoid_inbreeding, save_history=save_history, save_covariances=save_covariances, e_covmat=e_covmat, f_mat=f_mat, a_mat=a_mat, delta_mat=delta_mat, am_mat_list=am_mat_list, y_covmat=y_covmat, vg_mat=vg_mat, num_variants=num_variants)

    # Gather Ouput:
    # results = [am_data["results_summary"], am_data["covs"], am_data["phen_summary"]]
    # results = [f"Summary.it{iteration_num}" for iteration_num in results]
    # results = [f"Covs.it{iteration_num}" for iteration_num in results[1:]]
    # results = [f"PHEN.it{iteration_num}" for iteration_num in results[2:]]

    return results_summary, covariances, history, phenotype_summary