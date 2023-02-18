import numpy as np
import pandas as pd
from assort_mate import assort_mate
from reproduce import reproduce

def am_simulate(variant_info, num_generations, population_size, avoid_inbreeding, save_history, save_covariances, e_covmat, f_mat, a_mat, delta_mat, am_mat_list, y_covmat, vg_mat, num_variants):

    # Create observed and latent genotypes for generation 0:
    num_genotypes = num_variants * population_size
    obs_alleles = np.random.binomial(2, variant_info['MAF'], (int(population_size), num_variants))
    lat_alleles = np.random.binomial(2, variant_info['MAF'], (int(population_size), num_variants))

    # Check:  
    np.corrcoef(np.array(variant_info['MAF']).reshape(1, num_variants), obs_alleles.mean(axis=0)) # Looks good!

    # Create the PGS's for Gen 0:
    obs_pgs = obs_alleles @ variant_info[['Alpha1', 'Alpha2']]
    lat_pgs = lat_alleles @ variant_info[['Alpha1', 'Alpha2']]

    # Create observed and latent PGS matrices:
    obs_pgs_var = obs_pgs.cov()
    lat_pgs_var = lat_pgs.cov()

    # Create the components of Y:
    # Parental Effects:
    F = (np.random.multivariate_normal([0,0], y_covmat, population_size) @ f_mat.T) + (np.random.multivariate_normal([0,0], y_covmat, population_size) @ f_mat.T)
    F_y = F @ np.identity(2)

    # Environmental Effects: 
    E = np.random.multivariate_normal([0,0], e_covmat, population_size)
    E_y = E @ np.identity(2)

    # Genetic Effects:
    obs_pgs_y = obs_pgs @ delta_mat.T
    lat_pgs_y = lat_pgs @ a_mat.T

    # Create the phenotype from its component parts:
    Y = obs_pgs_y + lat_pgs_y + F_y + E_y

    # Check: 
    np.cov(E_y.T); e_covmat # Looks good!
    np.cov(F_y.T); 2 * (f_mat @ y_covmat @ f_mat.T) # Looks good!
    np.cov(obs_pgs_y.T); delta_mat @ vg_mat @ delta_mat.T # Looks good!
    np.cov(lat_pgs_y.T); a_mat @ vg_mat @ a_mat.T # Looks good!
    np.cov(Y.T); y_covmat # Looks good!

    # Create relative information for Gen0:
    ids = np.random.randint(10000000, 99999999, size=(population_size, 7)).reshape(population_size, 7)

    # Assign sexes to the people in Gen0, such that there are an equal number of males and females:
    sex = np.random.randint(2, size=population_size).reshape(-1, 1)

    # Create empty results matrix to be filled in as we go: 
    empty_mat = np.zeros((population_size, 28), dtype=float)

    phen_summary = np.column_stack((ids, sex, obs_pgs, lat_pgs, F, E, obs_pgs_y, lat_pgs_y, F_y, E_y, Y, empty_mat))

    phen_col_names = ['ID','Father_ID','Mother_ID','Fathers_Father_ID','Fathers_Mother_ID','Mothers_Father_ID','Mothers_Mother_ID','Sex','Obs_PGS_1','Obs_PGS_2','Lat_PGS_1','Lat_PGS_2','F1','F2','E1','E2','Obs_PGS_Y1','Obs_PGS_Y2','Lat_PGS_Y1','Lat_PGS_Y2','Fy1','Fy2','Ey1','Ey2','Y1','Y2','PGS_NT_Obs1','PGS_NT_Obs2','TP_Obs1','TP_Obs2','TM_Obs1','TM_Obs2','NTP_Obs1','NTP_Obs2','NTM_Obs2','NTM_Obs2','PGS_NT_Lat1','PGS_NT_Lat2','TP_Lat1','TP_Lat2','TM_Lat1','TM_Lat1','NTP_Lat1','NTP_Lat2','NTM_Lat1','NTM_Lat2','YP1','YP2','FP1','FP2','YM1','YM2','FM1','FM2']

    phen_summary = pd.DataFrame(phen_summary, columns = phen_col_names)

    # Calculate heritability information now that phen_summary has been created: 
    VA_obs = phen_summary[['Obs_PGS_1', 'Obs_PGS_2']].cov()
    VA_lat = phen_summary[['Lat_PGS_1', 'Lat_PGS_2']].cov()
    VY = phen_summary[['Y1', 'Y2']].cov()

    h2_obs = np.array(VA_obs)/ np.array(VY)
    h2_lat = np.array(VA_lat)/ np.array(VY)
    h2 = h2_obs + h2_lat

    # Save Gen0 info if desired:
    history = []
    if save_history:
        history.append({"Mates": [phen_summary], "Phen": [phen_summary], "obs_alleles": [obs_alleles], "lat_alleles": [lat_alleles]})

    # Save covariances if desired:
    covs = []
    if save_covariances:
        covs.append(None)

    # Save results for Gen0:
    results_summary = []
    results_summary.append({
        "Generation": 0, 
        "num_variants": num_variants,
        "Mate_Corr": am_mat_list[0], 
        "Population_Size": population_size, 

        "VA_Obs": VA_obs,
        "VA_Lat": VA_lat,

        "VY": VY,
        "VF": phen_summary[['F1', 'F2']].cov(),
        "VE": phen_summary[['E1', 'E2']].cov(),

        "h2": h2,
        "h2_Obs": h2_obs,
        "h2_Lat": h2_lat,

        "Y_Cov": float('nan'),
        "F_Cov": float('nan'),
        "E_Cov": float('nan'),
        
        "g": float('nan'),
        "h": float('nan'),
        "i": float('nan'),
        "w": float('nan'),
        "v": float('nan'),
        
        "Omega": float('nan'),
        "Gamma": float('nan'),
        "Theta_T": float('nan'),
        "Theta_NT": float('nan'),
    })

    # Now iteratively create the remaining generations:
    for current_gen in range(1, num_generations+1):
        
        # Gather mating information:
        mate_cor = am_mat_list[current_gen] 
        pheno_mate_cor = np.identity(4)
        pheno_mate_cor[0,1] = pheno_mate_cor[1,0] = pheno_mate_cor[3,2] = pheno_mate_cor[2,3] = phen_summary['Y1'].corr(phen_summary['Y2'])

        pheno_mate_cor[0:2, 2:4] = mate_cor
        pheno_mate_cor[2:4, 0:2] = mate_cor.T

        # Assortatively mate the individuals in this generation:
        female_mates, male_mates = assort_mate(phen_summary = phen_summary, pheno_mate_cor = pheno_mate_cor, population_size = population_size, avoid_inbreeding = avoid_inbreeding)

        # Have mates reproduce, creating new genotypes/ phenotypes for their offspring:
        phen_summary, obs_alleles, lat_alleles = reproduce(female_mates=female_mates, male_mates=male_mates, obs_pgs=obs_pgs, lat_pgs=lat_pgs, phen_summary=phen_summary, variant_info=variant_info, current_gen = current_gen, e_covmat=e_covmat, f_mat=f_mat, a_mat=a_mat, delta_mat=delta_mat, am_mat_list=am_mat_list, y_covmat=y_covmat, vg_mat=vg_mat, population_size=population_size, obs_alleles=obs_alleles, lat_alleles=lat_alleles)


        ###################

        # Haplotypic PGS Covariances:
        # Note that this script is not accounting for parental sex effects-- Hence, g and h are symmetric matrices. 
        pgs_covmat = np.cov(phen_summary.loc[:, ['Obs_PGS_TM_1', 'Obs_PGS_TP_1', 'Obs_PGS_NTM_1', 'Obs_PGS_NTP_1', 'Obs_PGS_TM_2', 'Obs_PGS_TP_2', 'Obs_PGS_NTM_2', 'Obs_PGS_NTP_2', 'Lat_PGS_TM_1', 'Lat_PGS_TP_1', 'Lat_PGS_NTM_1', 'Lat_PGS_NTP_1', 'Lat_PGS_TM_2', 'Lat_PGS_TP_2', 'Lat_PGS_NTM_2', 'Lat_PGS_NTP_2']].T)

        # Observed PGS covariance matrix:
        g_covmat = np.matrix([[
            pgs_covmat[ :4,  :4].mean(),
            pgs_covmat[ :4, 4:8].mean(),
            pgs_covmat[ :4, 4:8].mean(),
            pgs_covmat[4:8, 4:8].mean()]], 
            dtype = float).reshape(2,2)

        # Latent PGS covariance matrix:
        h_covmat = np.matrix([[
            pgs_covmat[ 8:12,  8:12].mean(),
            pgs_covmat[ 8:12, 12:16].mean(),
            pgs_covmat[ 8:12, 12:16].mean(),
            pgs_covmat[12:16, 12:16].mean()]], 
            dtype = float).reshape(2,2)

        # Observed PGS - Latent PGS covariance matrix:
        i_covmat = np.matrix([[
            pgs_covmat[ :4,  8:12].mean(),
            pgs_covmat[ :4, 12:16].mean(), # obs1 - latent 2
            pgs_covmat[4:8,  8:12].mean(), # obs2 - latent 1
            pgs_covmat[4:8, 12:16].mean()]], 
            dtype = float).reshape(2,2)

        # Calculate covariance matrices for each set of variables
        y_pgs_covmat = np.cov(phen_summary.loc[:, ['Y1', 'Y2', 'Y1_M','Y1_P','Y2_M', 'Y2_P', 'Obs_PGS_TM_1', 'Obs_PGS_TP_1', 'Obs_PGS_NTM_1', 'Obs_PGS_NTP_1', 'Obs_PGS_TM_2', 'Obs_PGS_TP_2', 'Obs_PGS_NTM_2', 'Obs_PGS_NTP_2', 'Lat_PGS_TM_1', 'Lat_PGS_TP_1', 'Lat_PGS_NTM_1', 'Lat_PGS_NTP_1', 'Lat_PGS_TM_2', 'Lat_PGS_TP_2', 'Lat_PGS_NTM_2', 'Lat_PGS_NTP_2']].T)

        # Covariance between parental phenotypes and observed PGSs:
        omega_covmat = np.matrix([[
            y_pgs_covmat[2:4, 6:10].mean(),
            y_pgs_covmat[2:4,10:14].mean(), # Y1 -> PGS2
            y_pgs_covmat[4:6, 6:10].mean(), # Y2 -> PGS1
            y_pgs_covmat[4:6, 6:10].mean()]], 
            dtype = float).reshape(2,2)

        # Covariance between parental phenotypes and latent PGSs:
        gamma_covmat = np.matrix([[
            y_pgs_covmat[2:4,14:18].mean(),
            y_pgs_covmat[2:4,18:22].mean(), # Y1 -> PGS2
            y_pgs_covmat[4:6,14:18].mean(), # Y2 -> PGS1
            y_pgs_covmat[4:6,18:22].mean()]], 
            dtype = float).reshape(2,2)

        # Covariance between offspring phenotype and transmitted PGSs:
        ThetaT_covmat = np.matrix([[
            y_pgs_covmat[0, 6:8].mean(),
            y_pgs_covmat[0, 10:12].mean(), 
            y_pgs_covmat[1, 6:8].mean(), 
            y_pgs_covmat[1, 10:12].mean() ]], 
            dtype = float).reshape(2,2)

        # Covariance between offspring phenotype and non-transmitted PGSs:
        ThetaNT_covmat = np.matrix([[
            y_pgs_covmat[0, 8:10].mean(),
            y_pgs_covmat[0, 12:14].mean(), 
            y_pgs_covmat[1, 8:10].mean(), 
            y_pgs_covmat[1, 12:14].mean() ]], 
            dtype = float).reshape(2,2)

        # All the other covariances:
        F_pgs_covmat = np.cov(phen_summary.loc[:, ['F1','F2','Obs_PGS_1', 'Obs_PGS_2','Lat_PGS_1','Lat_PGS_2']].T)

        w_covmat = F_pgs_covmat[:2, 2:4] # Cov between familial environment and observed PGS
        v_covmat = F_pgs_covmat[:2, 4:]  # Cov between familial environment and latent PGS
        F_covmat = F_pgs_covmat[:2, :2] # Familial environment covariance matrix
        E_covmat =  np.cov(phen_summary.loc[:, ['E1','E2']].T) # Unique environmental covariance matrix
        Y_covmat =  np.cov(phen_summary.loc[:, ['Y1_M','Y2_M','Y1_P','Y2_P','Y1','Y2']].T) # Phenotype covariance matrix

        VA_obs = phen_summary[['Obs_PGS_1', 'Obs_PGS_2']].cov() # Genetic variance due to observed PGSs
        VA_lat = phen_summary[['Lat_PGS_1', 'Lat_PGS_2']].cov() # Genetic variance due to latent PGSs
        VY = phen_summary[['Y1', 'Y2']].cov() # Phenotypic Variance
        h2_obs = np.array(VA_obs)/ np.array(VY) # Heritability due to observed PGSs
        h2_lat = np.array(VA_lat)/ np.array(VY) # Heritability due to latent PGSs
        h2 = h2_obs + h2_lat # Overall heritability

        # Collect results:
        results_summary.append({
            "Generation": current_gen, 
            "num_variants": num_variants,
            "Mate_Corr": mate_cor, 
            "Population_Size": population_size, 

            "VA_Obs": VA_obs,
            "VA_Lat": VA_lat,

            "VY": VY,
            "VF": F_covmat,
            "VE": E_covmat,

            "h2": h2,
            "h2_Obs": h2_obs,
            "h2_Lat": h2_lat,

            "Y_Cov": Y_covmat,
            "F_Cov": F_covmat,
            "E_Cov": E_covmat,
            
            "g": g_covmat,
            "h": h_covmat,
            "i": i_covmat,
            "w": w_covmat,
            "v": v_covmat,
            
            "Omega": omega_covmat,
            "Gamma": gamma_covmat,
            "Theta_T": ThetaT_covmat,
            "Theta_NT": ThetaNT_covmat,
        })

        # Save information for this generation, if desired:
        if save_history:
            history.append({"Mates": [phen_summary], "Phen": [phen_summary], "obs_alleles": [obs_alleles], "lat_alleles": [lat_alleles]})

        # Save covariances for this generation, if desired:
        if save_covariances:
            covs.append(phen_summary[['Obs_PGS_1', 'Obs_PGS_2', 'Lat_PGS_1','Lat_PGS_2', 'Obs_NT_PGS_1', 'Obs_NT_PGS_2', 'Lat_NT_PGS_1', 'Lat_NT_PGS_2', 'Obs_PGS_TM_1', 'Obs_PGS_TM_2', 'Obs_PGS_TP_1', 'Obs_PGS_TP_2', 'Obs_PGS_NTM_1', 'Obs_PGS_NTM_2', 'Obs_PGS_NTP_1', 'Obs_PGS_NTP_2', 'Lat_PGS_TM_1', 'Lat_PGS_TM_2', 'Lat_PGS_TP_1', 'Lat_PGS_TP_2', 'Lat_PGS_NTM_1', 'Lat_PGS_NTM_2', 'Lat_PGS_NTP_1', 'Lat_PGS_NTP_2', 'F1','F2','E1','E2', 'Obs_PGS_Y1','Obs_PGS_Y2','Lat_PGS_Y1','Lat_PGS_Y2','Fy1','Fy2','Ey1','Ey2','Y1','Y2','Y1_M', 'Y2_M', 'Y1_P', 'Y2_P', 'F1_M', 'F2_M', 'F1_P', 'F2_P']].cov())    

    return results_summary, covs, history, phen_summary 