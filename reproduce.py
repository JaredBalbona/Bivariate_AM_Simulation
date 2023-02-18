import numpy as np
import pandas as pd

def reproduce(female_mates, male_mates, obs_pgs, lat_pgs, phen_summary, variant_info, current_gen, e_covmat, f_mat, a_mat, delta_mat, am_mat_list, y_covmat, vg_mat, population_size, obs_alleles, lat_alleles):

    # Line up the male and female genotypes with the mating array:
    females_phen = female_mates[female_mates["n_offspring"] > 0]
    males_phen   = male_mates[male_mates["n_offspring"] > 0]

    merge_f = pd.merge(females_phen, phen_summary.reset_index(), on="ID", how="inner")
    merge_m = pd.merge(males_phen, phen_summary.reset_index(), on="ID", how="inner")

    indices_f = merge_f.loc[:, "index"].tolist()
    indices_m = merge_m.loc[:, "index"].tolist()

    obs_alleles_f1 = pd.DataFrame(obs_alleles).loc[indices_f]
    obs_alleles_m1 = pd.DataFrame(obs_alleles).loc[indices_m]

    lat_alleles_f1 = pd.DataFrame(lat_alleles).loc[indices_f]
    lat_alleles_m1 = pd.DataFrame(lat_alleles).loc[indices_m]

    # Create 1 row per new offspring:
    current_gen_offspring = females_phen.iloc[:, females_phen.columns.get_loc("n_offspring")] # Same for M and F
    offspring_repeated = np.repeat(range(females_phen.shape[0]), current_gen_offspring) # Also same for M and F

    obs_alleles_f = obs_alleles_f1.iloc[offspring_repeated, :].reset_index(drop=True)
    obs_alleles_m = obs_alleles_m1.iloc[offspring_repeated, :].reset_index(drop=True)
    lat_alleles_f = lat_alleles_f1.iloc[offspring_repeated, :].reset_index(drop=True)
    lat_alleles_m = lat_alleles_m1.iloc[offspring_repeated, :].reset_index(drop=True)

    #Create haplotypes for each of the offspring:
    def create_haplotypes(parental_alles):
        # Wildcard to randomize which allele is transmitted to the offspring:
        T_allele_wildcard = np.random.randint(2, size=(parental_alles.shape[0], parental_alles.shape[1]))

        # Separate parents by zygosity:
        heterozygotes = (parental_alles == 1) * 1
        homozygotes   = (parental_alles == 2) * 1

        # For heterozygotes, randomly choose which allele is transmitted vs. non-transmitted:
        het_T_allele  = T_allele_wildcard  * heterozygotes
        het_NT_allele = ((T_allele_wildcard == 0) * 1) * heterozygotes 

        # Return transmitted and non-transmitted haplotypes:
        return het_T_allele + homozygotes, het_NT_allele + homozygotes

    # Observed PGS vs. Latent PGS x Maternal vs. Paternal:
    obs_hap_tm, obs_hap_ntm = create_haplotypes(obs_alleles_f)
    obs_hap_tp, obs_hap_ntp = create_haplotypes(obs_alleles_m)
    lat_hap_tm, lat_hap_ntm = create_haplotypes(lat_alleles_f)
    lat_hap_tp, lat_hap_ntp = create_haplotypes(lat_alleles_m)

    # Create genotypes for the new offspring: 
    obs_alleles = obs_hap_tm + obs_hap_tp
    lat_alleles = lat_hap_tm + lat_hap_tp

    # Create 1 row per new offspring from the phenotype matrix:
    # Note: To clarify, my naming convention up to this point has been 'm' for male and 'f' for female--Not 'mother' and 'father' (which can be confusing).
    mothers_phen = females_phen.iloc[offspring_repeated, :].reset_index(drop=True)
    fathers_phen = males_phen.iloc[offspring_repeated, :].reset_index(drop=True)

    # Create observed haplotypic PGS's:
    obs_pgs = obs_alleles @ np.matrix(variant_info[["Alpha1","Alpha2"]], dtype = float)  
    obs_pgs_tm = obs_hap_tm @ np.matrix(variant_info[["Alpha1","Alpha2"]], dtype = float)
    obs_pgs_tp = obs_hap_tp @ np.matrix(variant_info[["Alpha1","Alpha2"]], dtype = float)
    obs_pgs_ntm = obs_hap_ntm @ np.matrix(variant_info[["Alpha1","Alpha2"]], dtype = float)
    obs_pgs_ntp = obs_hap_ntp @ np.matrix(variant_info[["Alpha1","Alpha2"]], dtype = float)

    obs_pgs_nt = obs_pgs_ntm + obs_pgs_ntp

    # Create latent haplotypic PGS's:
    lat_pgs = lat_alleles @ np.matrix(variant_info[["Alpha1","Alpha2"]], dtype = float)  
    lat_pgs_tm = lat_hap_tm @ np.matrix(variant_info[["Alpha1","Alpha2"]], dtype = float)
    lat_pgs_tp = lat_hap_tp @ np.matrix(variant_info[["Alpha1","Alpha2"]], dtype = float)
    lat_pgs_ntm = lat_hap_ntm @ np.matrix(variant_info[["Alpha1","Alpha2"]], dtype = float)
    lat_pgs_ntp = lat_hap_ntp @ np.matrix(variant_info[["Alpha1","Alpha2"]], dtype = float)

    lat_pgs_nt = lat_pgs_ntm + lat_pgs_ntp

    all_pgs = np.column_stack([obs_pgs, lat_pgs, obs_pgs_nt, lat_pgs_nt, obs_pgs_tm, obs_pgs_tp, obs_pgs_ntm, obs_pgs_ntp, lat_pgs_tm, lat_pgs_tp, lat_pgs_ntm, lat_pgs_ntp])

    # Check:
    np.sum(obs_pgs - (obs_pgs_tm + obs_pgs_tp)) # Basically 0!
    np.sum(lat_pgs - (lat_pgs_tm + lat_pgs_tp)) # Also basically 0!

    #Create components of Y (the above multiplied by their path coefficients):
    F = np.matrix(mothers_phen[["Y1","Y2"]], dtype = float) @ f_mat.T + np.matrix(fathers_phen[["Y1","Y2"]], dtype = float) @ f_mat.T
    E = np.random.multivariate_normal([0,0], e_covmat, population_size)

    obs_pgs_y = obs_pgs @ delta_mat.T
    lat_pgs_y = lat_pgs @ a_mat.T

    F_y = F @ np.identity(2)
    E_y = E @ np.identity(2)

    Y = (obs_pgs_y + lat_pgs_y + F_y + E_y)

    y_components = np.column_stack([F, E, obs_pgs_y, lat_pgs_y, F_y, E_y, Y])

    # Create parental y components:
    mother_Y = females_phen[['Y1','Y2']].iloc[offspring_repeated, :].reset_index(drop=True)
    father_Y =   males_phen[['Y1','Y2']].iloc[offspring_repeated, :].reset_index(drop=True)

    mother_F = females_phen[['F1','F2']].iloc[offspring_repeated, :].reset_index(drop=True)
    father_F =   males_phen[['F1','F2']].iloc[offspring_repeated, :].reset_index(drop=True)

    parental_y_components = np.column_stack([mother_Y, father_Y, mother_F, father_F])

    # Create relative ID information for this generation:
    ids = np.random.randint(10000000, 99999999, size=(population_size))
    mother_id = females_phen[['ID']].iloc[offspring_repeated, :].reset_index(drop=True)
    father_id = males_phen[['ID']].iloc[offspring_repeated, :].reset_index(drop=True)

    mothers_mother_id = females_phen[['Mother_ID']].iloc[offspring_repeated, :].reset_index(drop=True)
    fathers_mother_id =   males_phen[['Mother_ID']].iloc[offspring_repeated, :].reset_index(drop=True)
    mothers_father_id = females_phen[['Father_ID']].iloc[offspring_repeated, :].reset_index(drop=True)
    fathers_father_id =   males_phen[['Father_ID']].iloc[offspring_repeated, :].reset_index(drop=True)

    all_ids = np.column_stack([ids, mother_id, father_id, mothers_mother_id, mothers_father_id, fathers_mother_id, fathers_father_id])

    # Randomly assign sexes for this generation such that there are an equal number of males and females:
    sex = np.random.randint(2, size=population_size).reshape(-1, 1)

    # Create the column names for the final dataframe using a dictionary
    column_names = {
        'all_ids': ['ID', 'Mother_ID', 'Father_ID', 'Mothers_Mother_ID', 'Mothers_Father_ID', 'Fathers_Mother_ID', 'Fathers_Father_ID'],
        'sex': ['Sex'],
        'all_pgs': ['Obs_PGS_1', 'Obs_PGS_2', 'Lat_PGS_1','Lat_PGS_2', 'Obs_NT_PGS_1', 'Obs_NT_PGS_2', 'Lat_NT_PGS_1', 'Lat_NT_PGS_2', 'Obs_PGS_TM_1', 'Obs_PGS_TM_2', 'Obs_PGS_TP_1', 'Obs_PGS_TP_2', 'Obs_PGS_NTM_1', 'Obs_PGS_NTM_2', 'Obs_PGS_NTP_1', 'Obs_PGS_NTP_2', 'Lat_PGS_TM_1', 'Lat_PGS_TM_2', 'Lat_PGS_TP_1', 'Lat_PGS_TP_2', 'Lat_PGS_NTM_1', 'Lat_PGS_NTM_2', 'Lat_PGS_NTP_1', 'Lat_PGS_NTP_2'],
        'y_components': ['F1','F2','E1','E2', 'Obs_PGS_Y1','Obs_PGS_Y2','Lat_PGS_Y1','Lat_PGS_Y2','Fy1','Fy2','Ey1','Ey2','Y1','Y2'],
        'parental_y_components': ['Y1_M', 'Y2_M', 'Y1_P', 'Y2_P', 'F1_M', 'F2_M', 'F1_P', 'F2_P']
    }

    # Concatenate the dataframes using np.concatenate and the column names dictionary
    phen_summary = pd.DataFrame(np.concatenate([
        all_ids, sex, all_pgs, y_components, parental_y_components
    ], axis=1), columns=sum(column_names.values(), []))

    return phen_summary, obs_alleles, lat_alleles

