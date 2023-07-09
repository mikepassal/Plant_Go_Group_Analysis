# %%
import pickle 
import pandas as pd
import itertools
import tqdm as tq
import Name_resolver 
import numpy as np

# %%
genes_pd = pickle.load(open('/data/passala/Generated_Tables/Comparing_all_go_groups_across_species/OrthoDB_python_files/orthodb_v11_genes.p','rb'))
og2genes_pd = pd.read_csv('/data/passala/OrthoDB_data/embryophyta_groups.csv',names = ['Orthogroup','Gene','Species'])

# %%
og2genes_pd['Ortholevel'] = og2genes_pd['Orthogroup'].str.split('t').str[1]

# %%
species_with_nets = pd.read_csv('/data/passala/Generated_Tables/Reference_tables/Species_name_resolver.csv')
species_with_nets = species_with_nets[:18]
# # species_with_nets = species_with_nets.drop(index = [13,14])
# # species_with_nets.loc[0,'Taxa ID'] = 39947 
taxa_to_keep = species_with_nets['Taxa ID'].to_list()


# %%
og2genes_only_cococonet = og2genes_pd.loc[og2genes_pd['Species'].isin(taxa_to_keep)]


# %%
groups_present = og2genes_only_cococonet['Species'].unique()


# %%
species_to_run_on = species_with_nets['Taxa ID'].loc[species_with_nets['Taxa ID'].isin(groups_present)].to_list()

# %%
species_combinations = list(itertools.combinations(species_to_run_on,2))

# %%
ncbi_mapping = pd.read_csv('/data/passala/OrthoDB_data/NCBI_data/merged_ncbi_to_orthodb_fixed_non_genesymbol.csv')

# %%
og2genes_only_cococonet = og2genes_only_cococonet.merge(right = ncbi_mapping[['Orthodb Gene','Symbol']], right_on = 'Orthodb Gene',left_on='Gene')

# %%
og2genes_only_cococonet.to_csv('/data/passala/OrthoDB_data/NCBI_data/og_2_Genes_with_network_id.csv',index = False)

# %%
for combo in tq.tqdm(species_combinations,desc='outer',position = 0):
    species_1,species_2 = combo[0],combo[1]

    species_1_name = Name_resolver.species_name_resolver(species_1,desired_type='common')
    species_2_name = Name_resolver.species_name_resolver(species_2,desired_type='common')

    first_species_ortho_groups = og2genes_only_cococonet.loc[og2genes_only_cococonet['Species'] == species_1]
    second_species_ortho_groups = og2genes_only_cococonet.loc[og2genes_only_cococonet['Species'] == species_2]
    shared_orthogroups = np.intersect1d(first_species_ortho_groups['Orthogroup'].unique(),second_species_ortho_groups['Orthogroup'].unique())

    list_of_orthogene_pds = []
    for orthogroup in tq.tqdm(shared_orthogroups,desc ='inner_loop',position= 0,leave = False):
        species_1_genes = first_species_ortho_groups['Gene'].loc[first_species_ortho_groups['Orthogroup']== orthogroup].to_list()
        species_2_genes = second_species_ortho_groups['Gene'].loc[second_species_ortho_groups['Orthogroup'] == orthogroup].to_list()
        all_gene_combos = list(itertools.product(species_1_genes,species_2_genes))
        current_orthogroup_pd = pd.DataFrame(columns = [f'{species_1_name} OrthoGene',f'{species_2_name} OrthoGene'],data = all_gene_combos)
        current_orthogroup_pd['Orthogroup'] = orthogroup
        list_of_orthogene_pds.append(current_orthogroup_pd)

    final_species_lineup = pd.concat(list_of_orthogene_pds)
    ncbi_added_once = final_species_lineup.merge(right = ncbi_mapping[['Orthodb Gene','Symbol']], right_on = 'Orthodb Gene',left_on=f'{species_1_name} OrthoGene')
    ncbi_added_once_clean= ncbi_added_once.drop(columns = 'Orthodb Gene').rename(columns = {'Symbol':f'{species_1_name} Symbol'})
    ncbi_added_twice = ncbi_added_once_clean.merge(right = ncbi_mapping[['Orthodb Gene','Symbol']], right_on = 'Orthodb Gene',left_on=f'{species_2_name} OrthoGene')
    full_final = ncbi_added_twice.drop(columns = 'Orthodb Gene').rename(columns = {'Symbol':f'{species_2_name} Symbol'})
    final_species_file_name = f"/home/passala/passala/OrthoDB_data/V_11_pairwise_maps_fixed_problem_species/{species_1_name}_to_{species_2_name}_ortholog_NM.csv"
    full_final.to_csv(final_species_file_name, index=False)


