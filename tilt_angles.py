import os
import utils_mda.seq_manipulation as seq_manipulation
from importlib import reload
reload(seq_manipulation)
import MDAnalysis
from MDAnalysis import analysis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis.analysis.hbonds
from tqdm import tqdm
from pathlib import Path
import sys
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA


def hbond_per_res(u, peptide_name, num_peptides, membrane_type):
    ##define lipid acceptors/donors. This has already been done for POPE and POPG 
    lipid_acceptors  = ['031','032','021','O22','O11','O13','O14','OC3','OC2','N']
    lipid_donors =['HO2','HO3','HN1','HN2','HN3']
    # peptide = seq_manipulation.get_peptide_info(u, peptide_name, 4)
    # total = peptide['Aminos']*peptide['pepnum']
    hbonds = []
    _, seq_dict = seq_manipulation.get_aa_sequence(u, [l for l in peptide_name.split("_")])
    print(_)
    print(seq_dict)
    for ts in tqdm(u.trajectory[-1]):
        for res_id, res_name in seq_dict.items():
            timecount1 = u.trajectory.time
            output1 = np.mean(u.select_atoms(f"name CA and resid {res_id}").positions[:,[2]].astype(float))
            print(output1)

    df = pd.DataFrame(hbonds)    
    # df.rename(columns={ df.columns[0]: "Time (ns)", df.columns[1]: "Residue",
    #                    df.columns[2]: "Hbond_num" }, inplace = True)
    # df['Time (ns)'] = df['Time (ns)']/1000
    # ##sum hydrogen bonds for each residue, unstack for wideform
    # df = df.groupby(['Residue','Time (ns)']).sum().unstack()
    # df = df.reset_index()
    # df['Section_Number'] = df['Residue'].str.replace('([A-Z]+)', '').astype(float)
    # df = df.sort_values('Section_Number')
    # df = df.drop('Section_Number',1)
    # df = df.set_index('Residue')
    # df.columns = df.columns.droplevel()
    return df



def calc_and_write_to_file(path, membrane_type, sims_type, single=None):
    results_directory = f"hbonds/may"
    Path(results_directory).mkdir(parents=True, exist_ok=True)
    #create selection
    all_hbond={}
    if single == "yes":
        if os.path.isdir(path):
            peptide_name = os.path.basename(path)
            print(f"Starting calculations for single peptide -- {peptide_name}")
            peptide_path = path
            u = seq_manipulation.get_universe(peptide_path)
            pep_h_bond = hbond_per_res(u, peptide_name, 4, membrane_type)
            pep_h_bond.to_csv(f"{results_directory}/hbonds_{peptide_name}_{membrane_type}")
            print(f"{peptide_name} --- DONE")
            all_hbond["peptide_name"] = pep_h_bond
    else:
        print(f"Starting calculations for multiple peptideß")
        for directory in tqdm(os.listdir(path)):
            folder = os.path.join(path,directory)
            if os.path.isdir(folder) and "_pg" in os.path.basename(folder):
                peptide_path = folder
                peptide_name = os.path.basename(peptide_path)
                my_file = Path(f"{results_directory}/hbonds_{peptide_name}_{membrane_type}")
                print("MY FILE")
                print(my_file)
                if my_file.is_file():
                    print(f"Results for {peptide_name} already exists! Skipping peptide.")
                else:
                    print(f"Starting calculations for peptide -- {peptide_name}")
                    u = seq_manipulation.get_universe(peptide_path)
                    pep_h_bond = hbond_per_res(u, peptide_name, 4, membrane_type)
                    pep_h_bond.to_csv(f"{results_directory}/hbonds_{peptide_name}_{membrane_type}.csv")

                    print(f"{peptide_name} --- DONE")
                    all_hbond["peptide_name"] = pep_h_bond
        if all_hbond:
            df_pdf = pd.DataFrame(all_hbond)
            df_pdf.to_csv(f"{results_directory}/hydr_all_new_{sims_type}_{membrane_type}.csv")
        else:
            print("No simulation folders found in the path provided!")
        
# print(seq_manipulation.get_peptide_info(u, "pep", 4))
# if __name__== "__main__":
#     if sys.argv[1]:
#         single="yes"
#     else:
#         single="no"
#     p = Path(f"/Volumes/miru_back/my_MDs/finished_short/{str(sys.argv[1])}")
#     calc_and_write_to_file(p, "pg", "my_sims", single="yes")


for pep in [ "WF1a_WF2"]:
    for memb in [ "pg"]:
        folder_name  = "popg" if memb == "pg" else "pepg"
        p=f"/Volumes/miru_backup/jade_2synergy/{folder_name}/{pep}_{memb}"
        calc_and_write_to_file(p, memb, "synergy", single="yes")