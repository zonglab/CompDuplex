import os, sys
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt

def skip_header(PATH_in):
    with open(PATH_in, "r") as FILE_in:
        for i, line in enumerate(FILE_in):
            if line[0] != "#": return i

#############################################################################################################################
### config input parameters
#############################################################################################################################
WORK___DIR=sys.argv[1]
BULK__PATH=sys.argv[2]

#############################################################################################################################
### wait until the previous steps have finished
#############################################################################################################################
LIST_bcd = [f"{WORK___DIR}/1_split/complete/{i}{j}{k}.txt" for i in "ATCG" for j in "ATCG" for k in "ATCG"]
LIST_chr = [f"{WORK___DIR}/1_split/complete/chr{i}.txt" for i in list(range(1,23))+["X","Y"]]
while True:
    finish_1 = all([os.path.exists(i) for i in LIST_bcd])
    finish_2 = all([os.path.exists(i) for i in LIST_chr])
    if  finish_1 | finish_2: break
    os.system("sleep 30s")
os.system(f"cat {WORK___DIR}/1_split/*.vis > {WORK___DIR}/14_a4s2_gather.vis")

#############################################################################################################################
### read ref bulk germlines
#############################################################################################################################
PATH_germ_het = f"{BULK__PATH}/15_fin_hetero.vcf"
PATH_germ_hom = f"{BULK__PATH}/15_fin___homo.vcf"

df_germhom = pd.read_csv(PATH_germ_hom, sep="\t", skiprows=skip_header(PATH_germ_hom), header=None).loc[:,[0,1,3,4,2]]
df_germhet = pd.read_csv(PATH_germ_het, sep="\t", skiprows=skip_header(PATH_germ_het), header=None).loc[:,[0,1,3,4,2]]

df_germhom["bulk"] = "homo"
df_germhet["bulk"] = "hetero"
df_germlne = pd.merge(df_germhet, df_germhom, on=[0,1,3,4,"bulk"], how="outer").loc[:,[0,1,3,4,"bulk"]]
df_germlne.columns = ["chr","pos","ref","alt","bulk"]

#############################################################################################################################
### read detected mutations
#############################################################################################################################
PATH_det_muts = f"{WORK___DIR}/14_a4s2_gather.vis"
PATH_det_hets = f"{WORK___DIR}/15_1_germline_hets.bdy"
PATH_det_homs = f"{WORK___DIR}/15_2_germline_homs.bdy"
PATH_de_novos = f"{WORK___DIR}/15_3_not_germlines.bdy"
PATH__damages = f"{WORK___DIR}/16_1_damages_error.bdy"

df_det_totl = pd.read_csv(PATH_det_muts, skiprows=skip_header(PATH_det_muts), sep="\t", header=None).loc[:,[0,1,2,3,4,6,7,8,9]]
df_det_totl.columns = ["chr","pos","id","ref","alt","aMsN","dist","call","stats"]

#############################################################################################################################
### find ds mutations
#############################################################################################################################
df_det_muts = df_det_totl.loc[df_det_totl["call"] == "ds_Mut",["chr","pos","id","ref","alt","dist","stats"]]
df_det_muts.loc[:,["dist"]] = df_det_muts["dist"].str.replace("[\[\] ]", "", regex=True).str.split(",", expand=True).astype(int).min(axis=1)

#############################################################################################################################
### intersect with reference mutations
#############################################################################################################################
df_det_germ = pd.merge(df_det_muts, df_germlne, on=["chr","pos","ref","alt"], how="left")
df_det_germ.loc[df_det_germ["bulk"].isna(),"bulk"] = "de_novo"

#############################################################################################################################
### count DNA damages
#############################################################################################################################
df__damages = df_det_totl.loc[(df_det_totl["call"]=="Dmg_Tran")|(df_det_totl["call"]=="Dmg____0"),["chr","pos","id","ref","alt","aMsN","dist","stats"]]

df__damages.loc[:,["dist"]] = df__damages["dist"].str.replace("[\[\] ]", "", regex=True).str.split(",", expand=True).astype(int).min(axis=1)
dmg_strands = df__damages["aMsN"].str.replace("[\[\] ]", "", regex=True).str.split(",", expand=True).astype(int)
df__damages.loc[:,["stnd"]] = ((dmg_strands[0]==0)&(dmg_strands[2]>=2)).map({True:"rev", False:"for"})
df__damages = df__damages.loc[:,["chr","pos","ref","alt","stnd","dist"]]

#############################################################################################################################
### write output
#############################################################################################################################
df_det_germ = df_det_germ.sort_values(["chr","pos"])
df__damages = df__damages.sort_values(["chr","pos"])
print(sum(df_det_germ["bulk"]=="de_novo"), sum(df_det_germ["bulk"]== "hetero"), df__damages.shape[0])
df_det_germ.loc[df_det_germ["bulk"]== "hetero",].to_csv(PATH_det_hets, index=False, header=False, sep="\t")
df_det_germ.loc[df_det_germ["bulk"]==   "homo",].to_csv(PATH_det_homs, index=False, header=False, sep="\t")
df_det_germ.loc[df_det_germ["bulk"]=="de_novo",].to_csv(PATH_de_novos, index=False, header=False, sep="\t")
df__damages.to_csv(PATH__damages, index=False, header=False, sep="\t")