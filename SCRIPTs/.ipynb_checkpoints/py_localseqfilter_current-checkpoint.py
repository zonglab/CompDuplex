import sys
import numpy  as np
import pandas as pd

def skip_header(in_path, header="#"):
    with open(in_path, "r") as FILE_in:
        for i, line in enumerate(FILE_in):
            if line[0] != header: return i

PATH_vcf = sys.argv[1]
PATH_seq = sys.argv[2]
PATH_out = sys.argv[3]

df___vcf = pd.read_csv(PATH_vcf, skiprows=skip_header(PATH_vcf), sep="\t", header=None) 
df___seq = pd.read_csv(PATH_seq, skiprows=skip_header(PATH_seq), sep="\t", header=None) 
df___vcf[4] = df___vcf[4].str.replace(",<*>","")

df___vcf = df___vcf.sort_values([0,1])
df___seq = df___seq.sort_values([0,1])

df___seq.loc[:,"head"] = df___seq.loc[:,2].str.slice( 0,10).str.upper()
df___seq.loc[:,"tail"] = df___seq.loc[:,2].str.slice(10,20).str.upper()
locs = df___seq[2].str.len() >= 20
for NT in ["A", "T", "C", "G"]:
    locs = locs & (df___seq["head"].str.count(NT) < 8) & (df___seq["tail"].str.count(NT) < 8)
print(sum(locs)/len(locs))

df___vcf = df___vcf.loc[df___seq.index[locs],:]
df___vcf.to_csv(PATH_out, sep="\t", header=False, index=False)