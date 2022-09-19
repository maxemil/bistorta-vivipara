import seaborn as sns
import pandas as pd
import ete3
import networkx as nx
from matplotlib import pyplot as plt
from itertools import combinations

ncbi = ete3.ncbi_taxonomy.NCBITaxa()

def get_family(name):
    try:
        taxid = ncbi.get_name_translator([name])[name][0]
        for l in ncbi.get_lineage_translator(ncbi.get_lineage(taxid)):
            if ncbi.get_rank([l])[l] == 'family':
                return ncbi.get_taxid_translator([l])[l]
    except:
        return name

df = pd.read_csv('Supplementary_Item1_split.csv', sep='\t', header=0, dtype={'host': str})
df = df.fillna('')


df['host'] = df['host'].apply(lambda x: x.split('_')[0])
df['host'] = df['host'].apply(lambda x: get_family(x))
df = df.dropna()


num_otu_overlap = []
for x,y in combinations(set(df['host']), 2):
    s1 = set(df.loc[df['host'] == x, 'SH code (1.0)'])
    s2 = set(df.loc[df['host'] == y, 'SH code (1.0)'])
    num_overlap = len(s1.intersection(s2))
    if num_overlap > 0:
        num_otu_overlap.append((x,y,num_overlap))

df_num = pd.DataFrame(num_otu_overlap)
df_num.columns = ['source', 'target', 'weight']
df_num.to_csv('OTU_overlaps.csv', header=True, sep='\t', index=False)
