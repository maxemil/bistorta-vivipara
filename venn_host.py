import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import ete3

ncbi = ete3.ncbi_taxonomy.NCBITaxa()

def get_family(name):
    try:
        taxid = ncbi.get_name_translator([name])[name][0]
        for l in ncbi.get_lineage_translator(ncbi.get_lineage(taxid)):
            if ncbi.get_rank([l])[l] == 'family':
                return ncbi.get_taxid_translator([l])[l]
    except:
        return name

df = pd.read_csv('alpine_sequences_host_geo.csv', sep='\t', header=0)
df['host'] = df['host'].apply(lambda x: x.split('_')[0])
df['host'] = df['host'].apply(lambda x: get_family(x))

s = df['SH code (1.0)'].value_counts()
df = df[df.isin(s.index[s >= 2]).values]

bistorta = set(df.loc[df['host'] == 'Polygonaceae', 'SH code (1.0)'])
other = set(df.loc[df['host'] != 'Polygonaceae', 'SH code (1.0)'])

plt.figure(figsize=(4,4))
venn2([bistorta, other], ('Bistorta vivipara', 'Other hosts'), normalize_to=0.5)
plt.title("SH overlap")
plt.tight_layout()
plt.savefig('SH_overlap_Bistorta.pdf')


df = df[df['host'] != 'Polygonaceae']

genus = 'Ericaceae'

bistorta = set(df.loc[df['host'] == genus, 'SH code (1.0)'])
other = set(df.loc[df['host'] != genus, 'SH code (1.0)'])

plt.figure(figsize=(4,4))
venn2([bistorta, other], (genus, 'Other hosts'), normalize_to=0.5)
plt.title("SH overlap")
plt.tight_layout()
plt.savefig('SH_overlap_{}.pdf'.format(genus))
