import pandas as pd
import numpy as np
import seaborn as sns
import ete3
import matplotlib.pyplot as plt
from collections import defaultdict

ncbi = ete3.NCBITaxa()

df = pd.read_csv('Fallopia_all_counts_genus.csv', sep='\t')
df = df[['plot']]
df.columns = ['Genus']
df.to_csv("test.csv", index_label='OTU_ID', sep='\t')



ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
def taxdefault():
    return {s:'' for s in ranks}

all_tax = []
for t in df['Genus']:
    # t = df.loc['OTU0','Genus']
    taxdict = taxdefault()
    try:
        tid = ncbi.get_name_translator([t])[t]
        for l, lname in ncbi.get_taxid_translator(ncbi.get_lineage(tid[0])).items():
            rank = ncbi.get_rank([l])[l]
            if rank in ranks:
                taxdict[rank] = lname
    except:
        taxdict['genus'] = t
    all_tax.append(taxdict)


df = pd.DataFrame(all_tax)
df.index = pd.Series(df.index).apply(lambda x: "OTU{}".format(x))
df.columns = [f.capitalize() for f in df.columns]
df.to_csv("Fallopia_Taxonomy.csv", index_label='OTU_ID', sep='\t')

os.system('python /data/Software/FUNGuild/FUNGuild.py guild -taxa Fallopia_Taxonomy.csv')

gd = pd.read_csv('Fallopia_Taxonomy.guilds.txt', sep='\t')
gd['trophicMode'] = gd['trophicMode'].replace('Saportroph', 'Saprotroph')
gd['trophicMode'] = gd['trophicMode'].replace('Saprotroph; Symbiotroph', 'Saprotroph-Symbiotroph')


guilds = gd['trophicMode'].value_counts()
guilds.index = [i if not i == 'na' else 'Unassigned' for i in guilds.index]

df = pd.read_csv('Fallopia_all_counts_genus.csv', sep='\t')
df.index = df['plot']
df = df.drop('plot', axis=1)
df.rename(index={'Phialemoniopsis':'Thyridium'}, inplace=True)
df = df > 0

gd = gd[['Genus', 'trophicMode']]

df = df.loc[gd['Genus']]

gd.index = gd['Genus']
gd = pd.Series(gd['trophicMode'])
gd = gd.apply(lambda x: x if not x == 'na' else 'Unassigned')

dfg = df.mul(gd, axis=0)
dfga = dfg.apply(pd.Series.value_counts)
dfga['trophicMode'] = dfga.index
dm = pd.melt(dfga, id_vars=['trophicMode'])
dm['variable'] = dm['variable'].apply(lambda x: 'Vegetation' if 'Begleit' in x else x)
dm['variable'] = dm['variable'].apply(lambda x: 'R. japonica' if 'jap.' in x else x)
dm['variable'] = dm['variable'].apply(lambda x: 'R. bohemica' if 'boh.' in x else x)
dm['trophicMode'].replace('', np.nan, inplace=True)
dm = dm.dropna()
dm = dm.sort_values(by='value')

fig, ax = plt.subplots(figsize=(14,7))
# sns.set_theme(font_scale = 1.5)
g = sns.violinplot(y='trophicMode', x='value', data=dm, hue='variable', facet_kws={'legend_out': True})
leg = g.axes.get_legend()
leg.set_title('sample type')
new_labels = ['Vegetation (n=152)', 'R. japonica (n=104)', 'R. bohemica (n=48)']
for t, l in zip(leg.texts, new_labels):
    t.set_text(l)
g.set_xlabel('count')
g.set_ylabel('trophic mode')
plt.tight_layout()
plt.savefig('Fallopia_FunGuild.pdf')
