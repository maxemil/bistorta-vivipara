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

numall = []
numbistorta = []
for sh in set(df['SH code (1.0)']):
    hosts = set(df.loc[df['SH code (1.0)'] == sh, 'host'])
    numall.append(len(hosts))
    if 'Bistorta_vivipara' in hosts:
        numbistorta.append(len(hosts))


sh_bistorta = set()
for sh in set(df['SH code (1.0)']):
    hosts = set(df.loc[df['SH code (1.0)'] == sh, 'host'])
    if 'Bistorta_vivipara' in hosts:
        sh_bistorta.add(sh)


df_bis = df[df['SH code (1.0)'].apply(lambda x: x in sh_bistorta)]
df_bis['host'] = df_bis['host'].apply(lambda x: x.split('_')[0])
df_bis['host'] = df_bis['host'].apply(lambda x: get_family(x))
df_bis = df_bis.dropna()

df_bis = df_bis[df_bis['host'].apply(lambda x: x not in ['', 'x'])]
df_bis.loc[df_bis['host'] == 'orchid', 'host'] = 'Orchidaceae'

sh2host = {}
for sh in set(df_bis['SH code (1.0)']):
    hosts = len(set(df_bis.loc[df_bis['SH code (1.0)'] == sh, 'host']))
    sh2host[sh] = hosts    

sel_fam = ['Polygonaceae', 'Pinaceae', 'Fagaceae', 'Orchidaceae', 'Ericaceae',
'Betulaceae', 'Salicaceae']
# df_bis = df_bis[df_bis['host'].apply(lambda x: x in sel_fam)]

num_otu_overlap = []
for x,y in combinations(set(df_bis['host']), 2):
    s1 = set(df_bis.loc[df_bis['host'] == x, 'SH code (1.0)'])
    s2 = set(df_bis.loc[df_bis['host'] == y, 'SH code (1.0)'])
    num_overlap = len(s1.intersection(s2))
    if num_overlap > 0:
        num_otu_overlap.append((x,y,num_overlap))

df_num = pd.DataFrame(num_otu_overlap)
df_num.columns = ['source', 'target', 'weight']
df_num.to_csv('OTU_overlaps.csv', header=True, sep='\t', index=False)

G = nx.from_pandas_edgelist(df_num, source='source', target='target', edge_attr='weight')
pos=nx.circular_layout(G)
edges = G.edges()
weights = [G[u][v]['weight']/2 for u,v in edges]

SHnum = {}
for h in set(df_bis['host']):
    SHnum[h] = len(set(df_bis.loc[df_bis['host'] == h, 'SH code (1.0)']))
node_sizes = [SHnum[g]*20 for g in G.nodes()]

plt.figure(figsize=(10,10))
nx.draw(G, pos=pos, with_labels=True, width=weights, verticalalignment='top', node_size=node_sizes)
plt.savefig('circular_network.pdf')
# df_bis.loc[df_bis['host'].apply(lambda x: x not in sel_fam), 'host'] = 'others'

df_bis = df_bis[df_bis['host'] != 'Polygonaceae']

host_sh_count = df_bis.groupby(['host', 'SH code (1.0)']).size().reset_index().rename(columns={0:'count'})
min_host = host_sh_count['host'].value_counts()[host_sh_count['host'].value_counts() > 5].index
df_bis = df_bis[df_bis['host'].apply(lambda x: x in min_host)]

G = nx.from_pandas_edgelist(df_bis, source='host', target='SH code (1.0)')

cmap = sns.color_palette("YlOrBr", as_cmap=True)
colors = []
for n in G:
    if n in sh2host:
        colors.append(cmap(sh2host[n]/14))
    else:
        colors.append('blue')

top = nx.bipartite.sets(G)[0]
pos = nx.bipartite_layout(G, top)

fig, axs = plt.subplots(2, figsize=(10,16), gridspec_kw={'height_ratios': [20, 1]})
p = nx.draw(G, ax=axs[0], pos=pos, with_labels=True, node_color=colors)


cb = mpl.colorbar.ColorbarBase(ax=axs[1], orientation='horizontal', 
                               cmap=cmap, norm=mpl.colors.Normalize(0, 14))

# plt.tight_layout()
plt.savefig('Bistorta_sebacinales_network.pdf', dpi=600)
