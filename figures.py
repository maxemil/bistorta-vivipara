import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import upsetplot

sh = pd.read_csv('B_vivipara_matches_out_99_all_proc.csv', sep='\t', index_col=False)

sh['genbank'] =sh['seq_accno'].apply(lambda x: x.split('|')[0])
accessions = [line.strip() for line in open('accessions.txt')]
# sh = sh[sh['genbank'].apply(lambda x : x in accessions)]
# sh['region'] = sh['genbank'].apply(lambda x: 'Swabian alps' if x.startswith('KY') else 'Bavarian alps')
swalps = [x for x in accessions if x.startswith('KY')]
bvalps = [x for x in accessions if not x.startswith('KY')]
sh.loc[sh['genbank'].apply(lambda x: x in swalps), 'region'] = 'Swabian alps'
sh.loc[sh['genbank'].apply(lambda x: x in bvalps), 'region'] = 'Bavarian alps'
sh['region'] = sh['region'].fillna('Other')

sh.to_csv('Oberjoch_Alb_SH.csv', header=True, sep='\t', index=False)

order2sh = defaultdict(set)
for row in sh.iterrows():
    tax = row[1]['SH/compound taxonomy (1.0)']
    order = {t.split('__')[0]:t.split('__')[1] for t in tax.split(';')}['o']
    order2sh[order].add(row[1]['SH code (1.0)'])

sh2order = {s:t for t,sh in order2sh.items() for s in sh}
sh['order'] = sh['SH code (1.0)'].apply(lambda x: sh2order[x])

#number of sequences per order
sh.loc[sh['region'] == 'Bavarian alps', 'order'].value_counts()
sh.loc[sh['region'] == 'Swabian alps', 'order'].value_counts()
sh.loc[sh['region'] == 'Other', 'order'].value_counts()
#number of SHs per order
sh.loc[sh['region'] == 'Bavarian alps', ['SH code (1.0)', 'order']].drop_duplicates()['order'].value_counts()
sh.loc[sh['region'] == 'Swabian alps', ['SH code (1.0)', 'order']].drop_duplicates()['order'].value_counts()
sh.loc[sh['region'] == 'Other', ['SH code (1.0)', 'order']].drop_duplicates()['order'].value_counts()

shu = sh[['region', 'order', 'SH code (1.0)']].drop_duplicates()

contents = {'Other':shu.loc[shu['region'] == 'Other', 'SH code (1.0)'],
        'Swabian alps':shu.loc[shu['region'] == 'Swabian alps', 'SH code (1.0)'],
        'Bavarian alps':shu.loc[shu['region'] == 'Bavarian alps', 'SH code (1.0)']}

# fig, axs = plt.subplots(1, 2, figsize=(10,20))
fig = plt.figure(figsize=(15,15))

upset = upsetplot.from_contents(contents)
uplot = upsetplot.plot(upset)

plt.savefig('Upset_Bviv.pdf')

fig,ax = plt.subplots(figsize=(12,6))

shugroup = shu.groupby(['region', 'order']).size().reset_index().pivot(columns='order', index='region', values=0)
shugroup = shugroup.sort_values(axis=1, by='Other', ascending=False)
shugroup = shugroup.fillna(0)
shugroup = shugroup.div(shugroup.sum(axis=1), axis=0) * 100
colors = sns.color_palette("hls", shugroup.shape[1], as_cmap=True)
shugroup.plot(kind='bar', stacked=True, ax=ax, cmap=colors)

ax.set_ylabel("Relative abundance [%]")
ax.set_xlabel("Plot")
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.tight_layout()
plt.savefig('order_dist.pdf')



#### old stuff
counts = pd.DataFrame(index=order2sh.keys(), columns=['Bavarian alps', 'Shared', 'Swabian alps']).fillna(0)

for order, shs in order2sh.items():
    for s in shs:
        if len(set(sh.loc[sh['SH code (1.0)'] == s, 'region'])) == 1:
            counts.loc[order, set(sh.loc[sh['SH code (1.0)'] == s, 'region']).pop()] += 1
        else:
            counts.loc[order, 'Shared'] += 1

counts.sort_values(by='Bavarian alps', inplace=True, ascending=False)


# order2count = {k:len(v) for k,v in order2sh.items()}
# counts = pd.Series(order2count)
# counts.sort_values(ascending=False, inplace=True)
clrs_gray = ['gray' if x == 'Sebacinales' else 'darkgray' for x in counts.index]


#ax = sns.barplot(counts.index, counts.values, palette=clrs_gray)
ax = counts.plot(kind='bar', stacked=True)
ax.set_xlabel("Order", fontsize='x-large')
ax.set_ylabel("# of SHs", fontsize='x-large')
ax.tick_params(axis='y', labelsize='x-large')
plt.xticks(
    rotation=65, 
    horizontalalignment='right',
    fontweight='light',
    fontsize='x-large'
)
plt.tight_layout()
plt.savefig('Fungal_orders_oberjoch.pdf')
# plt.show()
