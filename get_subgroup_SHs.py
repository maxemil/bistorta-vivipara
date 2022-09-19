import pandas as pd


# df = pd.read_csv('matches_out_99.csv', sep='\t', index_col=False)

sebacinaceae = [line for line in open('Sebacinaceae')]
serendipidaceae = [line for line in open('Serendipidaceae')]
# acc_ser = set()
# for row in df.iterrows():
#     for s in serendipidaceae:
#         if s.startswith(row[1]['seq_accno'].split('|')[0]): 
#             acc_ser.add(row[1]['SH code (1.0)'])
#             serendipidaceae.remove(s)

ser_all = []
seb_all = []
count = 0
failed = 0
for line in open("mapping_to_initial_seqs.txt"):
    if count > 1154:
        break
    line = line.strip('\n').split('\t')
    if any([line[0] in s for s in serendipidaceae]):
        ser_all += line[1].split()
    elif any([line[0] in s for s in sebacinaceae]):
        seb_all += line[1].split()
    else:
        failed += 1
        print(line)
    count += 1

df = pd.read_csv('Supplementary Item 1 split.csv', sep='\t', index_col=False)

ser_shs = set()
for s in ser_all:
    try:
        ser_shs.add(df.loc[df['seq_accno'] == s, 'SH code (1.0)'].item())
    except:
        print(s)
seb_shs = set()
for s in seb_all:
    try:
        seb_shs.add(df.loc[df['seq_accno'] == s, 'SH code (1.0)'].item())
    except:
        print(s)
