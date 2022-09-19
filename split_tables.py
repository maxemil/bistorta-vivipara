with open('Table_S3_split.csv', 'w') as out:
    for line in open('Table_S3.csv'):
        line = line.strip('\n').split('\t')
        line[1] = line[1].replace('|', '\t')
        print("\t".join(line), file=out)
with open('Supplementary_Item1_split.csv', 'w') as out:
    for line in open('Supplementary_Item1.csv'):
        line = line.strip('\n').split('\t')
        line[1] = line[1].replace('|', '\t')
        print("\t".join(line), file=out)
