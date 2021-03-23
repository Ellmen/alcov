from Bio import SeqIO

cov = list(SeqIO.parse("sequence.gb", "genbank"))[0]

genes = {}

for f in cov.features:
    if f.type == 'gene':
        gene = f.qualifiers['gene'][0]
        start = int(f.location.start)
        end = int(f.location.end)
        genes[gene] = [start, end]


with open('sars_cov_2.py', 'w') as f:
    f.write('genes = {}\n\nseq = \'{}\''.format(genes, cov.seq))
