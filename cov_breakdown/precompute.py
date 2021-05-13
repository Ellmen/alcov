import json

from Bio import SeqIO


def process_reference():
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


def get_amplicons():
    cov = list(SeqIO.parse("sequence.gb", "genbank"))[0]
    with open('nCoV-2019.insert.bed', 'r') as f:
        inserts = [l.split('\t') for l in f.read().split('\n')[:-1]]

    for i in range(len(inserts)):
        insert = inserts[i]
        section = cov.seq[int(insert[1]):int(insert[2])]
        gc = (section.count('G') + section.count('C')) / len(section)
        inserts[i].append(gc)

    with open('artic_amplicons.py', 'w') as f:
        f.write('inserts = {}'.format(inserts))


def get_mutations():
   # TODO figure out capitalization
   with open('mutations.json', 'r') as f:
        raw_mutations = json.loads(f.read())

   muts = list(set([m['mutation'].upper() for m in raw_mutations]))
   lins = list(set([m['pangolin_lineage'].upper() for m in raw_mutations]))
   mut_lins = {mut: {lin: 0 for lin in lins} for mut in muts}
   for raw_m in raw_mutations:
       mut = raw_m['mutation'].upper()
       lin = raw_m['pangolin_lineage'].upper()
       prev = raw_m['prevalence']
       mut_lins[mut][lin] = prev

   # print(mut_lins['S:N501Y'])
   with open('mutations.py', 'w') as f:
        f.write('mutations = {}'.format(mut_lins))


if __name__ == '__main__':
    # process_reference()
    # get_amplicons()
    get_mutations()
