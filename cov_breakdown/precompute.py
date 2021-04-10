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
    with open('nCoV-2019.insert.bed', 'r') as f:
        inserts = [l.split('\t') for l in f.read().split('\n')[:-1]]

    with open('artic_amplicons.py', 'w') as f:
        f.write('inserts = {}'.format(inserts))


if __name__ == '__main__':
    # process_reference()
    get_amplicons()
