import re

from .sars_cov_2 import genes, seq


codons = { 
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
} 

nts = 'ACGT'


def aa(mut):
    gene = mut[:mut.find(':')]
    if gene == 'DEL':
        nt_idx, length = re.findall(r'\d+', mut)
        nt_idx = int(nt_idx)
        length = int(length)
        return ['{}{}-'.format(seq[nt_idx+i], nt_idx+i) for i in range(length)]
    aa_idx = int(re.findall(r'\d+', mut)[-1])
    nt_idx = genes[gene][0] + (aa_idx - 1) * 3
    codon = seq[nt_idx:nt_idx+3]
    if mut[mut.find(':')+1:].startswith('DEL') or mut[-1] == '-':
        return ['{}{}-'.format(seq[nt_idx],nt_idx+1) for nt_idx in range(nt_idx,nt_idx+3)]
    new_acid = mut[-1]
    acid = codons[codon]
    nt_muts = []

    for i in range(4):
        if codons[nts[i] + codon[1] + codon[2]] == new_acid:
            nt_muts.append('{}{}{}'.format(codon[0], nt_idx + 1, nts[i]))

    for i in range(4):
        if codons[codon[0] + nts[i] + codon[2]] == new_acid:
            nt_muts.append('{}{}{}'.format(codon[1], nt_idx + 2, nts[i]))

    for i in range(4):
        if codons[codon[0] + codon[1] + nts[i]] == new_acid:
            nt_muts.append('{}{}{}'.format(codon[2], nt_idx + 3, nts[i]))

    return nt_muts


def nt(mut):
    base = mut[0]
    new_base = mut[-1]
    nt_idx = int(re.findall(r'\d+', mut)[0]) - 1
    for gene in genes:
        l = genes[gene]
        if l[0] < nt_idx and l[1] > nt_idx:
            nt_offset = (nt_idx - l[0]) % 3
            # Slightly roundabout way of getting start of codon
            aa_idx = int((nt_idx - l[0]) / 3) + 1
            if new_base == '-':
                return '{}:DEL{}'.format(gene, aa_idx)
            codon_start = l[0] + (aa_idx - 1) * 3
            codon = list(seq[codon_start:codon_start+3])
            acid = codons[''.join(codon)]
            codon[nt_offset] = new_base
            new_acid = codons[''.join(codon)]
            return '{}:{}{}{}'.format(gene, acid, aa_idx, new_acid)

