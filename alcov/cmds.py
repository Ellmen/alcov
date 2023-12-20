from .convert_mutations import aa as convert_aa, nt as convert_nt


def aa(mut):
    nt_muts = convert_aa(mut)
    for nt_mut in nt_muts:
        print('{} causes {}'.format(nt_mut, mut))


def nt(mut):
    aa_mut = convert_nt(mut)
    print('{} causes {}'.format(mut, aa_mut))

