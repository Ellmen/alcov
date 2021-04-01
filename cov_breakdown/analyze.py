import pysam

from .convert_mutations import aa, nt


b117_mutations = [
    'A28271-',
]

b117_aa_mutations = [
    'S:N501Y',
    'S:A570D'
]


def parse_snv(snv):
    pos = int(snv[1:-1])
    old_bp = snv[0]
    new_bp = snv[-1]
    return old_bp, pos, new_bp


def mut_in_col(pileupcolumn, mut):
    muts = 0
    not_muts = 0
    if mut == '-':
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del:
                muts += 1
            else:
                not_muts += 1
    else:
        for pileupread in pileupcolumn.pileups:
            qpos = pileupread.query_position
            if qpos is None:
                continue
            base = pileupread.alignment.query_sequence[qpos]
            if base == mut:
                muts += 1
            else:
                not_muts += 1
    return muts, not_muts


def print_mut_result(mut_result):
    name, muts, not_muts = mut_result
    new_base = name[-1]
    if new_base == '-':
        print('{}:'.format(name))
    else:
        print('{} ({}):'.format(name, nt(name)))
    total = muts + not_muts
    if total == 0:
        print('No coverage of {}'.format(name))
    else:
        print('{} are {}, {} are wildtype ({:.2f}% of {} total)'.format(
            muts,
            new_base,
            not_muts,
            muts/total*100,
            total
        ))


def find_mutants(bam_path):
    mutations = b117_mutations + sum([aa(mut) for mut in b117_aa_mutations], [])
    samfile = pysam.Samfile(bam_path, "rb")

    parsed_muts = [parse_snv(mut) for mut in mutations]
    mut_results = [] # TODO: fill with 0s so uncovered mutations are still printed

    for pileupcolumn in samfile.pileup():
        pos = pileupcolumn.pos + 1
        for m in parsed_muts:
            if pos == m[1]:
                muts, not_muts = mut_in_col(pileupcolumn, m[2])
                mut_results.append(['{}{}{}'.format(m[0], m[1], m[2]), muts, not_muts])
    samfile.close()

    for result in mut_results:
        print_mut_result(result)
