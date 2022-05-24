from .artic_amplicons import primers
from .convert_mutations import aa, nt
from .mutations import mutations as mut_lins


def parse_snv(snv):
    pos = int(snv[1:-1])
    old_bp = snv[0]
    new_bp = snv[-1]
    return old_bp, pos, new_bp


def snv_name(snv):
    return '{}{}{}'.format(*snv)

def parse_mutation(mut):
    if ':' in mut:
        muts = aa(mut)
    else:
        muts = [mut]
    return [parse_snv(m) for m in muts]


def get_lins(mut_lins):
    some_mut = next(iter(mut_lins.keys()))
    return list(mut_lins[some_mut].keys())


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
                not_muts += 1
                continue
            base = pileupread.alignment.query_sequence[qpos]
            if base == mut:
                muts += 1
            else:
                not_muts += 1
    return muts, not_muts


def print_mut_results(mut_results, min_depth):
    cov = 0
    mut_cov = 0
    for name in mut_results:
        muts, not_muts = mut_results[name]
        # new_base = name[-1]
        # print('{}:'.format(name))
        total = muts + not_muts
        if total >= min_depth:
            cov += 1
            if muts > 0:
                mut_cov += 1

    # for name in mut_results:
    #     muts, not_muts = mut_results[name]
    #     new_base = name[-1]
    #     print('{}:'.format(name))
    #     total = muts + not_muts
    #     if total == 0:
    #         print('No coverage of {}'.format(name))
    #     else:
    #         print('{} are {}, {} are wildtype ({:.2f}% of {} total)'.format(
    #             muts,
    #             new_base,
    #             not_muts,
    #             muts/total*100,
    #             total
    #         ))

    print('{}/{} mutations covered'.format(cov, len(mut_results)))
    print('{}/{} mutations detected'.format(mut_cov, len(mut_results)))

def plot_mutations(sample_results, sample_names, min_depth, img_path=None):
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set_theme()
    names = list(sample_results[0].keys())
    names.sort(key=mut_idx)
    sample_counts = [[mut_results[mut] for mut in names] for mut_results in sample_results]
    num_mutations = len(names)
    mut_fractions = [[] for _ in range(num_mutations)]
    for i in range(num_mutations):
        for counts in sample_counts:
            count = counts[i]
            total = count[0] + count[1]
            fraction = count[0]/total if total >= min_depth else -1
            # mut_fractions[i].append(round(fraction, 2))
            mut_fractions[i].append(round(fraction, 4))
    no_reads = np.array([[f == -1 for f in fractions] for fractions in mut_fractions])
    ax = sns.heatmap(
        mut_fractions,
        annot=True,
        mask=no_reads,
        cmap=sns.cm.rocket_r,
        xticklabels=sample_names,
        yticklabels=names,
        vmin=0,
        vmax=1,
    )
    plt.xlabel('Sample')
    plt.xticks(rotation=30)
    plt.ylabel('Mutation')
    # plt.show()
    plt.tight_layout()
    if img_path is not None:
        plt.savefig(img_path, dpi=300)
    else:
        plt.show()


def find_mutants_in_bam(bam_path, mutations):
    import pysam

    samfile = pysam.Samfile(bam_path, "rb")

    # parsed_muts = [parse_snv(mut) for mut in mutations]
    parsed_muts = {}
    for mut in mutations:
        parsed_muts[mut] = parse_mutation(mut)
    mut_results = {mut: {snv_name(m): [0,0] for m in parsed_muts[mut]} for mut in mutations}

    for pileupcolumn in samfile.pileup(stepper='nofilter'):
        pos = pileupcolumn.pos + 1
        for mut in parsed_muts:
            for m in parsed_muts[mut]:
                if pos == m[1]:
                    muts, not_muts = mut_in_col(pileupcolumn, m[2])
                    mut_results[mut][snv_name(m)] = [muts, not_muts]
    samfile.close()

    for mut in mut_results:
        max_freq = -1
        max_muts = [0,0]
        for m in mut_results[mut]:
            muts, not_muts = mut_results[mut][m]
            freq = muts/(muts+not_muts) if muts + not_muts > 0 else 0
            if freq > max_freq:
                max_freq = freq
                max_muts = [muts, not_muts]
        mut_results[mut] = max_muts

    return mut_results


def in_primer(pos):
    for primer in primers:
        start = int(primer[1])
        end = int(primer[2])
        if start <= pos and end >= pos:
            return True
    return False


def find_mutants_in_grom(grom_path, aa_mutations, filter_primers=True):
    with open(grom_path, 'r') as f:
        mut_lines = [line.split(',') for line in f.readlines()[1:]]
    cov_path = grom_path.replace('.mapped.csv', '.coverage.csv')
    with open(cov_path, 'r') as f:
        cov_lines = [line.split(',') for line in f.readlines()[1:]]
    loc_cov = [int(cl[1]) for cl in cov_lines]
    mut_results = {}
    for line in mut_lines:
        pos = int(line[0])
        if in_primer(pos):
            continue
        mut = line[2]
        mut_name = '{}@{}'.format(mut, pos)
        if mut_name not in aa_mutations:
            continue
        freq = float(line[3])
        # if freq < 0.05:
        #     freq = 0
        coverage = int(line[4])
        muts = round(freq * coverage)
        not_muts = coverage - muts
        mut_results[mut_name] = [muts, not_muts]

    # aa_mutations = [m for m in mut_lins.keys()]
    for mut_name in aa_mutations:
        if mut_name not in mut_results:
            mut, pos = mut_name.split('@')
            pos = int(pos)
            coverage = loc_cov[pos]
            mut_results[mut_name] = [0,coverage]

    return mut_results


def mut_idx(mut):
    # Sort by genomic index of mutations
    # snvs = parse_mutation(mut)
    # if len(snvs) == 0:
    #     return -1
    # return snvs[0][1]
    return int(mut.split('@')[1])


# def find_mutants(file_path, mutations_path, min_depth, not_in): #TODO: not in lineage
def find_mutants(file_path, mutations_path, min_depth, save_img):
    """
    Accepts either a bam file or a tab delimited  txt file like
    s1.bam  Sample 1
    s2.bam  Sample 2
    """

    sample_results = []
    sample_names = []
    # lineages = list(mut_lins['S:N501Y'].keys()) # arbitrary
    lineages = get_lins(mut_lins)
    # lineages = list(mut_lins['aa:S:N501Y'].keys()) # arbitrary
    # lineages = ['EPI_ISL_7190366', 'EPI_ISL_7877191', 'EPI_ISL_2793160', 'EPI_ISL_10389336']
    # print(lineages)
    if mutations_path in lineages:
        lin = mutations_path
        print('Searcing for {} mutations'.format(lin))
        # mutations = [mut for mut in mut_lins if mut_lins[mut][lin] > 0 and mut_idx(mut) != -1]
        mutations = [mut for mut in mut_lins if mut_lins[mut][lin] > 0]
        # Unique
        # mutations = [mut for mut in mutations if all(mut_lins[mut][l] == 0 for l in lineages if l != lin)]
        # mutations = [mut for mut in mutations if sum(mut_lins[mut][l] for l in lineages) == 1]
        mutations.sort(key=mut_idx)
    else:
        print('Searching for mutations in {}'.format(mutations_path))
        with open(mutations_path, 'r') as f:
            mutations = [mut for mut in f.read().split('\n') if len(mut)]

    if file_path.endswith('.bam'):
        sample_results.append(find_mutants_in_bam(file_path, mutations))
        sample_names.append('')
    if file_path.endswith('.csv'):
        sample_results.append(find_mutants_in_grom(file_path, mutations))
        sample_names.append('')
    else:
        with open(file_path, 'r') as f:
            samples = [line.split('\t') for line in f.read().split('\n')]
        for sample in samples:
            if sample[0].endswith('.bam'): # Mostly for filtering empty
                print('{}:'.format(sample[1]))
                sample_results.append(find_mutants_in_bam(sample[0], mutations))
                sample_names.append(sample[1])
                print_mut_results(sample_results[-1], min_depth)
                print()

    mutants_name = mutations_path.replace('.txt', '').replace('.', '')
    img_path = file_path.replace('.txt', '_{}_mutants.png'.format(mutants_name)) if save_img else None
    plot_mutations(sample_results, sample_names, min_depth, img_path)
