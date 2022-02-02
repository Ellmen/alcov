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


def print_mut_results(mut_results):
    for name in mut_results:
        muts, not_muts = mut_results[name]
        new_base = name[-1]
        print('{}:'.format(name))
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


def plot_mutations(sample_results, sample_names, min_depth):
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set_theme()
    names = sample_results[0].keys()
    sample_counts = [[mut_results[mut] for mut in names] for mut_results in sample_results]
    num_mutations = len(names)
    mut_fractions = [[] for _ in range(num_mutations)]
    for i in range(num_mutations):
        for counts in sample_counts:
            count = counts[i]
            total = count[0] + count[1]
            fraction = count[0]/total if total >= min_depth else -1
            mut_fractions[i].append(round(fraction, 2))
    no_reads = np.array([[f == -1 for f in fractions] for fractions in mut_fractions])
    ax = sns.heatmap(
        mut_fractions,
        annot=True,
        mask=no_reads,
        cmap=sns.cm.rocket_r,
        xticklabels=sample_names,
        yticklabels=names,
        vmin=0,
        vmax=1
    )
    plt.xlabel('Sample')
    plt.ylabel('Mutation')
    plt.show()


def find_mutants_in_bam(bam_path, mutations):
    import pysam

    samfile = pysam.Samfile(bam_path, "rb")

    # parsed_muts = [parse_snv(mut) for mut in mutations]
    parsed_muts = {}
    for mut in mutations:
        parsed_muts[mut] = parse_mutation(mut)
    mut_results = {mut: {snv_name(m): [0,0] for m in parsed_muts[mut]} for mut in mutations}

    for pileupcolumn in samfile.pileup():
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

    print_mut_results(mut_results)

    return mut_results


def mut_idx(mut):
    # Sort by genomic index of mutations
    snvs = parse_mutation(mut)
    if len(snvs) == 0:
        return -1
    return snvs[0][1]


# def find_mutants(file_path, mutations_path, min_depth, not_in): #TODO: not in lineage
def find_mutants(file_path, mutations_path, min_depth):
    """
    Accepts either a bam file or a tab delimited  txt file like
    s1.bam  Sample 1
    s2.bam  Sample 2
    """

    sample_results = []
    sample_names = []
    lineages = list(mut_lins['S:N501Y'].keys()) # arbitrary
    if mutations_path in lineages:
        lin = mutations_path
        print('Searcing for {} mutations'.format(lin))
        mutations = [mut for mut in mut_lins if mut_lins[mut][lin] > 0 and mut_idx(mut) != -1]
        # Unique
        mutations = [mut for mut in mutations if all(mut_lins[mut][l] == 0 for l in lineages if l != lin)]
        mutations.sort(key=mut_idx)
    else:
        print('Searching for mutations in {}'.format(mutations_path))
        with open(mutations_path, 'r') as f:
            mutations = [mut for mut in f.read().split('\n') if len(mut)]

    if file_path.endswith('.bam'):
        sample_results.append(find_mutants_in_bam(file_path, mutations))
        sample_names.append('')
    else:
        with open(file_path, 'r') as f:
            samples = [line.split('\t') for line in f.read().split('\n')]
        for sample in samples:
            if sample[0].endswith('.bam'): # Mostly for filtering empty
                print('{}:'.format(sample[1]))
                sample_results.append(find_mutants_in_bam(sample[0], mutations))
                sample_names.append(sample[1])
                print()
    plot_mutations(sample_results, sample_names, min_depth)
