from .convert_mutations import aa, nt


b117_mutations = [
    'A28271-',
    'G28280C',
    'A28281T',
    'T28282A',
    'S:N501Y',
    # 'S:A570D',
    # 'S:E484K'
]


def parse_mutations(mutations):
    nts = [mut for mut in mutations if ':' not in mut]
    aas = [mut for mut in mutations if ':' in mut]
    return nts + sum([aa(mut) for mut in aas], [])


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


def print_mut_results(mut_results):
    for name in mut_results:
        muts, not_muts = mut_results[name]
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


def plot_mutations(sample_results, sample_names):
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set_theme()
    names = sample_results[0].keys()
    sample_counts = [[mut_results[mut] for mut in names] for mut_results in sample_results]
    num_mutations = len(names)
    mut_fractions = [[] for _ in range(num_mutations)]
    min_reads = 5
    for i in range(num_mutations):
        for counts in sample_counts:
            count = counts[i]
            total = count[0] + count[1]
            fraction = count[0]/total if total >= min_reads else -1
            mut_fractions[i].append(fraction)
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

    parsed_muts = [parse_snv(mut) for mut in mutations]
    mut_results = {mut: [0,0] for mut in mutations}

    for pileupcolumn in samfile.pileup():
        pos = pileupcolumn.pos + 1
        for m in parsed_muts:
            if pos == m[1]:
                muts, not_muts = mut_in_col(pileupcolumn, m[2])
                mut_results['{}{}{}'.format(m[0], m[1], m[2])] = [muts, not_muts]
    samfile.close()

    print_mut_results(mut_results)

    return mut_results


def find_mutants(file_path, mutations_path=None):
    """
    Accepts either a bam file or a tab delimited  txt file like
    s1.bam  Sample 1
    s2.bam  Sample 2
    """

    sample_results = []
    sample_names = []
    if mutations_path:
        print('Searching for mutations in {}'.format(mutations_path))
        with open(mutations_path, 'r') as f:
            mutations = parse_mutations([mut for mut in f.read().split('\n') if len(mut)])
    else:
        print('Searching for B.1.1.7 mutations...')
        mutations = parse_mutations(b117_mutations)

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
    plot_mutations(sample_results, sample_names)
