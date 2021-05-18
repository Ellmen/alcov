from .convert_mutations import aa, nt
from .mutations import mutations as mut_lins


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
        total = muts + not_muts
        if total == 0:
            continue
        else:
            print('{} ({}):'.format(name, nt(name)))
            print('{} are {}, {} are wildtype ({:.2f}% of {} total)'.format(
                muts,
                new_base,
                not_muts,
                muts/total*100,
                total
            ))


def plot_lineages(sample_results, sample_names):
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set_theme()
    names = sample_results[0].keys()
    num_lineagess = len(names)
    lin_fractions = np.array([[lin_results[lin] for lin in names] for lin_results in sample_results]).T
    ax = sns.heatmap(
        lin_fractions,
        annot=True,
        # mask=no_reads,
        cmap=sns.cm.rocket_r,
        xticklabels=sample_names,
        yticklabels=names,
        vmin=0,
        vmax=1
    )
    plt.xlabel('Sample')
    plt.ylabel('Lineage')
    plt.tight_layout()
    plt.show()


def find_mutants_in_bam(bam_path):
    import numpy as np
    import pysam
    from sklearn.linear_model import LinearRegression
    from sklearn.linear_model import Ridge

    samfile = pysam.Samfile(bam_path, "rb")

    aa_mutations = [m for m in mut_lins.keys() if m[:5] not in ['ORF1a', 'ORF1b']] # TODO: parse orf1a/b
    # aa_mutations = [m for m in mut_lins.keys() if m[0] in ['S']] # Only spike
    aa_blacklist = ['S:D614G'] # all lineages contain this now
    aa_mutations = [m for m in aa_mutations if m not in aa_blacklist]
    mutations = parse_mutations(aa_mutations)
    lineages = list(mut_lins['S:N501Y'].keys()) # arbitrary

    parsed_muts = [parse_snv(mut) for mut in mutations]
    mut_results = {mut: [0,0] for mut in mutations}

    for pileupcolumn in samfile.pileup():
        pos = pileupcolumn.pos + 1
        for m in parsed_muts:
            if pos == m[1]:
                muts, not_muts = mut_in_col(pileupcolumn, m[2])
                mut_results['{}{}{}'.format(m[0], m[1], m[2])] = [muts, not_muts]
    samfile.close()
    covered_nt_muts = [m for m in mutations if sum(mut_results[m]) > 40]
    covered_muts = []
    for m in aa_mutations:
        nt_muts = aa(m)
        if any([nt_mut in covered_nt_muts for nt_mut in nt_muts]):
            # Pick only highest prevalence SNP if there are multiple
            most_common_snp = None
            snp_prevalence = -1
            for nt_mut in nt_muts:
                if mut_results[nt_mut][0] > snp_prevalence:
                    most_common_snp = nt_mut
                    snp_prevalence = mut_results[nt_mut][0]
            covered_muts.append(most_common_snp)
    covered_lineages = set()
    for m in covered_muts:
        aa_m = nt(m)
        for l in mut_lins[aa_m]:
            if mut_results[m][0] > 0 and mut_lins[aa_m][l] > 0.5:
                covered_lineages.add(l)
    covered_lineages = [l for l in lineages if l in covered_lineages]
    # TODO: filter uncovered lineages
    Y = np.array([mut_results[m][0]/sum(mut_results[m]) if sum(mut_results[m]) > 0 else 0 for m in covered_muts])
    lin_mut_profiles = [[round(mut_lins[nt(mut)][lin]) for mut in covered_muts] for lin in lineages]
    # Merge indistinguishable lineages
    merged_lmps = []
    merged_lins = []
    for i in range(len(lineages)):
        lmp = lin_mut_profiles[i]
        lin = lineages[i]
        if lmp in merged_lmps:
            lmp_idx = merged_lmps.index(lmp)
            merged_lins[lmp_idx] += ' or '
            merged_lins[lmp_idx] += lin
        else:
            merged_lmps.append(lmp)
            merged_lins.append(lin)
    # X = np.array([[round(mut_lins[nt(mut)][lin]) for lin in lineages] for mut in covered_muts])
    X = np.array(merged_lmps).T
    # print([nt(m) for m in covered_muts])
    # print(Y)
    reg = LinearRegression(fit_intercept=0, positive=True).fit(X, Y)

    print_mut_results(mut_results)
    # sample_results = {lineages[i]: round(reg.coef_[i], 3) for i in range(len(lineages))}
    sample_results = {merged_lins[i]: round(reg.coef_[i], 3) for i in range(len(merged_lins))}

    return sample_results


def find_lineages(file_path):
    """
    Accepts either a bam file or a tab delimited  txt file like
    s1.bam  Sample 1
    s2.bam  Sample 2
    """

    sample_results = []
    sample_names = []

    if file_path.endswith('.bam'):
        sample_results.append(find_mutants_in_bam(file_path))
        sample_names.append('')
    else:
        with open(file_path, 'r') as f:
            samples = [line.split('\t') for line in f.read().split('\n')]
        for sample in samples:
            if sample[0].endswith('.bam'): # Mostly for filtering empty
                print('{}:'.format(sample[1]))
                sample_results.append(find_mutants_in_bam(sample[0]))
                sample_names.append(sample[1])
                print()
    plot_lineages(sample_results, sample_names)
