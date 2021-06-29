from math import ceil, floor

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
    names = set()
    for sr in sample_results:
        for key in sr.keys():
            names.add(key)
    names = [n for n in names]
    names.sort()
    # names = sample_results[0].keys()
    lin_fractions = np.array([[lin_results[lin] if lin in lin_results else 0 for lin in names] for lin_results in sample_results]).T
    # fig, ax = plt.subplots(figsize=(len(sample_names)/2,len(names)/2))
    ax = sns.heatmap(
        lin_fractions,
        # annot=True,
        # mask=no_reads,
        cmap=sns.cm.rocket_r,
        xticklabels=sample_names,
        yticklabels=names,
        vmin=0,
        vmax=1,
        square=True,
        cbar=False,
    )
    plt.xlabel('Frequency in sample')
    plt.ylabel('SARS-CoV-2 lineage')
    # plt.subplots_adjust(bottom=0.3, left=0.3)
    # ax.figure.tight_layout()
    # mng = plt.get_current_fig_manager()
    # mng.frame.Maximize(True)
    plt.tight_layout()
    plt.show()


def plot_lineages_timeseries(sample_results, sample_names):
    from datetime import date
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set_theme()

    names = set()
    for sr in sample_results:
        for key in sr.keys():
            names.add(key)
    names = [n for n in names]
    names.sort()
    # names = list(sample_results[0].keys())
    # names.sort()
    if '_' in sample_names[0]:
        locations = [dt.split('_')[0] for dt in sample_names]
    # loc_data = {}
    locations = []
    dates = []
    for i in range(len(sample_names)):
        loc, dt = sample_names[i].split('_')
        dt = date.fromisoformat(dt)
        locations.append(loc)
        dates.append(dt)
    loc_set = sorted(list(set(locations)))
    ncols = 3
    nrows = ceil(len(loc_set) / ncols)
    fig, axes = plt.subplots(nrows, ncols)
    plt.tight_layout()
    non_zero = set()
    for sr in sample_results:
        for name in sr:
            if sr[name] > 0:
                non_zero.add(name)
    for sr in sample_results:
        for name in list(sr.keys()):
            if name not in non_zero:
                del sr[name]
    names = sorted(list(non_zero))
    for r in range(nrows):
        for c in range(ncols):
            i = r * ncols + c
            if i < len(loc_set):
                loc = loc_set[i]
                loc_results = []
                loc_dates = []
                for i in range(len(sample_results)):
                    sample_loc = locations[i]
                    if sample_loc == loc:
                        loc_results.append(sample_results[i])
                        loc_dates.append(dates[i])
                d = {name: [sr[name] if name in sr else 0 for sr in loc_results] for name in names}
                df = pd.DataFrame(data=d, index=loc_dates)
                # df = df.loc[:, df.any()] # Delete all-zero lineages
                axes[r,c].set_title(loc)
                axes[r,c].set_xlim([min(dates),max(dates)])
                axes[r,c].set_ylim([0,1])
                df.plot.area(ax=axes[r,c], fontsize=8, legend=False, rot=30)
            else:
                fig.delaxes(axes[r,c])
    leg = fig.legend(labels=names, loc="lower right", ncol=3)
    for line in leg.get_lines():
        line.set_linewidth(6)
    plt.show()


def show_lineage_predictions(sample_results, X, Y, covered_muts):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set_theme()
    merged_lins = list(sample_results.keys())
    d = {ml: [] for ml in merged_lins}
    d['Observed'] = []
    idx = []
    for i in range(len(covered_muts)):
        if Y[i] < 0.02:
            # Skip non-present mutations
            continue
        cm = nt(covered_muts[i])
        idx.append('{}\nobserved'.format(cm))
        idx.append('{}\npredicted'.format(cm))
        for j in range(len(merged_lins)):
            ml = merged_lins[j]
            lin_contribution = X[i][j]*sample_results[ml]
            d[ml] += [0, lin_contribution]
        d['Observed'] += [Y[i], 0]
    df = pd.DataFrame(data=d, index=idx)
    df = df.loc[:, df.any()] # Delete all-zero lineages
    df.plot.bar(stacked=True, rot=0)
    plt.show()


def do_regression(lmps, Y):
    # Perform linear regression and redo if sum(frequencies) > 1
    import numpy as np
    from sklearn.linear_model import LinearRegression
    from sklearn.linear_model import Lasso
    X = np.array(lmps).T
    reg = LinearRegression(fit_intercept=0, positive=True).fit(X, Y)
    score = reg.score(X, Y)
    best_score = score
    best_reg = reg
    if sum(reg.coef_) > 1:
        # If there is no valid solution, drop uncovered mutations one at a time
        # This allows for lineages with one missing mutation
        valid = False
        best_score = 0.8 # minimum valid score
        for i in range(len(X)):
            if Y[i] == 0:
                new_X = np.concatenate((X[:i], X[i+1:]))
                new_Y = np.concatenate((Y[:i], Y[i+1:]))
                new_reg = LinearRegression(fit_intercept=0, positive=True).fit(new_X, new_Y)
                new_score = new_reg.score(X, Y)
                if sum(new_reg.coef_) <= 1 and new_score > best_score:
                    valid = True
                    best_score = new_score
                    best_reg = new_reg
        if not valid:
            print('Warning: solutions sums to > 1')
    return X, best_reg


def find_mutants_in_bam(bam_path, return_data=False):
    import numpy as np
    import pysam

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
    covered_nt_muts = [m for m in mutations if sum(mut_results[m]) > 5]
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
    if len(covered_muts) == 0:
        print('No coverage')
        return None
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
    X, reg = do_regression(merged_lmps, Y)

    print_mut_results(mut_results)
    sample_results = {merged_lins[i]: round(reg.coef_[i], 2) for i in range(len(merged_lins))}

    if return_data:
        return sample_results, X, Y, covered_muts
    return sample_results


def find_lineages(file_path, ts=False, csv=False):
    """
    Accepts either a bam file or a tab delimited  txt file like
    s1.bam  Sample 1
    s2.bam  Sample 2
    """

    sample_results = []
    sample_names = []

    if file_path.endswith('.bam'):
        sr, X, Y, covered_muts = find_mutants_in_bam(file_path, return_data=True)
        show_lineage_predictions(sr, X, Y, covered_muts)
        sample_results.append(sr)
        sample_names.append('')
    else:
        with open(file_path, 'r') as f:
            samples = [line.split('\t') for line in f.read().split('\n')]
        for sample in samples:
            if sample[0].endswith('.bam'): # Mostly for filtering empty
                print('{}:'.format(sample[1]))
                sample_result = find_mutants_in_bam(sample[0])
                if sample_result is not None:
                    sample_results.append(sample_result)
                    sample_names.append(sample[1])
                print()
    if ts:
        plot_lineages_timeseries(sample_results, sample_names)
    else:
        plot_lineages(sample_results, sample_names)
