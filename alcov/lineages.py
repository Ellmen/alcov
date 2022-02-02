from math import ceil, floor

from .analyze import find_mutants_in_bam, print_mut_results
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


def write_csv(sample_results, sample_names):
    lin_names = set()
    for sr in sample_results:
        for key in sr.keys():
            lin_names.add(key)
    lin_names = [n for n in lin_names]
    lin_names.sort()
    csv_headers = ['Sample name'] + [n + ' %' for n in lin_names]
    csv_rows = []
    for i in range(len(sample_names)):
        sr = sample_results[i]
        csv_rows.append([sample_names[i]] + [str(round(sr[n], 3)) if n in sr else '0' for n in lin_names])
    with open('sample_lineages.csv', 'w') as f:
        f.write('\n'.join(','.join(row) for row in [csv_headers] + csv_rows))


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
        annot=True,
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
        Y_p = 0
        for j in range(len(merged_lins)):
            ml = merged_lins[j]
            lin_contribution = X[i][j]*sample_results[ml]
            Y_p += lin_contribution
        if Y[i] < 0.02 and Y_p < 0.02:
            # Skip non-present mutations
            # pass
            continue
        cm = covered_muts[i]
        idx.append('{}\nobserved'.format(cm))
        idx.append('{}\npredicted'.format(cm))
        for j in range(len(merged_lins)):
            ml = merged_lins[j]
            lin_contribution = X[i][j]*sample_results[ml]
            d[ml] += [0, lin_contribution]
        d['Observed'] += [Y[i], 0]
    df = pd.DataFrame(data=d, index=idx)
    df = df.loc[:, df.any()] # Delete all-zero lineages
    if len(df) == 0:
        print('No mutations to show')
        return
    df.plot.bar(stacked=True, rot=0)
    plt.tight_layout()
    plt.show()


def show_lineage_pie(sample_results):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set_theme()
    merged_lins = list(sample_results.keys())
    freqs = [sample_results[ml] for ml in merged_lins]
    if sum(freqs) < 1:
        # merged_lins.append('Other')
        # freqs.append(1 - sum(freqs))
        freqs = [f/sum(freqs) for f in freqs]
    df = pd.DataFrame(data={'Fraction': freqs}, index=merged_lins)
    df = df[df.Fraction > 0] # Delete all-zero lineages
    print(df)
    print(sample_results)
    df.plot.pie(y='Fraction', legend=False, autopct='%1.1f%%', ylabel='')
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


def find_lineages_in_bam(bam_path, return_data=False, min_depth=40, lineages=[], unique=True):
    import numpy as np
    import pysam

    samfile = pysam.Samfile(bam_path, "rb")

    aa_mutations = [m for m in mut_lins.keys()]
    # aa_mutations = [m for m in mut_lins.keys() if m[0] in ['S']] # Only spike
    # aa_mutations = [m for m in mut_lins.keys() if m[0] in ['N']] # Only N
    aa_blacklist = ['S:D614G'] # all lineages contain this now
    aa_mutations = [m for m in aa_mutations if m not in aa_blacklist]
    if len(lineages) == 0:
        lineages = list(mut_lins['S:N501Y'].keys()) # arbitrary
    if unique:
        aa_mutations = [mut for mut in aa_mutations if sum(mut_lins[mut][l] for l in lineages) == 1]
    mutations = parse_mutations(aa_mutations)
    vocs = ['B.1.1.7', 'B.1.617.2', 'P.1', 'B.1.351']
    vois = ['B.1.525', 'B.1.526', 'B.1.617.1', 'C.37']
    # if only_vocs:
    # lineages = vocs + vois
    # lineages = ['Omicron', 'BA.2', 'Delta']

    mut_results = find_mutants_in_bam(bam_path, aa_mutations)

    covered_muts = [m for m in aa_mutations if sum(mut_results[m]) > min_depth]
    if len(covered_muts) == 0:
        print('No coverage')
        return None
    covered_lineages = set()
    for m in covered_muts:
        for l in mut_lins[m]:
            if mut_results[m][0] > 0 and mut_lins[m][l] > 0.5:
                covered_lineages.add(l)
    covered_lineages = [l for l in lineages if l in covered_lineages]
    Y = np.array([mut_results[m][0]/sum(mut_results[m]) if sum(mut_results[m]) > 0 else 0 for m in covered_muts])
    lin_mut_profiles = [[round(mut_lins[mut][lin]) for mut in covered_muts] for lin in lineages]
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
    sample_results = {merged_lins[i]: round(reg.coef_[i], 3) for i in range(len(merged_lins))}
    # Normalize frequencies
    freq_sum = sum(sample_results[lin] for lin in sample_results)
    for lin in sample_results:
        sample_results[lin] = sample_results[lin] / freq_sum

    if return_data:
        return sample_results, X, Y, covered_muts
    return sample_results


def find_lineages(file_path, lineages_path, ts, csv, min_depth, show_stacked, unique):
    """
    Accepts either a bam file or a tab delimited  txt file like
    s1.bam  Sample 1
    s2.bam  Sample 2
    """

    sample_results = []
    sample_names = []

    lineages = []
    if lineages_path is not None:
        with open(lineages_path, 'r') as f:
            lineages = f.read().splitlines()
    if file_path.endswith('.bam'):
        sr, X, Y, covered_muts = find_lineages_in_bam(file_path, True, min_depth, lineages, unique)
        if show_stacked:
            show_lineage_predictions(sr, X, Y, covered_muts)
            show_lineage_pie(sr)
        sample_results.append(sr)
        sample_names.append('')
    else:
        with open(file_path, 'r') as f:
            samples = [line.split('\t') for line in f.read().split('\n')]
        for sample in samples:
            if sample[0].endswith('.bam'): # Mostly for filtering empty
                print('{}:'.format(sample[1]))
                sample_result = find_lineages_in_bam(sample[0], False, min_depth, lineages, unique)
                if sample_result is not None:
                    sample_results.append(sample_result)
                    sample_names.append(sample[1])
                print()
    if ts:
        plot_lineages_timeseries(sample_results, sample_names)
    else:
        plot_lineages(sample_results, sample_names)
    if csv:
        write_csv(sample_results, sample_names)
