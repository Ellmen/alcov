from .artic_amplicons import inserts
from .convert_mutations import aa, nt


inserts = inserts[1:98]

def plot_depths(sample_results, sample_names):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set_theme()
    import seaborn as sns
    # sns.set_theme(style="whitegrid")
    # samples = sum([98*[name] for name in sample_names], [])
    # samples = [s[71:85] for s in samples]
    samples = sum([len(inserts)*[name] for name in sample_names], [])
    amplicons = sum([list(amplicon.keys()) for amplicon in sample_results], [])
    # pools = ['Pool 1' if int(a) % 2 == 0 else 'Pool 2' for a in amplicons]
    pools = ['Pool 1' if True else 'Pool 2' for a in amplicons]
    depths = sum([[np.log(max(a, 1)) for a in amplicon.values()] for amplicon in sample_results], [])
    d = {
        'Sample': samples,
        'Amplicon number': amplicons,
        'Pool': pools,
        'Log depth': depths
    }
    df = pd.DataFrame(data=d)
    g = sns.FacetGrid(df, row="Sample", hue="Pool", height=1.7, aspect=8)
    g.map(sns.barplot, "Amplicon number", "Log depth", order=[str(i) for i in range(1,99)], hue_order=['Pool 1', 'Pool 2'])
    # plt.locator_params(axis='x', nbins=20)
    plt.locator_params(axis='x')
    plt.savefig('SA_coverage.png')
    plt.show()


def plot_depths_gc(sample_results, sample_names):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set_theme()
    import seaborn as sns
    # sns.set_theme(style="whitegrid")
    depths = sample_results[0]
    samples = sum([98*[name] for name in sample_names], [])
    amplicons = sum([list(amplicon.keys()) for amplicon in sample_results], [])
    gcs = [inserts[int(a)-1][6] for a in amplicons]
    pools = ['Pool 1' if int(a) % 2 == 0 else 'Pool 2' for a in amplicons]
    depths = sum([[np.log(max(a, 1)) for a in amplicon.values()] for amplicon in sample_results], [])
    d = {
        'Sample': samples,
        'Amplicon number': [int(a) for a in amplicons],
        'GC content': gcs,
        'Pool': pools,
        'Log depth': depths
    }
    df = pd.DataFrame(data=d)
    df = df.loc[df['Log depth'] != 0]
    g = sns.FacetGrid(df, col="Sample")
    g.map(sns.regplot, "GC content", "Log depth")
    # g.map(sns.regplot, "Amplicon number", "Log depth")
    plt.show()


def plot_amplified_fraction(sample_results, sample_names):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set_theme()
    import seaborn as sns
    # sns.set_theme(style="whitegrid")
    depths = sample_results[0]
    d = {
        'Sample': sum([98*[name] for name in sample_names], []),
        'Amplicon': sum([list(amplicon.keys()) for amplicon in sample_results], []),
        'Depth': sum([[min(d, 1) for d in amplicon.values()] for amplicon in sample_results], [])
    }
    df = pd.DataFrame(data=d)
    g = sns.FacetGrid(df, row="Sample")
    g.map(sns.barplot, "Amplicon", "Depth")
    plt.locator_params(axis='x', nbins=20)
    plt.show()


def find_depths_in_bam(bam_path, max_depth=50000):
    import pysam

    samfile = pysam.Samfile(bam_path, "rb")

    amp_mids = {int((int(i[1]) + int(i[2])) / 2): i[3] for i in inserts}
    amplified = {i[3]: 0 for i in inserts}

    for pileupcolumn in samfile.pileup(max_depth=max_depth):
        pos = pileupcolumn.pos
        if pos in amp_mids:
            depth = pileupcolumn.get_num_aligned()
            amplified[amp_mids[pos]] = depth
    samfile.close()

    # SAMPLE_NUM = '2'

    # with open('s{}_mixed.csv'.format(SAMPLE_NUM), 'w') as f:
    #     f.write('\n'.join(['{},{}'.format(insert_name, amplified[insert_name]) for insert_name in amplified.keys()]))

    return amplified


def amplicon_coverage(file_path):
    """
    Accepts either a bam file or a tab delimited  txt file like
    s1.bam  Sample 1
    s2.bam  Sample 2
    """
    sample_results = []
    sample_names = []
    if file_path.endswith('.bam'):
        sample_results.append(find_depths_in_bam(file_path))
        sample_names.append('')
    else:
        with open(file_path, 'r') as f:
            samples = [line.split('\t') for line in f.read().split('\n')]
        for sample in samples:
            if sample[0].endswith('.bam'): # Mostly for filtering empty
                sample_results.append(find_depths_in_bam(sample[0]))
                sample_names.append(sample[1])
    plot_depths(sample_results, sample_names)


def gc_depth(file_path):
    """
    Accepts either a bam file or a tab delimited  txt file like
    s1.bam  Sample 1
    s2.bam  Sample 2
    """
    sample_results = []
    sample_names = []
    if file_path.endswith('.bam'):
        sample_results.append(find_depths_in_bam(file_path))
        sample_names.append('')
    else:
        with open(file_path, 'r') as f:
            samples = [line.split('\t') for line in f.read().split('\n')]
        for sample in samples:
            if sample[0].endswith('.bam'): # Mostly for filtering empty
                sample_results.append(find_depths_in_bam(sample[0]))
                sample_names.append(sample[1])
    plot_depths_gc(sample_results, sample_names)
