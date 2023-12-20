from .amplicon_coverage import amplicon_coverage, gc_depth
from .analyze import find_mutants
from .cmds import aa, nt
from .lineages import find_lineages
# from .consensus import consensus_from_bam


class CovBreakdown(object):
    """Identify frequencies of concerning mutations from aligned reads"""

    def __str__(self):
        return "Identify frequencies of concerning mutations from aligned reads"

    def aa(self, mut):
        aa(mut)

    def nt(self, mut):
        nt(mut)

    def find_mutants(self, samples_path, mutations_path=None, min_depth=40, save_img=False, csv=False):
        find_mutants(samples_path, mutations_path=mutations_path, min_depth=min_depth, save_img=save_img, csv=csv)

    def find_lineages(self, samples_path, lineages_path=None, ts=False, csv=False, min_depth=40, show_stacked=False, unique=F$
        find_lineages(samples_path, lineages_path=lineages_path, ts=ts, csv=csv, min_depth=min_depth, show_stacked=show_stack$

    def amplicon_coverage(self, samples_path):
        amplicon_coverage(samples_path)

    def gc_depth(self, samples_path):
        gc_depth(samples_path)

