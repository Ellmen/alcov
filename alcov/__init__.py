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

    def find_mutants(self, samples_path, mutations_path=None):
        find_mutants(samples_path, mutations_path)

    def find_lineages(self, samples_path):
        find_lineages(samples_path)

    # def consensus(self, file_path):
    #     consensus_from_bam(file_path)

    def amplicon_coverage(self, samples_path):
        amplicon_coverage(samples_path)

    def gc_depth(self, samples_path):
        gc_depth(samples_path)
