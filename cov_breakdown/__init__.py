from .amplicon_coverage import amplicon_coverage
from .analyze import find_mutants
from .cmds import aa, nt
# from .consensus import consensus_from_bam


class CovBreakdown(object):
    """Identify frequencies of concerning mutations from aligned reads"""

    def __str__(self):
        return "Identify frequencies of concerning mutations from aligned reads"

    def aa(self, mut):
        aa(mut)

    def nt(self, mut):
        nt(mut)

    def find_mutants(self, file_path):
        find_mutants(file_path)

    # def consensus(self, file_path):
    #     consensus_from_bam(file_path)

    def amplicon_coverage(self, file_path):
        amplicon_coverage(file_path)
