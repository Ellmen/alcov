import pysam


# TODO: Abstract this process by parsing associated mutations (also clean up)
b117_mutations = [
    'A28271-',
    'A23063T'
]


def find_mutants(bam_path):
    samfile = pysam.Samfile(bam_path, "rb")

    dels = 0
    not_dels = 0
    snps = 0
    not_snps = 0

    for pileupcolumn in samfile.pileup():
        pos = pileupcolumn.pos
        if pos == 28270:
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del:
                    dels += 1
                else:
                    not_dels += 1
        elif pos == 23062:
            for pileupread in pileupcolumn.pileups:
                qpos = pileupread.query_position
                if qpos is None:
                    continue
                base = pileupread.alignment.query_sequence[qpos]
                if base == 'T':
                    snps += 1
                else:
                    not_snps += 1
    samfile.close()

    print('A28271-:')
    total_dels = dels + not_dels
    if total_dels == 0:
        print('No coverage of A28271-')
    else:
        print('{} are deletions, {} are wildtype ({:.2f}% of {} total)'.format(dels, not_dels, dels/total_dels*100, total_dels))

    print('A23063T:')
    total_snps = snps + not_snps
    if total_snps == 0:
        print('No coverage of A23063T')
    else:
        print('{} are T, {} are wildtype ({:.2f}% of {} total)'.format(snps, not_snps, snps/total_snps*100, total_snps))
