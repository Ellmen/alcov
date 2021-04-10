# Note: Still in active development

def consensus_from_bam(bam_path):
    import pysam

    samfile = pysam.Samfile(bam_path, "rb")

    seq = ['N' for _ in range(29903)]

    for pileupcolumn in samfile.pileup():
        pos = pileupcolumn.pos
        bases = {}
        for pileupread in pileupcolumn.pileups:
            qpos = pileupread.query_position
            if qpos is None:
                continue
            base = pileupread.alignment.query_sequence[qpos]
            seq[pos] = base
    samfile.close()

    print(''.join(seq))
