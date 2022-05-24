import json

from Bio import SeqIO


def process_reference():
    cov = list(SeqIO.parse("sequence.gb", "genbank"))[0]

    genes = {}

    for f in cov.features:
        if f.type == 'gene':
            gene = f.qualifiers['gene'][0]
            start = int(f.location.start)
            end = int(f.location.end)
            genes[gene] = [start, end]

    with open('sars_cov_2.py', 'w') as f:
        f.write('genes = {}\n\nseq = \'{}\''.format(genes, cov.seq))


def get_amplicons():
    cov = list(SeqIO.parse("sequence.gb", "genbank"))[0]
    with open('nCoV-2019.insert.bed', 'r') as f:
        inserts = [l.split('\t') for l in f.read().split('\n')[:-1]]
    with open('SARS-CoV-2.scheme.bed', 'r') as f:
        primers = [l.split('\t') for l in f.readlines()]

    for i in range(len(inserts)):
        insert = inserts[i]
        section = cov.seq[int(insert[1]):int(insert[2])]
        gc = (section.count('G') + section.count('C')) / len(section)
        inserts[i].append(gc)

    with open('artic_amplicons.py', 'w') as f:
        f.write('inserts = {}\nprimers = {}'.format(inserts, primers))


def fix_mut_name(old_mut_name):
    mut_name = old_mut_name.upper()
    if '/' in mut_name:
        # If multiple amino acid dels, only take the first
        mut_name = mut_name[:mut_name.find('/')]
    if not mut_name.startswith('ORF'):
        return mut_name
    orf = 'ORF'
    col_idx = mut_name.find(':')
    mid = mut_name[3:col_idx]
    end = mut_name[col_idx:]
    return orf + mid.lower() + end


def get_mutations():
    # TODO figure out capitalization
    with open('mutations.json', 'r') as f:
        raw_mutations = json.loads(f.read())

    muts = list(set([fix_mut_name(m['mutation']) for m in raw_mutations]))
    lins = list(set([m['pangolin_lineage'].upper() for m in raw_mutations]))
    mut_lins = {mut: {lin: 0 for lin in lins} for mut in muts}
    for raw_m in raw_mutations:
        mut = fix_mut_name(raw_m['mutation'])
        lin = raw_m['pangolin_lineage'].upper()
        prev = raw_m['prevalence']
        mut_lins[mut][lin] = prev

    with open('mutations.py', 'w') as f:
         f.write('mutations = {}'.format(mut_lins))


def get_who_mutations():
    with open('clusters.json', 'r') as f:
        lineages = json.loads(f.read())['clusters']
    vocs = []
    for lin in lineages:
        if 'who_name' in lin and len(lin['who_name']) > 0:
            vocs.append(lin)
    mutations = {}
    who_names = [voc['who_name'][0] for voc in vocs] + ['BA.1', 'BA.2']
    print(who_names)
    who_names = [wn for wn in who_names if wn != 'Omicron']
    print(who_names)
    for voc in vocs:
        who_name = voc['who_name'][0]
        adn = 'alt_display_name'
        if adn in voc and len(voc[adn]) and voc[adn][0] == 'BA.2':
            who_name = 'BA.2'
        # if adn in voc and len(voc[adn]) and voc[adn][0] == 'BA.1':
        if who_name == 'Omicron':
            who_name = 'BA.1'
        muts = voc['mutations']['nonsynonymous']
        muts += voc['mutations']['synonymous']
        for mut in muts:
            if 'gene' in mut:
                if mut['gene'] == 'ORF9b':
                    continue
                if mut['right'] == '-':
                    mut_name = '{}:DEL{}'.format(mut['gene'], mut['pos'])
                else:
                    mut_name = '{}:{}{}{}'.format(
                        mut['gene'],
                        mut['left'],
                        mut['pos'],
                        mut['right']
                    )
            else:
                mut_name = '{}{}{}'.format(
                    mut['left'],
                    mut['pos'],
                    mut['right']
                )
            # TODO: mutations after deletions don't work
            blacklist = [
                'ORF1a:P2287S',
                # Omicron post-deletion
                'ORF1a:L2084I',
                'ORF1a:I3758V',
                'S:Y145D',
                'S:L212I',
                # BA.2 post-deletion
                'S:A27S',
                # Possibly in Omicron
                'ORF1a:DEL3677',
                'S:R346K',
                # Delta post-deletion
                'S:R158G'
            ]
            if mut_name in blacklist:
                continue
            if mut_name not in mutations:
                mutations[mut_name] = {who_name: 0 for who_name in who_names}
            mutations[mut_name][who_name] = 1

    with open('mutations.py', 'w') as f:
         f.write('mutations = {}'.format(mutations))


def get_constellations():
    import json
    from os import getenv, listdir, mkdir, system
    from os.path import isfile, join

    def parse_mut(mut):
        mut = mut.upper()
        if mut[0].isnumeric():
            mut = 'ORF' + mut
        mut = mut.replace('SPIKE', 'S')
        mut = mut.replace('NUC:', '')
        mut = mut.replace('ORF1AB', 'ORF1ab')
        mut = mut.replace('ORF1A', 'ORF1a')
        mut = mut.replace('ORF1B', 'ORF1b')
        mut = mut.replace('ORF3A', 'ORF3a')
        mut = mut.replace('ORF7A', 'ORF7a')
        mut = mut.replace('ORF7b', 'ORF7b')
        return mut

    home = getenv("HOME")
    const_dir = '{}/Code/constellations'.format(home) # Fill in
    data_path = '{}/constellations/definitions'.format(const_dir)

    fns = [f for f in listdir(data_path) if isfile(join(data_path, f)) and f.endswith('.json')]
    fns.sort()
    constellations = {}
    mutations = {}
    for fn in fns:
        with open('{}/{}'.format(data_path, fn), 'r') as f:
            voc = json.loads(f.read())
        label = voc['label']
        constellations[label] = [parse_mut(mut) for mut in voc['sites']]
    with open('constellations.py', 'w') as f:
         f.write('constellations = {}'.format(constellations))
    vocs = list(constellations.keys())
    # Remove non-specific lineages
    constellations['Omicron (XE-like)'] += constellations['XE-parent1']
    constellations['Omicron (XE-like)'] += constellations['XE-parent2']
    for voc in constellations:
        if 'Omicron' in voc:
            constellations[voc] += constellations['Omicron (Unassigned)']
    vocs.remove('Omicron (Unassigned)')
    vocs.remove('XE-parent1')
    vocs.remove('XE-parent2')
    for voc in vocs:
        for mut_name in constellations[voc]:
            # if mut_name in blacklist:
            #     continue
            if '+' in mut_name: # TODO: support
                continue
            if mut_name.startswith('NSP'): # TODO: support
                continue
            if mut_name not in mutations:
                mutations[mut_name] = {voc: 0 for voc in vocs}
            mutations[mut_name][voc] = 1
    with open('mutations.py', 'w') as f:
         f.write('mutations = {}'.format(mutations))


def get_challenge_mutations():
    with open('challenge_mutations.json', 'r') as f:
        lineages = json.loads(f.read())
    mutations = {}
    # lins = list(lineages.keys())
    lins = {'EPI_ISL_7190366': 'BA.2', 'EPI_ISL_7877191': 'BA.1', 'EPI_ISL_2793160': 'Delta', 'EPI_ISL_10389336': 'Recombinant'}
    for lin in lins:
    # for lin in lineages:
        for mut_name, loc in lineages[lin]:
            # if 'ins' in mut_name: # TODO: support
            #     continue
            mut_name = '{}@{}'.format(mut_name, loc)
            if mut_name not in mutations:
                mutations[mut_name] = {lins[lin]: 0 for lin in lins}
            mutations[mut_name][lins[lin]] = 1
            # if mut_name not in mutations:
            #     mutations[mut_name] = {lin: 0 for lin in lineages}
            # mutations[mut_name][lin] = 1
    with open('mutations.py', 'w') as f:
         f.write('mutations = {}'.format(mutations))


if __name__ == '__main__':
    # process_reference()
    # get_amplicons()
    # get_mutations()
    # get_who_mutations()
    # get_constellations()
    get_challenge_mutations()
