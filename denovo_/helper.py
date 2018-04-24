from __future__ import division
import sys
import alleles


def check_fast(var_key, data, config, multiallelic_calls, mother_var_data, father_var_data, gnomad_exomes_file, gnomad_genomes_file, control_data, mother_bam, father_bam):

    csn_key = (data['gene'], data['csn'])

    # Check if variant is multiallelic
    if config['REMOVE_MULTI_ALLELE_CALLS'] and var_key[:3] in multiallelic_calls:
        return {'filter': 'multi_allele_call'}

    # Check if variant is called in either parent
    if var_key in mother_var_data or var_key in father_var_data:
        return {'filter': 'called_in_parent'}

    # Check if variant is "low" quality (as flagged by postCAVA.py)
    if data['quality'] == 'low':
        return {'filter': 'low_quality'}

    # Check if variant is outside splice site boundary
    if not within_splice_site_boundary(data['csn'], config['SPLICE_SITE_BOUNDARY']):
        return {'filter': 'outside_splice_site_boundary'}

    # Check TR in the child
    if data['TR'] < config['CHILD_MIN_TR']:
        return {'filter': 'low_child_tr ({})'.format(data['TR'])}

    # Check TC In the child
    if data['TC'] < config['CHILD_MIN_TC']:
        return {'filter': 'low_child_tc ({})'.format(data['TC'])}

    # Check TR/TC in the child
    if data['TR'] / data['TC'] < config['CHILD_MIN_TR_PER_TC']:
        return {'filter': 'low_child_tr_per_tc ({})'.format(data['TR'] / data['TC'])}

    # Calculate control frequency
    control_freq = control_data[csn_key] if csn_key in control_data else 0.0

    # Check control variant frequency
    if control_freq > config['CONTROL_MAX_FREQUENCY']:
        return {'filter': 'high_control_frequency ({})'.format(control_freq)}

    # Calculate gnomAD frequencies
    gnomad_exomes_freq, pop_gnomad_exomes = read_gnomad_data(gnomad_exomes_file, var_key, csn_key)
    gnomad_genomes_freq, pop_gnomad_genomes = read_gnomad_data(gnomad_genomes_file, var_key, csn_key)

    # Check gnomAD exomes variant frequency
    if gnomad_exomes_freq != 'NA':
        if gnomad_exomes_freq > config['GNOMAD_MAX_FREQUENCY']:
            return {'filter': 'high_gnomad_exomes_frequency ({})'.format(gnomad_exomes_freq)}

    # Check gnomAD genomes variant frequency
    if gnomad_genomes_freq != 'NA':
        if gnomad_genomes_freq > config['GNOMAD_MAX_FREQUENCY']:
            return {'filter': 'high_gnomad_genomes_frequency ({})'.format(gnomad_genomes_freq)}

    # Count alleles in parents
    parent_alleles = count_parent_alleles(mother_bam, father_bam, var_key)

    # Check TC and TR in the mother
    if parent_alleles['mother_tc'] < config['PARENT_MIN_COVERAGE']:
        return {'filter': 'low_mother_tc ({})'.format(parent_alleles['mother_tc'])}
    if parent_alleles['mother_tr'] > config['PARENT_MAX_ALT_ALLELE_COUNT']:
        return {'filter': 'high_mother_tr ({})'.format(parent_alleles['mother_tr'])}

    # Check TC and TR in the father
    if parent_alleles['father_tc'] < config['PARENT_MIN_COVERAGE']:
        return {'filter': 'low_father_tc ({})'.format(parent_alleles['father_tc'])}
    if parent_alleles['father_tr'] > config['PARENT_MAX_ALT_ALLELE_COUNT']:
        return {'filter': 'high_father_tr ({})'.format(parent_alleles['father_tr'])}

    return {
        'filter': None,
        'control_freq': control_freq,
        'gnomad_exomes_freq': gnomad_exomes_freq,
        'gnomad_genomes_freq': gnomad_genomes_freq,
        'parent_alleles': parent_alleles,
        'pop_gnomad_exomes': pop_gnomad_exomes,
        'pop_gnomad_genomes': pop_gnomad_genomes
    }


def check_slow(var_key, data, config, multiallelic_calls, mother_var_data, father_var_data, gnomad_exomes_file, gnomad_genomes_file, control_data, mother_bam, father_bam):

    csn_key = (data['gene'], data['csn'])

    filter = []

    # Check if variant is multiallelic
    if config['REMOVE_MULTI_ALLELE_CALLS'] and var_key[:3] in multiallelic_calls:
        filter.append('multi_allele_call')

    # Check if variant is called in either parent
    if var_key in mother_var_data or var_key in father_var_data:
        filter.append('called_in_parent')

    # Check if variant is "low" quality (as flagged by postCAVA.py)
    if data['quality'] == 'low':
        filter.append('low_quality')

    # Check if variant is outside splice site boundary
    if not within_splice_site_boundary(data['csn'], config['SPLICE_SITE_BOUNDARY']):
        filter.append('outside_splice_site_boundary')

    # Check TR in the child
    if data['TR'] < config['CHILD_MIN_TR']:
        filter.append('low_child_tr ({})'.format(data['TR']))

    # Check TC In the child
    if data['TC'] < config['CHILD_MIN_TC']:
        filter.append('low_child_tc ({})'.format(data['TC']))

    # Check TR/TC in the child
    if data['TR'] / data['TC'] < config['CHILD_MIN_TR_PER_TC']:
        filter.append('low_child_tr_per_tc ({})'.format(data['TR'] / data['TC']))

    # Calculate control frequency
    control_freq = control_data[csn_key] if csn_key in control_data else 0.0

    # Check control variant frequency
    if control_freq > config['CONTROL_MAX_FREQUENCY']:
        filter.append('high_control_frequency ({})'.format(control_freq))

    # Calculate gnomAD frequencies
    gnomad_exomes_freq, pop_gnomad_exomes = read_gnomad_data(gnomad_exomes_file, var_key, csn_key)
    gnomad_genomes_freq, pop_gnomad_genomes = read_gnomad_data(gnomad_genomes_file, var_key, csn_key)

    # Check gnomAD exomes variant frequency
    if gnomad_exomes_freq != 'NA':
        if gnomad_exomes_freq > config['GNOMAD_MAX_FREQUENCY']:
            filter.append('high_gnomad_exomes_frequency ({})'.format(gnomad_exomes_freq))

    # Check gnomAD genomes variant frequency
    if gnomad_genomes_freq != 'NA':
        if gnomad_genomes_freq > config['GNOMAD_MAX_FREQUENCY']:
            filter.append('high_gnomad_genomes_frequency ({})'.format(gnomad_genomes_freq))

    # Count alleles in parents
    parent_alleles = count_parent_alleles(mother_bam, father_bam, var_key)

    # Check TC and TR in the mother
    if parent_alleles['mother_tc'] < config['PARENT_MIN_COVERAGE']:
        filter.append('low_mother_tc ({})'.format(parent_alleles['mother_tc']))
    if parent_alleles['mother_tr'] > config['PARENT_MAX_ALT_ALLELE_COUNT']:
        filter.append('high_mother_tr ({})'.format(parent_alleles['mother_tr']))

    # Check TC and TR in the father
    if parent_alleles['father_tc'] < config['PARENT_MIN_COVERAGE']:
        filter.append('low_father_tc ({})'.format(parent_alleles['father_tc']))
    if parent_alleles['father_tr'] > config['PARENT_MAX_ALT_ALLELE_COUNT']:
        filter.append('high_father_tr ({})'.format(parent_alleles['father_tr']))

    return {
        'filter': ','.join(filter) if len(filter) > 0 else None,
        'control_freq': control_freq,
        'gnomad_exomes_freq': gnomad_exomes_freq,
        'gnomad_genomes_freq': gnomad_genomes_freq,
        'parent_alleles': parent_alleles,
        'pop_gnomad_exomes': pop_gnomad_exomes,
        'pop_gnomad_genomes': pop_gnomad_genomes
    }


def count_parent_alleles(mother_bam, father_bam, var_key):

    ret = {}
    ret['mother_tc'], ret['mother_tr'] = alleles.count(mother_bam, var_key)
    ret['father_tc'], ret['father_tr'] = alleles.count(father_bam, var_key)
    return ret


def frequencies(gnomad_file, control_data, var_key, csn_key):

    ret = {}
    ret['gnomad'] = read_gnomad_data(gnomad_file, var_key, csn_key)
    if csn_key in control_data:
        ret['control'] = control_data[csn_key]
    else:
        ret['control'] = 0.0
    return ret


def find_multiallelic_calls(var_data):

    counts = {}
    for var_key in var_data.keys():
        k = var_key[:3]
        if k not in counts:
            counts[k] = 1
        else:
            counts[k] += 1
    return [x for x in counts.keys() if counts[x] > 1]


def within_splice_site_boundary(csn, cutoff):

    csn_c = csn.split('_p.', 1)[0] if '_p.' in csn else csn
    cutidx = [csn_c.find(x) for x in ['A>', 'C>', 'G>', 'T>', 'del', 'ins', 'dup'] if x in csn_c]
    coords = csn_c[csn_c.find("c.") + 2:min(cutidx)]

    if "_" in coords:
        [coord_part1, coord_part2] = coords.split('_', 1)
    else:
        coord_part1, coord_part2 = coords, None

    # Split first coordinate to exon and intron part
    _, intron_coordinate1 = _split_to_exon_and_intron_coordinates(coord_part1)

    if intron_coordinate1 is not None:
        if abs(int(intron_coordinate1)) > cutoff:
            return False

    # Split second coordinate to exon and intron part (if any)
    if coord_part2 is not None:
        _, intron_coordinate2 = _split_to_exon_and_intron_coordinates(coord_part2)
    else:
        intron_coordinate2 = None

    if intron_coordinate2 is not None:
        if abs(int(intron_coordinate2)) > cutoff:
            return False

    return True


def _split_to_exon_and_intron_coordinates(coord_part):

    exonpart, intronpart = coord_part, None
    for i in range(1, len(coord_part)):
        if coord_part[i] in ['-', '+']:
            exonpart = coord_part[:i]
            intronpart = coord_part[i:]
            break

    return exonpart, intronpart


def output_header(outfile, maxentscan_columns, exac_columns):

    header = [
        'CHROM',
        'POS',
        'REF',
        'ALT',
        'GENE',
        'CSN',
        'CLASS',
        'ALTANN',
        'ALTCLASS',
        'gnomAD_exomes_frequency',
        'Popuplation_gnomAD_exomes',
        'gnomAD_genomes_frequency',
        'Popuplation_gnomAD_genomes',
        'Control_frequency',
        'Child_TR',
        'Child_TC',
        'Mother_TR',
        'Mother_TC',
        'Father_TR',
        'Father_TC'
    ]
    if maxentscan_columns:
        header += [
            'MaxEntScan_RefKnown5',
            'MaxEntScan_RefKnown3',
            'MaxEntScan_AltKnown',
            'MaxEntScan_AltHighest5',
            'MaxEntScan_RefHighest5',
            'MaxEntScan_AltHighest3',
            'MaxEntScan_RefHighest3',
            'MaxEntScan_Boundary5',
            'MaxEntScan_Boundary3',
            'MaxEntScan_PI5',
            'MaxEntScan_PI3',
            'MaxEntScan_RefKnown',
            'MaxEntScan_SpliceSiteScore',
            'MaxEntScan_SpliceSiteType',
            'MaxEntScan_PercentReduction',
            'MaxEntScan_MAX5',
            'MaxEntScan_MAX3'
        ]

    if exac_columns:
        header += [
            'ExAC_N_missense',
            'ExAC_Exp_missense',
            'ExAC_Z_missense',
            'ExAC_N_lof',
            'ExAC_Exp_lof',
            'ExAC_pLI',
        ]

    outfile.write('\t'.join(header)+'\n')


def output_simple(out, var_key, data, result):

    out.write(
        '\t'.join(
            [
                var_key[0],
                var_key[1],
                var_key[2],
                var_key[3],
                data['gene'],
                data['csn'],
                result['filter']
            ]
        ) + '\n'
    )


def output(outfile, var_key, data, res, maxentscan_scores, exac_values):

    (chrom, pos, ref, alt) = var_key

    record = [
        chrom,
        pos,
        ref,
        alt,
        data['gene'],
        data['csn'],
        data['class_'],
        data['altann'],
        data['altclass'],
        res['gnomad_exomes_freq'],
        res['pop_gnomad_exomes'],
        res['gnomad_genomes_freq'],
        res['pop_gnomad_genomes'],
        res['control_freq'],
        data['TR'],
        data['TC'],
        res['parent_alleles']['mother_tr'],
        res['parent_alleles']['mother_tc'],
        res['parent_alleles']['father_tr'],
        res['parent_alleles']['father_tc']
    ]

    if maxentscan_scores is not None:
        maxentscan_columns = [
            'RefKnown5',
            'RefKnown3',
            'AltKnown',
            'AltHighest5',
            'RefHighest5',
            'AltHighest3',
            'RefHighest3',
            'Boundary5',
            'Boundary3',
            'PI5',
            'PI3',
            'RefKnown',
            'SpliceSiteScore',
            'SpliceSiteType',
            'PercentReduction',
            'MAX5',
            'MAX3'
        ]
        if maxentscan_scores != {}:
            for c in maxentscan_columns:
                record.append(maxentscan_scores[c])
        else:
            record += ['.'] * len(maxentscan_columns)

    if exac_values is not None:
        exac_columns = [
            'N_missense',
            'Exp_missense',
            'Z_missense',
            'N_lof',
            'Exp_lof',
            'pLI',
        ]
        if exac_values != {}:
            for c in exac_columns:
                record.append(exac_values[c])
        else:
            record += ['.'] * len(exac_columns)

    record = map(str, record)
    outfile.write('\t'.join(record)+'\n')


def welcome(version):

    print '\n{} DeNovoFilter {} {}'.format('='*3, version, '='*80)


def goodbye(counter_denovo, counter_filtered, output_prefix):

    print '\nNumber of de novo candidates: {}  ({}_denovo_candidates.txt)'.format(counter_denovo, output_prefix)
    print 'Number of variants filtered out: {}  ({}_filtered.txt)'.format(counter_filtered, output_prefix)
    print '\n{}\n'.format('=' * 103)


def init_progress():

    print ''
    sys.stdout.write('\rProcessing variants ... 0.0%')
    sys.stdout.flush()


def print_progress(counter, total):

    x = round(100 * counter / total, 1)
    x = min(x, 100.0)
    sys.stdout.write('\rProcessing variants ... {}%'.format(x))
    sys.stdout.flush()


def finalize_progress():

    sys.stdout.write('\rProcessing variants ... 100.0%')
    sys.stdout.flush()
    print ' - Done.'


def read_gnomad_data(tabix_file, var_key, csn_key):

    delta = 100

    chrom = var_key[0]
    pos = int(var_key[1])
    gene = csn_key[0]
    csn = csn_key[1]

    if chrom not in tabix_file.contigs:
        return 'NA'

    for line in tabix_file.fetch(chrom, pos - delta, pos + delta):
        line = line.strip()
        cols = line.split('\t')

        if cols[6] != 'PASS':
            continue

        n_alts = len(cols[4].split(','))

        flags = {}
        for x in cols[7].split(';'):
            if '=' not in x:
                continue
            k = x[:x.find('=')]
            v = x[x.find('=') + 1:]
            flags[k] = v


        for i in range(n_alts):

            if not (flags['GENE'].split(',')[i] == gene and flags['CSN'].split(',')[i] == csn):
                    continue

            popfreqs = []
            pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
            for pop in pops:

                key = 'GC_' + pop

                if chrom in ['X', 'chrX']:
                    popfreqs.append(float(frequencies_combined(flags[key + '_Male'], flags[key + '_Female'], n_alts)[i]))

                elif chrom in ['Y', 'chrY']:
                    popfreqs.append(float(flags['AF_' + pop].split(',')[i]))

                else:
                    popfreqs.append(float(frequencies(flags[key], n_alts)[i]))

            if len(popfreqs) > 0:
                return max(popfreqs), pops[popfreqs.index(max(popfreqs))]

    return 0.0


def frequencies(gc, n_alts):

    ret = []

    GCv = map(int, gc.split(','))

    if sum(GCv) == 0:
        return [0.0] * n_alts

    genotypes = []
    for a1 in range(n_alts+1):
        for a2 in range(a1+1):
            genotypes.append((a1,a2))

    for i in range(n_alts):
        N = 0
        for j in range(len(genotypes)):
            if i+1 in genotypes[j]:
                N += GCv[j]
        ret.append(100 * N / sum(GCv))

    return ret


def frequencies_combined(gc_male, gc_female, n_alts):

    ret = []

    GCv_male = map(int, gc_male.split(','))
    GCv_female = map(int, gc_female.split(','))

    total = sum(GCv_male) + sum(GCv_female)

    if total == 0:
        return [0.0] * n_alts

    genotypes = []
    for a1 in range(n_alts+1):
        for a2 in range(a1+1):
            genotypes.append((a1,a2))

    for i in range(n_alts):
        N_male = 0
        N_female = 0
        for j in range(len(genotypes)):
            if i+1 in genotypes[j]:
                N_male += GCv_male[j]
                N_female += GCv_female[j]
        ret.append(100 * (N_male + N_female) / total)

    return ret

