from __future__ import division
import alleles
import sys


def is_candidate(var_key, data, config, multiallelic_calls, mother_var_data, father_var_data, freqs, mother_bam, father_bam):

    ret = True
    reason = []

    # Check if variant is multiallelic
    if config['REMOVE_MULTI_ALLELE_CALLS'] and var_key[:3] in multiallelic_calls:
        reason.append('multiallelic_call')
        ret = False

    # Check if variant is called in either parent
    if var_key in mother_var_data:
        reason.append('found_in_mother')
        ret = False
    if var_key in father_var_data:
        reason.append('found_in_father')
        ret = False

    # Check if variant is "low" quality (as flagged by postCAVA.py)
    if data['quality'] == 'low':
        reason.append('low_quality')
        ret = False

    # Check if variant is outside splice site boundary
    if not within_splice_site_boundary(data['csn'], config['SPLICE_SITE_BOUNDARY']):
        reason.append('outside_splice_site_boundary')
        ret = False

    # Check gnomAD variant frequency
    if freqs['gnomad'] > config['GNOMAD_MAX_FREQUENCY']:
        reason.append('high_gnomad_frequency')
        ret = False

    # Check control variant frequency
    if freqs['control'] > config['CONTROL_MAX_FREQUENCY']:
        reason.append('high_control_frequency')
        ret = False

    # Check TR in the child
    if data['TR'] < config['CHILD_MIN_TR']:
        reason.append('low_child_tr')
        ret = False

    # Check TC In the child
    if data['TC'] < config['CHILD_MIN_TC']:
        reason.append('low_child_tc')
        ret = False

    # Check TR/TC in the child
    if data['TR'] / data['TC'] < config['CHILD_MIN_TR_PER_TC']:
        reason.append('low_child_tr_per_tc')
        ret = False

    # Check TC and TR in the mother
    mother_tc, mother_tr = alleles.count(mother_bam, var_key)
    if mother_tc < config['PARENT_MIN_COVERAGE']:
        reason.append('low_mother_tc')
        ret = False
    if mother_tr >= config['PARENT_MAX_ALT_ALLELE_COUNT']:
        reason.append('high_mother_tr')
        ret = False

    # Check TC and TR in the father
    father_tc, father_tr = alleles.count(father_bam, var_key)
    if father_tc < config['PARENT_MIN_COVERAGE']:
        reason.append('low_father_tc')
        ret = False
    if father_tr >= config['PARENT_MAX_ALT_ALLELE_COUNT']:
        reason.append('high_father_tr')
        ret = False

    return ret, reason


def find_multiallelic_calls(var_data):

    counts = {}
    for var_key in var_data.keys():
        k = var_key[:3]
        if k not in counts:
            counts[k] = 1
        else:
            counts[k] += 1
    return [x for x in counts.keys() if counts[k] > 1]


def within_splice_site_boundary(csn, cutoff):

    csn_c = csn.split('_p.', 1)[0] if '_p.' in csn else csn
    cutidx = [csn_c.find(x) for x in ['A>', 'C>', 'G>', 'T>', 'del', 'ins', 'dup'] if x in csn_c]
    coords = csn_c[csn_c.find("c.") + 2:min(cutidx)]

    if "_" in coords:
        [coord_part1, coord_part2] = coords.split('_', 1)
    else:
        coord_part1, coord_part2 = coords, None

    # Split first coordinate to exon and intron part
    coordinate1, intron_coordinate1 = _split_to_exon_and_intron_coordinates(coord_part1)

    if intron_coordinate1 is not None:
        if abs(int(intron_coordinate1)) > cutoff:
            return False

    # Split second coordinate to exon and intron part (if any)
    if coord_part2 is not None:
        coordinate2, intron_coordinate2 = _split_to_exon_and_intron_coordinates(coord_part2)
    else:
        coordinate2 = intron_coordinate2 = None

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


def output_header(outfile):

    header = [
        'CHROM',
        'POS',
        'REF',
        'ALT',
        'GENE'
        'CSN',
        'CLASS',
        'ALTANN',
        'ALTCLASS',
        'GNOMAD_freq',
        'CONTROL_freq',
        'TR',
        'TC'
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
        'MaxEntScan_MAX3',
        'Info'
    ]
    outfile.write('\t'.join(header)+'\n')


def output(outfile, var_key, data, freqs, maxentscan_scores, reason):

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
        freqs['gnomad'],
        freqs['control'],
        data['TR'],
        data['TC']
    ]

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
    if maxentscan_scores is not None:
        for c in maxentscan_columns:
            record.append(maxentscan_scores[c])
    else:
        record += ['.'] * len(maxentscan_columns)

    record += [reason]

    record = map(str, record)
    outfile.write('\t'.join(record)+'\n')


def welcome(version):

    print '\n{} DeNovoFilter {} {}'.format('='*3, version, '='*80)


def goodbye(counter_denovo, counter_filtered, output_prefix):

    print '\nNumber of de novo candidates: {}  ({}_denovo_candidates.txt)'.format(counter_denovo, output_prefix)
    print 'Number of variants filtered out: {}  ({}_filtered.txt)'.format(counter_filtered, output_prefix)
    print '\n{}\n'.format('=' * 103)


def print_info(options, config):

    print '\nInput files:'
    print '-' * 60
    print 'Variants data file of the child: {}'.format(options.child_var)
    print 'Variants data file of the mother: {}'.format(options.mother_var)
    print 'Variants data file of the father: {}'.format(options.father_var)
    print '-' * 60
    print 'BAM file of the child: {}'.format(options.child_bam)
    print 'BAM file of the mother: {}'.format(options.mother_bam)
    print 'BAM file of the father: {}'.format(options.father_bam)
    print '-' * 60
    print 'gnomAD data file: {}'.format(config['GNOMAD_DATA_FILE'])
    print 'Control data file: {}'.format(config['CONTROL_DATA_FILE'])
    print 'MaxEntScan data file: {}'.format(config['MAXENTSCAN_DATA_FILE'])


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