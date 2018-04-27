from __future__ import division
import sys
import alleles
import parsers
import pysam
import gnomad
import maxentscan


def read_data(options, config):

    sys.stdout.write('\nReading input data ... ')
    sys.stdout.flush()

    ret = {}

    # Read variant data of the three individuals
    ret['child_var'] = parsers.read_variant_file(options.child_var)
    ret['mother_var'] = parsers.read_variant_file(options.mother_var)
    ret['father_var'] = parsers.read_variant_file(options.father_var)

    # Create list of multiallelic calls in the child
    ret['multiallelic_calls'] = find_multiallelic_calls(ret['child_var'])

    # Connect to BAM files of the mother and the father
    ret['mother_bam'] = pysam.AlignmentFile(options.mother_bam, "rb")
    ret['father_bam'] = pysam.AlignmentFile(options.father_bam, "rb")

    # Create GnomadDBReader objects for both gnomAD exomes and genomes
    ret['gnomad_exomes_reader'] = gnomad.GnomadDBReader(config['GNOMAD_EXOMES_DATA_FILE'])
    ret['gnomad_genomes_reader'] = gnomad.GnomadDBReader(config['GNOMAD_GENOMES_DATA_FILE'], exomes=False)

    # Read control data
    ret['control_csnkey'] = parsers.read_custom_database_file_by_csnkey(config['CONTROL_DATA_FILE'])
    ret['control_varkey'] = parsers.read_custom_database_file_by_varkey(config['CONTROL_DATA_FILE'])

    # Read MaxEntScan data
    fn = config['MAXENTSCAN_DATA_FILE']
    ret['maxentscan'] = maxentscan.MaxEntScanData(fn) if fn != '' else None

    # Read ExAC data
    fn = config['EXAC_DATA_FILE']
    ret['exac'] = parsers.read_exac_data_file(fn) if fn != '' else None

    print ' - Done.'

    return ret


def count_parent_alleles(mother_bam, father_bam, var_key):

    ret = {}
    ret['mother_tc'], ret['mother_tr'] = alleles.count(mother_bam, var_key)
    ret['father_tc'], ret['father_tr'] = alleles.count(father_bam, var_key)
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


def output_header_simplified(outfile):

    header = [
        'CHROM',
        'POS',
        'REF',
        'ALT',
        'GENE',
        'CSN',
        'Filter'
    ]
    outfile.write('\t'.join(header)+'\n')


def output_header(outfile, maxentscan_columns, exac_columns, output_filter_column):

    header = [
        'CHROM',
        'POS',
        'REF',
        'ALT',
        'GENE',
        'CSN',
        'CLASS',
        'ALTANN',
        'ALTCLASS'
    ]

    if output_filter_column:
        header.append('Filter')

    header += [
        'gnomAD_exomes_frequency',
        'Population_gnomAD_exomes',
        'gnomAD_genomes_frequency',
        'Population_gnomAD_genomes',
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


def output_simplified(out, var_key, variant, result):

    out.write('\t'.join(list(var_key[:4])+[variant['gene'], variant['csn'], result['filter']]) + '\n')


def output(outfile, var_key, data, result, maxentscan_scores, exac_values, output_filter_column):

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
        data['altclass']
    ]

    if output_filter_column:
        record.append(result['filter'])

    record += [
        result['gnomad_exomes_freq'],
        result['pop_gnomad_exomes'],
        result['gnomad_genomes_freq'],
        result['pop_gnomad_genomes'],
        result['control_freq'],
        data['TR'],
        data['TC'],
        result['parent_alleles']['mother_tr'],
        result['parent_alleles']['mother_tc'],
        result['parent_alleles']['father_tr'],
        result['parent_alleles']['father_tc']
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

    record = [str(round(x, 2)) if type(x)==float else str(x) for x in record]
    outfile.write('\t'.join(record)+'\n')


def welcome(version):

    print '\n{} DeNovoFilter {} {}'.format('='*3, version, '='*80)


def goodbye(counter_denovo, counter_filtered, output_prefix):

    print '\nNumber of de novo candidates: {}  ({}_denovo_candidates.txt)'.format(counter_denovo, output_prefix)
    print 'Number of variants filtered out: {}  ({}_filtered_out.txt)'.format(counter_filtered, output_prefix)
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

