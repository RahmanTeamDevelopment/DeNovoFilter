from __future__ import division
from collections import OrderedDict
import sys



def parse_vcf_record(line):

    cols = line.split()

    info = dict(x.split('=', 1) for x in cols[7].split(';') if '=' in x)

    if float(info['TC']) == 0:
        return

    by_alt = dict((key, info[key].split(',')) for key in info)

    chrom = cols[0] if not cols[0].startswith('chr') else cols[0][3:]

    ret = {}
    for i, alt in enumerate(cols[4].split(",")):

        var_key = (chrom, cols[1], cols[3], alt)

        by_transcript = dict((key, by_alt[key][i].split(':')) for key in ['GENE', 'CSN', 'CLASS', 'ALTANN', 'ALTCLASS'])

        if by_alt['TYPE'][i] == 'Substitution':
            qual_flag = 'high' if float(cols[5]) >= 100 else 'low'
        else:
            prop = float(by_alt['TR'][i]) / float(info['TC'])
            qual_flag = 'high' if prop > 0.2 and cols[6] == 'PASS' else 'low'

        ret[var_key] = []
        for j in range(len(by_transcript['GENE'])):

            ret[var_key].append(
                {
                    'quality': qual_flag,
                    'TR': int(by_alt['TR'][i]),
                    'TC': int(info['TC']),
                    'NF': int(by_alt['NF'][i]),
                    'NR': int(by_alt['NR'][i]),
                    'gene': by_transcript['GENE'][j],
                    'csn': by_transcript['CSN'][j],
                    'class_': by_transcript['CLASS'][j],
                    'altann': by_transcript['ALTANN'][j],
                    'altclass': by_transcript['ALTCLASS'][j]
                }
            )

    return ret


def read_vcf_file(fn):

    ret = OrderedDict()
    for line in open(fn):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        parsed_record = parse_vcf_record(line)
        if parsed_record is not None:
            ret.update()

    return ret


def read_config_file(fn):

    sys.stdout.write('\nReading configuration file ... ')
    sys.stdout.flush()

    ret = {
        'GNOMAD_MAX_FREQUENCY': 0.1,
        'CONTROL_MAX_FREQUENCY': 0.1,
        'SPLICE_SITE_BOUNDARY': 10,
        'CHILD_MIN_TR': 3,
        'CHILD_MIN_TC': 15,
        'CHILD_MIN_TR_PER_TC': 0.2,
        'PARENT_MIN_COVERAGE': 6,
        'PARENT_MAX_ALT_ALLELE_COUNT': 1,
        'REMOVE_MULTI_ALLELE_CALLS': 'true',
        'GNOMAD_EXOMES_DATA_FILE': '',
        'GNOMAD_GENOMES_DATA_FILE': '',
        'CONTROL_DATA_FILE': '',
        'MAXENTSCAN_DATA_FILE': '',
        'EXAC_DATA_FILE': ''
    }

    if fn is not None:
        with open(fn) as f:
            for line in f:
                line = line.strip()
                if line == '' or line[0] == '#':
                    continue

                [key, value] = line.split('=')
                key = key.strip().upper()
                value = value.strip()
                if key in ret:
                    ret[key] = value

    for k in ['GNOMAD_MAX_FREQUENCY', 'CONTROL_MAX_FREQUENCY', 'CHILD_MIN_TR_PER_TC']:
        ret[k] = float(ret[k])

    for k in [
        'SPLICE_SITE_BOUNDARY',
        'CHILD_MIN_TR',
        'CHILD_MIN_TC',
        'PARENT_MIN_COVERAGE',
        'PARENT_MAX_ALT_ALLELE_COUNT'
    ]:
        ret[k] = int(ret[k])

    ret['REMOVE_MULTI_ALLELE_CALLS'] = (ret['REMOVE_MULTI_ALLELE_CALLS'].upper() == 'TRUE')

    print ' - Done.'

    return ret



def read_custom_database_file_by_csnkey(fn):

    ret = {}
    num_of_samples = None
    with open(fn) as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue
            if line.startswith('##SAMPLES='):
                num_of_samples = int(line[line.find('=')+1:])
                continue
            if line[0] == '#':
                continue

            cols = line.split('\t')

            key = (cols[6], cols[9])
            ret[key] = 100 * int(cols[20]) / num_of_samples

    return ret



def read_custom_database_file_by_varkey(fn):

    ret = {}
    num_of_samples = None
    with open(fn) as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue
            if line.startswith('##SAMPLES='):
                num_of_samples = int(line[line.find('=')+1:])
                continue
            if line[0] == '#':
                continue

            cols = line.split('\t')

            key = tuple(cols[:4])
            ret[key] = 100 * int(cols[20]) / num_of_samples

    return ret


def read_exac_data_file(fn):

    ret = {}
    with open(fn) as f:
        for line in f:
            line = line.strip()
            if 'N_missense' in line or line == '':
                continue
            cols = line.split('\t')

            ret[cols[0]] = {
                'N_missense': int(cols[1]),
                'Exp_missense': float(cols[2]),
                'Z_missense': float(cols[3]),
                'N_lof': int(cols[4]),
                'Exp_lof': float(cols[5]),
                'pLI': float(cols[6]),
            }
    return ret