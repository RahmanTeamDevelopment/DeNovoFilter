from __future__ import division
from collections import OrderedDict
import sys


def read_variant_file(fn):

    ret = OrderedDict()

    idx = {}

    with open(fn) as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue

            if line[0] == '#':
                header = line[1:].split('\t')
                for t in ['QUALFLAG', 'TR', 'TC', 'GENE', 'CSN', 'CLASS', 'ALTANN', 'ALTCLASS']:
                    idx[t] = header.index(t)

            cols = line.split('\t')
            key = tuple(cols[:4])
            if key not in ret:
                ret[key] = []
            ret[key].append(
                {
                    'quality': cols[idx['QUALFLAG']],
                    'TR': int(cols[idx['TR']]),
                    'TC': int(cols[idx['TC']]),
                    'gene': cols[idx['GENE']],
                    'csn': cols[idx['CSN']],
                    'class_': cols[idx['CLASS']],
                    'altann': cols[idx['ALTANN']],
                    'altclass': cols[idx['ALTCLASS']]
                }
            )

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