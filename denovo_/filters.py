from __future__ import division
import helper


class Filters(object):


    def __init__(self, detailed):

        self.detailed = detailed
        self.filter = []


    def apply_filters(self, var_key, variant, config, data):

        csn_key = (variant['gene'], variant['csn'])

        self.filter = []

        try:

            # Check if variant is multiallelic
            self._apply(
                config['REMOVE_MULTI_ALLELE_CALLS'] and var_key[:3] in data['multiallelic_calls'],
                'multi_allele_call'
            )

            # Check if variant is called in either parent
            self._apply(
                var_key in data['mother_var'] or var_key in data['father_var'],
                'called_in_parent'
            )

            # Check if variant is "low" quality (as flagged by postCAVA)
            self._apply(
                variant['quality'] == 'low',
                'low_quality'
            )

            # Check if variant is outside splice site boundary
            self._apply(
                not helper.within_splice_site_boundary(variant['csn'], config['SPLICE_SITE_BOUNDARY']),
                'outside_splice_site_boundary'
            )

            # Check TR in the child
            self._apply(
                variant['TR'] < config['CHILD_MIN_TR'],
                'low_child_tr ({})'.format(variant['TR'])
            )

            # Check TC in the child
            self._apply(
                variant['TC'] < config['CHILD_MIN_TC'],
                'low_child_tc ({})'.format(variant['TC'])
            )

            # Check TR/TC in the child
            self._apply(
                variant['TR'] / variant['TC'] < config['CHILD_MIN_TR_PER_TC'],
                'low_child_tr_per_tc ({})'.format(round(variant['TR'] / variant['TC'], 2))
            )

            # Check control variant frequency
            control_freq = data['control'][csn_key] if csn_key in data['control'] else 0.0
            self._apply(
                control_freq > config['CONTROL_MAX_FREQUENCY'],
                'high_control_frequency ({})'.format(control_freq)
            )

            # Check gnomAD exomes variant frequency
            gnomad_exomes_freq, pop_gnomad_exomes = data['gnomad_exomes_reader'].get_max_frequency(var_key, csn_key)
            self._apply(
                gnomad_exomes_freq > config['GNOMAD_MAX_FREQUENCY'],
                'high_gnomad_exomes_frequency ({})'.format(round(gnomad_exomes_freq, 2))
            )

            # Check gnomAD genomes variant frequency
            gnomad_genomes_freq, pop_gnomad_genomes = data['gnomad_genomes_reader'].get_max_frequency(var_key, csn_key)
            self._apply(
                gnomad_genomes_freq > config['GNOMAD_MAX_FREQUENCY'],
                'high_gnomad_genomes_frequency ({})'.format(round(gnomad_genomes_freq, 2))
            )

            # Count alleles in parents
            parent_alleles = helper.count_parent_alleles(data['mother_bam'], data['father_bam'], var_key)

            # Check TC and TR in the mother
            self._apply(
                parent_alleles['mother_tc'] < config['PARENT_MIN_COVERAGE'],
                'low_mother_tc ({})'.format(parent_alleles['mother_tc'])
            )
            self._apply(
                parent_alleles['mother_tr'] > config['PARENT_MAX_ALT_ALLELE_COUNT'],
                'high_mother_tr ({})'.format(parent_alleles['mother_tr'])
            )

            # Check TC and TR in the father
            self._apply(
                parent_alleles['father_tc'] < config['PARENT_MIN_COVERAGE'],
                'low_father_tc ({})'.format(parent_alleles['father_tc'])
            )
            self._apply(
                parent_alleles['father_tr'] > config['PARENT_MAX_ALT_ALLELE_COUNT'],
                'high_father_tr ({})'.format(parent_alleles['father_tr'])
            )

            return {
                'filter': ','.join(self.filter) if len(self.filter) > 0 else '.',
                'control_freq': control_freq,
                'gnomad_exomes_freq': gnomad_exomes_freq,
                'gnomad_genomes_freq': gnomad_genomes_freq,
                'parent_alleles': parent_alleles,
                'pop_gnomad_exomes': pop_gnomad_exomes,
                'pop_gnomad_genomes': pop_gnomad_genomes
            }

        except ValueError as err:

            return {'filter': str(err)}


    def _apply(self, condition, txt):

        if condition:
            if self.detailed:
                self.filter.append(txt)
            else:
                raise ValueError(txt)















