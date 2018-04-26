from __future__ import division
import pysam


class GnomadDBReader(object):


    def __init__(self, format_fn, exomes=True):

        self.exomes = exomes

        self.tabix_files = {}
        if exomes:
            self.tabix_files['_'] = pysam.Tabixfile(format_fn)
        else:
            for chrom in map(str, range(1,23)) + ['X']:
                self.tabix_files[chrom] = pysam.Tabixfile(format_fn.format(chrom))


    def get_max_frequency(self, var_key, csn_key, variant_frequency=True):

        chrom = var_key[0]
        pos = int(var_key[1])
        ref = var_key[2]
        alt = var_key[3]
        gene = csn_key[0]
        csn = csn_key[1]

        for (n_alts, flags) in self._read_variants_in_vicinity(chrom, pos):

            for i in range(n_alts):

                if csn == '.':
                    if not (flags['pos'] == pos and flags['ref'] == ref and flags['alts'][i] == alt):
                        continue
                else:
                    if not (flags['GENE'].split(',')[i] == gene and flags['CSN'].split(',')[i] == csn):
                        continue

                pops = self._extract_pops(flags)
                if variant_frequency:
                    popfreqs = [self._variant_frequency(pop, chrom, flags, n_alts, i) for pop in pops]
                else:
                    popfreqs = [float(flags['AF_' + pop].split(',')[i]) for pop in pops]
                maxfreq = max(popfreqs)

                return maxfreq, pops[popfreqs.index(maxfreq)]

        return 0.0, '.'


    def _read_variants_in_vicinity(self, chrom, pos, delta=100):

        ret = []

        tabix_file = self.tabix_files['_'] if self.exomes else self.tabix_files[chrom]

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

            flags['pos'] = int(cols[1])
            flags['ref'] = cols[3]
            flags['alts'] = cols[4].split(',')

            ret.append((n_alts, flags))

        return ret


    def _variant_frequency(self, pop, chrom, flags, n_alts, alt_idx):

        key = 'GC_' + pop

        if chrom in ['X', 'chrX']:
            return self._gc_mf_to_variant_freqs(flags[key + '_Male'], flags[key + '_Female'], n_alts, alt_idx)

        elif chrom in ['Y', 'chrY']:
            return float(flags['AF_' + pop].split(',')[alt_idx])

        else:
            return self._gc_to_variant_freqs(flags[key], n_alts, alt_idx)


    def _gc_to_variant_freqs(self, gc, n_alts, alt_idx):

        ret = []

        gc_list = map(int, gc.split(','))

        if sum(gc_list) == 0:
            return 0.0

        genotypes = self._generate_genotypes(n_alts)

        for i in range(n_alts):
            N = 0
            for j in range(len(genotypes)):
                if i + 1 in genotypes[j]:
                    N += gc_list[j]
            ret.append(100 * N / sum(gc_list))

        return ret[alt_idx]


    def _gc_mf_to_variant_freqs(self, gc_male, gc_female, n_alts, alt_idx):

        ret = []

        gc_list_male = map(int, gc_male.split(','))
        gc_list_female = map(int, gc_female.split(','))

        total = sum(gc_list_male) + sum(gc_list_female)

        if total == 0:
            return 0.0

        genotypes = self._generate_genotypes(n_alts)

        for i in range(n_alts):
            N_male = 0
            N_female = 0
            for j in range(len(genotypes)):
                if i + 1 in genotypes[j]:
                    N_male += gc_list_male[j]
                    N_female += gc_list_female[j]
            ret.append(100 * (N_male + N_female) / total)

        return ret[alt_idx]


    def _generate_genotypes(self, n_alts):

        ret = []
        for a1 in range(n_alts + 1):
            for a2 in range(a1 + 1):
                ret.append((a1, a2))
        return ret


    def _extract_pops(self, flags):

        return [p for p in ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS'] if 'AF_' + p in flags]