from collections import OrderedDict


class MaxEntScanData(object):


    def __init__(self, fn):

        self.data = self._parse_data_file(fn)


    def _parse_data_file(self, fn):

        ret = {}

        header = [
            'CHROM', 'POS', 'REF', 'ALT', 'ENST', 'GENE', 'LOC', 'CSN', 'CLASS', 'RefKnown5', 'RefKnown3',
            'AltKnown', 'AltHighest5', 'RefHighest5', 'AltHighest3', 'RefHighest3', 'Boundary5', 'Boundary3'
        ]

        for line in open(fn):
            line = line.strip()
            if line == '' or line[0] == '#' or line.upper().startswith('CHROM'):
                continue
            cols = line.split('\t')

            values = OrderedDict()
            for i in range(len(header)):
                values[header[i]] = cols[i]

            key = tuple(cols[:4])
            ret[key] = values

        return ret


    def get_scores(self, key):

        if not (len(key[2]) == 1 and len(key[3]) == 1):
            return {}

        if key not in self.data:
            return {}

        ret = {}

        for k in [
            'RefKnown5',
            'RefKnown3',
            'AltKnown',
            'AltHighest5',
            'RefHighest5',
            'AltHighest3',
            'RefHighest3',
            'Boundary5',
            'Boundary3'
        ]:
            ret[k] = self.data[key][k]

        ret['PI5'] = self._pi5(key)
        ret['PI3'] = self._pi3(key)
        ret.update(self._process_maxentscan_data(key))
        return ret


    def _pi5(self, key):

        d = self.data[key]

        Alt5 = '.' if d['AltKnown'] != '.' else d['AltHighest5']

        if Alt5 == '.':
            return '.'

        AltHighest5 = float(d['AltHighest5'])
        RefHighest5 = float(d['RefHighest5'])

        if AltHighest5 <= 0:
            return 0
        if RefHighest5 <= 0:
            return 100

        tmp = int(round(100 * (AltHighest5 - min(RefHighest5, AltHighest5)) / RefHighest5))
        return min(100, tmp)


    def _pi3(self, key):

        d = self.data[key]

        Alt3 = '.' if d['AltKnown'] != '.' else d['AltHighest3']

        if Alt3 == '.':
            return '.'

        AltHighest3 = float(d['AltHighest3'])
        RefHighest3 = float(d['RefHighest3'])

        if AltHighest3 <= 0:
            return 0
        if RefHighest3 <= 0:
            return 100

        tmp = int(round(100 * (AltHighest3 - min(RefHighest3, AltHighest3)) / RefHighest3))
        return min(100, tmp)


    def _process_maxentscan_data(self, key):

        ret = {}
        ret['RefKnown'] = self._ref_known(key)
        ret['SpliceSiteScore'] = ret['RefKnown']
        ret['SpliceSiteType'] = 'spliceSiteRegion' if ret['SpliceSiteScore'] != '.' else '.'
        ret['PercentReduction'] = self._percent_reduction(ret['RefKnown'], key)
        d = self.data[key]
        if d['AltKnown'] != '.':
            ret['MAX5'] = '.'
            ret['MAX3'] = '.'
        else:
            ret['MAX5'] = d['AltHighest5']
            ret['MAX3'] = d['AltHighest3']
        return ret


    def _ref_known(self, key):

        d = self.data[key]

        if d['AltKnown'] == '.':
            return '.'

        if '+' in d['CSN']:
            return d['RefKnown5']

        if '-' in d['CSN']:
            return d['RefKnown3']

        cval = int(d['CSN'][2:d['CSN'].find('>') - 1])
        if int(d['Boundary5']) - 2 <= cval <= int(d['Boundary5']):
            return d['RefKnown5']
        else:
            if int(d['Boundary3']) <= cval <= int(d['Boundary3']) + 2:
                return d['RefKnown3']
            else:
                return 'check'


    def _percent_reduction(self, ref_known, key):

        d = self.data[key]

        if d['AltKnown'] == '.':
            return '.'

        y = float(ref_known)

        if float(d['AltKnown']) < 0:
            z = 0
        else:
            if float(d['AltKnown']) > float(ref_known):
                z = float(ref_known)
            else:
                z = float(d['AltKnown'])

        x = float(ref_known) - z

        return round(100 * x / y, 2)








