#!env/bin/python

from optparse import OptionParser
from denovo_ import main

# Version
_version = '0.1.0'

# Command line argument parsing
descr = 'DeNovoFilter v'+_version
parser = OptionParser(usage='DeNovoFilter-{}/denovo <options>'.format(_version), version=_version, description=descr)
parser.add_option("--child_var", default=None, dest='child_var', action='store',help="Variants data file of the child")
parser.add_option("--mother_var", default=None, dest='mother_var', action='store',help="Variants data file of the mother")
parser.add_option("--father_var", default=None, dest='father_var', action='store',help="Variants data file of the father")
parser.add_option("--child_bam", default=None, dest='child_bam', action='store',help="BAM file of the child")
parser.add_option("--mother_bam", default=None, dest='mother_bam', action='store',help="BAM file of the mother")
parser.add_option("--father_bam", default=None, dest='father_bam', action='store',help="BAM file of the father")
parser.add_option("--config", default=None, dest='config', action='store',help="Configuration file")
parser.add_option("--full_details", default=False, dest='full_details', action='store_true', help="Output full details")
parser.add_option('-o', "--output", default=None, dest='output', action='store', help="Output filename prefix")
(options, args) = parser.parse_args()

main.run(options, _version)

