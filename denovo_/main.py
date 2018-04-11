import parsers
import helper
import maxentscan
import pysam
import alleles


def run(options, version):

    # Print welcome message and information
    helper.welcome(version)

    # Read configuration file
    config = parsers.read_config_file(options.config)

    # Print info
    helper.print_info(options, config)

    # Read variant data of the three individuals
    child_var_data = parsers.read_variant_file(options.child_var)
    mother_var_data = parsers.read_variant_file(options.mother_var)
    father_var_data = parsers.read_variant_file(options.father_var)

    # Create list of multiallelic calls in the child
    multiallelic_calls = helper.find_multiallelic_calls(child_var_data)

    # Connect to BAM files of the mother and the father
    mother_bam = pysam.AlignmentFile(options.mother_bam, "rb")
    father_bam = pysam.AlignmentFile(options.father_bam, "rb")

    # Read gnomAD and control database files
    gnomad_data, gnomad_use_csn = parsers.read_gnomad_database_file(config['GNOMAD_DATA_FILE'])
    control_data = parsers.read_custom_database_file(config['CONTROL_DATA_FILE'])

    # Read MaxEntScan data
    maxentscan_data = maxentscan.MaxEntScanData(config['MAXENTSCAN_DATA_FILE']) if config['MAXENTSCAN_DATA_FILE'] != '' else None

    # Initialize output files
    out_denovo = open('{}_denovo_candidates.txt'.format(options.output), 'w')
    helper.output_header(out_denovo, config['MAXENTSCAN_DATA_FILE'] != '')
    out_filtered = open('{}_filtered.txt'.format(options.output), 'w')
    helper.output_header(out_filtered, config['MAXENTSCAN_DATA_FILE'] != '')

    # Initialize counters
    counter = 0
    counter_denovo = 0
    counter_filtered = 0

    # Initialize progress info
    helper.init_progress()

    # Iterate through the variants called in the child
    for var_key, data in child_var_data.iteritems():

        counter += 1

        # Print progress info
        helper.print_progress(counter, len(child_var_data))

        # Variant frequencies
        csn_key = (data['gene'], data['csn'])
        freqs = {
            'gnomad': gnomad_data[csn_key]['freq'] if gnomad_use_csn else gnomad_data[var_key]['freq'],
            'control': control_data[var_key]
        }

        # Allele counts in the parents
        parent_alleles = {}
        parent_alleles['mother_tc'], parent_alleles['mother_tr'] = alleles.count(mother_bam, var_key)
        parent_alleles['father_tc'], parent_alleles['father_tr'] = alleles.count(father_bam, var_key)

        # MaxEntScan scores of the variant
        maxentscan_scores = maxentscan_data.get_scores(var_key) if maxentscan_data is not None else None

        # Check if variant is a de novo candidate and output to file
        is_candidate, reason = helper.is_candidate(
                var_key,
                data,
                config,
                multiallelic_calls,
                mother_var_data,
                father_var_data,
                freqs,
                parent_alleles
        )
        if is_candidate:
            helper.output(out_denovo, var_key, data, freqs, parent_alleles, maxentscan_scores, '.')
            counter_denovo += 1
        else:
            helper.output(out_filtered, var_key, data, freqs, parent_alleles, maxentscan_scores, ';'.join(reason))
            counter_filtered += 1

    # Finalize progress info
    helper.finalize_progress()

    # Close output files
    out_denovo.close()
    out_filtered.close()

    # Print goodbye message and information
    helper.goodbye(counter_denovo, counter_filtered)
