import parsers
import helper
import maxentscan
import pysam


def run(options, version):

    # Print welcome message and information
    helper.welcome(version)

    # Read configuration file
    config = parsers.read_config_file(options.config)

    # Read variant data of the three individuals
    child_var_data = parsers.read_variant_file(options.child_var)
    mother_var_data = parsers.read_variant_file(options.mother_var)
    father_var_data = parsers.read_variant_file(options.father_var)

    # Create list of multiallelic calls in the child
    multiallelic_calls = helper.find_multiallelic_calls(child_var_data)

    # Connect to BAM files of the mother and the father
    mother_bam = pysam.AlignmentFile(options.mother_bam, "rb")
    father_bam = pysam.AlignmentFile(options.father_bam, "rb")

    # Connect to gnomAD and control database files
    gnomad_file = pysam.Tabixfile(config['GNOMAD_DATA_FILE'])
    control_data = parsers.read_custom_database_file(config['CONTROL_DATA_FILE'])

    # Read MaxEntScan data
    maxentscan_data = maxentscan.MaxEntScanData(config['MAXENTSCAN_DATA_FILE']) if config['MAXENTSCAN_DATA_FILE'] != '' else None

    # Initialize output files
    out_denovo = open('{}_denovo_candidates.txt'.format(options.output), 'w')
    helper.output_header(out_denovo, config['MAXENTSCAN_DATA_FILE'] != '')
    out_filtered = open('{}_filtered_out.txt'.format(options.output), 'w')

    # Initialize counters
    counter = 0
    counter_denovo = 0
    counter_filtered = 0

    # Initialize progress info
    helper.init_progress()

    # Iterate through the variants called in the child
    for var_key, data_list in child_var_data.iteritems():

        for data in data_list:

            counter += 1

            # Print progress info
            helper.print_progress(counter, len(child_var_data))

            res = helper.check(
                var_key,
                data,
                config,
                multiallelic_calls,
                mother_var_data,
                father_var_data,
                gnomad_file,
                control_data,
                mother_bam,
                father_bam
            )

            if type(res) == tuple:

                # MaxEntScan scores of the variant
                maxentscan_scores = maxentscan_data.get_scores(var_key) if maxentscan_data is not None else None

                helper.output(out_denovo, var_key, data, result[0], result[1], result[2], maxentscan_scores, '.')
                counter_denovo += 1

            else:

                out_filtered.write(
                    '\t'.join(
                        [
                            var_key[0],
                            var_key[1],
                            var_key[2],
                            var_key[3],
                            data['gene'],
                            data['csn'],
                            res
                        ]
                    ) + '\n'
                )
                counter_filtered += 1

    # Finalize progress info
    helper.finalize_progress()

    # Close output files
    out_denovo.close()
    out_filtered.close()

    # Print goodbye message and information
    helper.goodbye(counter_denovo, counter_filtered, options.output)
