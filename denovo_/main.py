import parsers
import helper
import filters


def run(options, version):

    # Print welcome message and information
    helper.welcome(version)

    # Read configuration file
    config = parsers.read_config_file(options.config)

    # Read input data
    data = helper.read_data(options, config)

    # Initialize output files
    out_included = open('{}_denovo_candidates.txt'.format(options.output), 'w')
    out_excluded = open('{}_filtered_out.txt'.format(options.output), 'w')

    # Decide whether to write MaxEntScan and ExAC columns
    write_maxentscan = config['MAXENTSCAN_DATA_FILE'] != ''
    write_exac = config['EXAC_DATA_FILE'] != ''

    # Write output file headers
    helper.output_header(out_included, write_maxentscan, write_exac)
    if options.full_details:
        helper.output_header(out_excluded, write_maxentscan, write_exac)
    else:
        helper.output_header_simplified(out_excluded)

    # Initialize Filters object
    filt = filters.Filters(options.full_details)

    # Initialize counters
    counter = 0
    counter_included = 0
    counter_excluded = 0

    # Initialize progress info
    helper.init_progress()

    # Iterate through the variants called in the child
    for var_key, variants in data['child_var'].iteritems():

        for variant in variants:
            counter += 1

            # Print progress info
            helper.print_progress(counter, len(data['child_var']))

            # MaxEntScan scores of the variant
            maxentscan_scores = data['maxentscan'].get_scores(var_key) if data['maxentscan'] is not None else None

            # ExAC column values of the variant
            if data['exac'] is None:
                exac_values = None
            else:
                gene = data['gene']
                exac_values = data['exac'][gene] if gene in data['exac'] else {}

            # Apply filters to the variant
            result = filt.apply_filters(var_key, variant, config, data)

            # Output result
            if result['filter'] == '.':
                helper.output(out_included, var_key, variant, result, maxentscan_scores, exac_values)
                counter_included += 1
            else:
                if options.full_details:
                    helper.output(out_excluded, var_key, variant, result, maxentscan_scores, exac_values)
                else:
                    helper.output_simplified(out_excluded, var_key, variant, result)
                counter_excluded += 1

    # Finalize progress info
    helper.finalize_progress()

    # Close output files
    out_included.close()
    out_excluded.close()

    # Print goodbye message and information
    helper.goodbye(counter_included, counter_excluded, options.output)


