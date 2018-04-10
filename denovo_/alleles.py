
def count(samfile, var_key):

    [chrom, pos, ref, alt] = var_key

    if is_substitution(ref, alt) or is_deletion(ref, alt):
        start, end = pos, pos
    elif is_insertion(ref, alt):
        start, end = pos, pos + 1
    else:
        start, end = pos, pos + len(ref) - 1

    tc = 0
    tr = 0
    for read in samfile.fetch(chrom, start, end + 1):

        # Filter out optical and PCR duplicate reads
        if read.is_duplicate:
            continue

        tc += 1

        # If substitution
        if is_substitution(ref, alt):
            if count_as_substitution(read, pos, alt):
                tr += 1
            continue

        # If deletion
        if is_deletion(ref, alt):
            if count_as_deletion(read, pos, ref):
                tr += 1
            continue

        # If insertion
        if is_insertion(ref, alt):
            if count_as_insertion(read, pos):
                tr += 1
            continue

        # If complex indel
        if count_as_complex(read, pos, ref):
            tr += 1

    return tc, tr


def count_as_substitution(read, pos, alt):

    for x in read.get_aligned_pairs(with_seq=True):
        if x[1] == pos:
            if x[0] is None:
                return False
            return x[2] == alt

    return False


def count_as_deletion(read, pos, ref):

    deletion_start = pos + 1
    deletion_end = pos + len(ref) - 1

    for d in deletion_blocks(read):
        if d[0] <= deletion_start and deletion_end <= d[1]:
            return True
    return False


def count_as_insertion(read, pos):

    for d in insertion_blocks(read):
        if d[0] == pos:
            return True
    return False


def count_as_complex(read, pos, ref):

    if count_as_deletion(read, pos, ref):
        return True

    for p in range(pos, pos + len(ref) - 1):
        if count_as_insertion(read, p):
            return True

    return False


def deletion_blocks(read):

    ret = []
    blocks = read.get_blocks()
    for i in range(1, len(blocks)):
        if blocks[i][0] > blocks[i - 1][1]:
            ret.append((blocks[i - 1][1], blocks[i][0] - 1))
    return ret


def insertion_blocks(read):

    ret = []
    blocks = read.get_blocks()
    for i in range(1, len(blocks)):
        if blocks[i][0] == blocks[i - 1][1]:
            ret.append((blocks[i - 1][1] - 1, blocks[i][0]))
    return ret


def is_substitution(ref, alt):

    return len(ref) == 1 and len(alt) == 1


def is_deletion(ref, alt):

    return ref[0] == alt[0] and len(ref) > 1 and len(alt) == 1


def is_insertion(ref, alt):

    return ref[0] == alt[0] and len(ref) == 1 and len(alt) > 1

