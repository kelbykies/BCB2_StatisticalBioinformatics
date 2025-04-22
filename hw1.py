
def extract_info(input_file):
    """ Reads in the sequences and corresponding quality scores from the provided sequence file.
    """

    seq = ""
    quality = ""

    try:
        file = open(input_file, "r")
    except:
        print("File doesn't exist.")


    array = file.readlines()

    for i in range(0, len(array)):
        if array[i].startswith('@'):
            if not array[i-1].startswith('+'):
                seq = seq + array[i+1][:-1]
        if array[i].startswith('+'):
            quality = quality + array[i+1][:-1]

    base_quality_pairs = []
    for i in range(0, len(seq)):
        tpm = [seq[i], quality[i]]
        base_quality_pairs.append(tpm)

    file.close()
    return base_quality_pairs

print(extract_info('example_fastq'))