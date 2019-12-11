nnk_table: {str: str} = {"CGT": "R", "CGG": "R", "AGG": "R",
                         "CTT": "L", "CTG": "L", "TTG": "L",
                         "TCT": "S", "TCG": "S", "AGT": "S",
                         "GCT": "A", "GCG": "A",
                         "GGT": "G", "GGG": "G",
                         "CCT": "P", "CCG": "P",
                         "ACT": "T", "ACG": "T",
                         "CAG": "Q",
                         "TAG": "q",
                         "GTT": "V", "GTG": "V",
                         "AAT": "N",
                         "GAT": "D",
                         "TGT": "C",
                         "GAG": "E",
                         "CAT": "H",
                         "ATT": "I",
                         "AAG": "K",
                         "ATG": "M",
                         "TTT": "F",
                         "TGG": "W",
                         "TAT": "Y"}

def verify_file_is_not_empty(file_path):
    import logging
    logger = logging.getLogger('main')
    # make sure that there are results and the file is not empty
    with open(file_path) as f:
        if len(f.read(10).strip()) == 0:
            # TODO: write error to a global error file
            msg = f'Input file is empty {file_path}'
            logger.error(msg)
            raise RuntimeError(msg)


def load_msa(msa_path):
    """
    :param msa_path: a path to an MSA file (in FASTA format)
    :return: a dictionary that maps each header (without ">" and rstriped()) to its corresponding sequence (rstriped())
             an int that represent the number of sequences
             an int that represent the length of the alignment
    """
    header_to_sequence = {}
    with open(msa_path) as f:
        for header in f:
            if not header.startswith('>'):
                raise TypeError('Illegal fasta file')
            # returns header without ">" !
            header_to_sequence[header[1:].rstrip()] = f.readline().rstrip()

    return header_to_sequence, len(header_to_sequence), len(header_to_sequence[header.rstrip()])
