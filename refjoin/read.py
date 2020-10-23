# Copyright 2019 by Kent Kawashima.  All rights reserved.

from collections import namedtuple
import numpy as np

Sequence = namedtuple('Sequence', ['name', 'description', 'sequence'])

def fasta_alignment_to_seqobj_list(path):
    """ Read FASTA file for aligned sequences to a list of Sequence objects. 
    This function assumes sequences are aligned so that sequence lengths are 
    the same across all samples. In addition, an empty FASTA file also raises
    an error.

    Parameter
    ---------
    path: str
        A path to FASTA file for aligned sequences.

    Return
    ------
    seq_list: list
        A list of Sequence objects. A Sequence represents an allele in the 
        alignment with attributes of its name, description and sequence.
        
    """
    name = ''
    desc = ''
    seq = ''
    seq_list = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    seq_list.append(Sequence(name, desc, seq))
                try:
                    name, desc = line[1:].split(maxsplit=2)
                except ValueError:
                    name = line[1:]
                seq = ''
            else:
                seq += line.upper()
        if seq:
            seq_list.append(Sequence(name, desc, seq))

    # If sequence length in the alignment is different, 
    # raise Error
    seq_len_list = [len(seq.sequence) for seq in seq_list]
    if min(seq_len_list) != max(seq_len_list):
        raise ValueError('Different sequence size is found: '\
                        f'{min(seq_len_list)} != {max(seq_len_list)} in {path}.')
    if len(seq_list) == 0:
        raise AssertionError(f'{path} is empty.')

    return seq_list
