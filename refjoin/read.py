# Copyright 2019 by Kent Kawashima.  All rights reserved.

from collections import namedtuple
import numpy as np

Sequence = namedtuple('Sequence', ['name', 'description', 'sequence'])

def fasta_to_seqobj_list(path):
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
                seq += line
        if seq:
            seq_list.append(Sequence(name, desc, seq))
    return seq_list
