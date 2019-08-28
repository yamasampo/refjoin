# Copyright 2019 by Kent Kawashima.  All rights reserved.

import numpy as np
from .constants import CODON_INTCODE, INTCODE_CODON
from .read import Sequence, fasta_to_seqobj_list

def transpose_codon_aln_matrix(seq_list):
    # Assumes all sequences have the same length
    t_seq_list = []
    for i in range(0, len(seq_list[0].sequence), 3):
        
        col_seq = [
            CODON_INTCODE[s.sequence[i: i+3]]
            if s.sequence[i: i+3] in CODON_INTCODE.keys() else 65
            for s in seq_list
        ]
        t_seq_list.append(col_seq)
    return np.array(t_seq_list)

def refjoin_codon_aln_matrices(tmat1, tmat2, ref1_pos, ref2_pos, gap_code=196):
    t_seq_list = []
    # Test if ungapped ref sequence are identical
    refseq1 = [v for v in tmat1[:, ref1_pos] if v != gap_code]
    refseq2 = [v for v in tmat2[:, ref2_pos] if v != gap_code]
    if refseq1 != refseq2:
        raise ValueError(f'Ungapped reference sequence in tmat1 (pos: {ref1_pos}) '
                         f'does not match the ungapped reference sequence in tmat2 '
                         f'(pos: {ref2_pos}): {refseq1} != {refseq2}')
    # Join row by row
    c1 = 0
    c2 = 0
    gapped1 = tmat1[:, ref1_pos]
    gapped2 = tmat2[:, ref2_pos]
    len1 = len(tmat1[0])
    len2 = len(tmat2[0])
    joined_t_seq_list = []
    while c1 < len(tmat1) and c2 < len(tmat2):
        if gapped1[c1] == gap_code and gapped2[c2] == gap_code:
            joined1 = tmat1[c1].tolist() + [gap_code]*len2
            joined_t_seq_list.append(joined1)
            joined2 = [gap_code]*len1 + tmat2[c2].tolist()
            joined_t_seq_list.append(joined2)
            c1 += 1
            c2 += 1
        elif gapped1[c1] == gap_code and gapped2[c2] != gap_code:
            joined = tmat1[c1].tolist() + [gap_code]*len2
            joined_t_seq_list.append(joined)
            c1 += 1
        elif gapped1[c1] != gap_code and gapped2[c2] == gap_code:
            joined = [gap_code]*len1 + tmat2[c2].tolist()
            joined_t_seq_list.append(joined)
            c2 += 1
        elif gapped1[c1] == gapped2[c2]:
            joined = tmat1[c1].tolist() + tmat2[c2].tolist()
            joined_t_seq_list.append(joined)
            c1 += 1
            c2 += 1
    if c1 < len(tmat1):
        for i in range(c1, len(tmat1)):
            joined = tmat1[i].tolist() + [gap_code]*len2
            joined_t_seq_list.append(joined)
    elif c2 < len(tmat2):
        for i in range(c2, len(tmat2)):
            joined = [gap_code]*len1 + tmat2[i].tolist()
            joined_t_seq_list.append(joined)
    return np.array(joined_t_seq_list).T

def refjoin_codon_alignments(aln_path1, aln_path2, ref1_pos, ref2_pos, output_path):
    aln1_seq_list = fasta_to_seqobj_list(aln_path1)
    aln2_seq_list = fasta_to_seqobj_list(aln_path2)
    
    tmat1 = transpose_codon_aln_matrix(aln1_seq_list)
    tmat2 = transpose_codon_aln_matrix(aln2_seq_list)
    
    joined_aln_mat = refjoin_codon_aln_matrices(tmat1, tmat2, ref1_pos, ref2_pos)
    
    seq_list = []
    for i, s in enumerate(aln1_seq_list + aln2_seq_list):
        seq = ''.join([INTCODE_CODON[v] for v in joined_aln_mat[i]])
        seq_list.append(Sequence(s.name, s.description, seq))
    
    # Write to file
    with open(output_path, 'w') as f:
        for s in seq_list:
            if s.description:
                print(f'>{s.name} {s.description}', file=f)
            else:
                print(f'>{s.name}', file=f)
            print(f'{s.sequence}', file=f)
            
    return fasta_to_seqobj_list(output_path)
