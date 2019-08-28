# Copyright 2019 by Kent Kawashima.  All rights reserved.

import numpy as np
from .read import Sequence, fasta_to_seqobj_list

# Create a list of columns with keeping sample orders
# Transpose sequence matrix
def transpose_aln_matrix(seq_list):
    # Assumes all sequences have the same length
    t_seq_list = []
    for i in range(len(seq_list[0].sequence)):
        col_seq = [s.sequence[i] for s in seq_list]
        t_seq_list.append(col_seq)
        
    return np.array(t_seq_list)
    

def refjoin_aln_matrices(
    tmat1, # Output from transpose_aln_matrix function
    tmat2, 
    ref1_pos, # Position of reference sequence in sample 1 alignment. If the first sample is reference, ref1_pos=0. 
    ref2_pos): # 0 based index. if -1 will work is not tested yet
    
    t_seq_list = []
    # Test if ungapped ref sequence are identical
    # Check if reference sequences have the same sequence other than gap
    # If we remove the gaps, ref seqs should be identical. length is not necessary
    # This "reference" does not mean Drosophila reference genome, 
    # does mean a template sequence shared between alignments.
    refseq1 = ''.join(tmat1[:, ref1_pos].tolist()).replace('-', '')
    refseq2 = ''.join(tmat2[:, ref2_pos].tolist()).replace('-', '')
    if refseq1 != refseq2:
        raise ValueError(f'Ungapped reference sequence in tmat1 (pos: {ref1_pos}) '
                         f'does not match the ungapped reference sequence in tmat2 '
                         f'(pos: {ref2_pos}): {refseq1} != {refseq2}')
    # Join row by row
    
    c1 = 0 # position counter
    c2 = 0
    gapped1 = tmat1[:, ref1_pos] # ref seq 1 with gaps
    gapped2 = tmat2[:, ref2_pos] # ref seq 2 with gaps
    len1 = len(tmat1[0]) # sample length of alignment
    len2 = len(tmat2[0]) # sample length of alignment
    joined_t_seq_list = []
    
    while c1 < len(tmat1) and c2 < len(tmat2):
        print(c1, c2)
        # If both ref seqs have gaps, then sample columns are separately appended to the output.
        if gapped1[c1] == '-' and gapped2[c2] == '-':
            joined1 = tmat1[c1].tolist() + ['-']*len2 # for sample 1 column
            joined_t_seq_list.append(joined1)
            joined2 = ['-']*len1 + tmat2[c2].tolist() # for sample 2 column
            joined_t_seq_list.append(joined2)
            c1 += 1
            c2 += 1
            
        # The case where ref seq1 has gap but ref seq2 does not.
        elif gapped1[c1] == '-' and gapped2[c2] != '-':
            joined = tmat1[c1].tolist() + ['-']*len2 # Insert a new column of all gaps into sample 2
            joined_t_seq_list.append(joined)
            c1 += 1 # position count of sample where ref seq has gap is incremented
            
        elif gapped1[c1] != '-' and gapped2[c2] == '-':
            joined = ['-']*len1 + tmat2[c2].tolist()
            joined_t_seq_list.append(joined)
            c2 += 1
            
        # if characters of ref seqs at the position is identical
        # then concatenate columns of the positions of two sample alignments
        # Even if there is gap in other samples, the postion will enter here
        elif gapped1[c1] == gapped2[c2]:
            joined = (tmat1[c1].tolist() # ['G', 'G', 'G', 'G']
                + tmat2[c2].tolist()) # ['G', 'G', 'G', 'G']
            joined_t_seq_list.append(joined)
            c1 += 1
            c2 += 1
            
    if c1 < len(tmat1):
        for i in range(c1, len(tmat1)):
            joined = tmat1[i].tolist() + ['-']*len2
            joined_t_seq_list.append(joined)
    elif c2 < len(tmat2):
        for i in range(c2, len(tmat2)):
            joined = ['-']*len1 + tmat2[i].tolist()
            joined_t_seq_list.append(joined)
            
    # Return a list of sequences (row)
    return np.array(joined_t_seq_list).T

def refjoin_alignments(aln_path1, aln_path2, ref1_pos, ref2_pos, output_path):
    # Read into a list of Sequence objects
    aln1_seq_list = fasta_to_seqobj_list(aln_path1)
    aln2_seq_list = fasta_to_seqobj_list(aln_path2)
    
    # Extract sequence and create transposed alignment sequence matrix
    tmat1 = transpose_aln_matrix(aln1_seq_list)
    tmat2 = transpose_aln_matrix(aln2_seq_list)
    
    # Join two alignments
    joined_aln_mat = refjoin_aln_matrices(tmat1, tmat2, ref1_pos, ref2_pos)
    
    # Connect sequence name and descriptions back to output object
    # Because while joining process, sample order is kept, 
    # seq name and descriptions can be just concatenated in the original order.
    seq_list = []
    for i, s in enumerate(aln1_seq_list + aln2_seq_list):
        seq = ''.join(joined_aln_mat[i].tolist())
        seq_list.append(Sequence(s.name, s.description, seq))
    
    # Write to file
    with open(output_path, 'w') as f:
        for s in seq_list:
            if s.description:
                print(f'>{s.name} {s.description}', file=f)
            else:
                print(f'>{s.name}', file=f)
            print(f'{s.sequence}', file=f)
            
    # Read back the written file into a list of Sequence objects
    return fasta_to_seqobj_list(output_path)
