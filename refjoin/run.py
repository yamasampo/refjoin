""" Python script calling refjoin to join two sequence alignments. """

import os
import argparse
import configparser
from warnings import warn
from datetime import datetime
import alignmentrs as al
from refjoin import __version__
from refjoin.codon_alignment import refjoin_codon_alignments
from refjoin.nucleotide_alignment import refjoin_alignments

def main(
        mode:str,
        aln_dir1:str, 
        aln_dir2:str,
        file_pair_list_path:str,
        key_seq_prefix1:str,
        key_seq_prefix2:str,
        out_aln_dir:str,
        out_suffix:str,
        out_info_dir:str
        ):

    start = datetime.now().isoformat()
    
    # Specify in which frame alignments are joined.
    if mode == 'codon':
        refjoin_func = refjoin_codon_alignments
    elif mode == 'nucleotide':
        refjoin_func = refjoin_alignments
    else:
        raise ValueError(f'Unknown mode {mode} is found. '\
                          'Please specify "codon" or "nucleotide."')
    
    # Make output directory
    os.mkdir(out_aln_dir)
    os.mkdir(out_info_dir)

    # Output read control file
    out_ctl_path = os.path.join(out_info_dir, '0.ctl_read.txt')
    not_found_item_path = os.path.join(out_info_dir, '0.notfound_aln_files.list')
    to_control_file(
        out_ctl_path, start, 
        mode=mode,
        aln_dir1=aln_dir1, 
        aln_dir2=aln_dir2,
        file_pair_list_path=file_pair_list_path,
        key_seq_prefix1=key_seq_prefix1,
        key_seq_prefix2=key_seq_prefix2,
        out_aln_dir=out_aln_dir,
        out_suffix=out_suffix,
        out_info_dir=out_info_dir
    )

    # Call refjoin to all given pairs of alignments
    apply_refjoin_alignments(
        refjoin_func,
        aln_dir1,
        aln_dir2,
        file_pair_list_path, 
        key_seq_prefix1,
        key_seq_prefix2,
        out_aln_dir,
        out_suffix,
        not_found_item_path
    )

def apply_refjoin_alignments(
        refjoin_func,
        aln_dir1:str,
        aln_dir2:str,
        file_pair_list_path:str, 
        key_seq_prefix1:str,
        key_seq_prefix2:str,
        out_aln_dir:str,
        out_suffix:str,
        not_found_item_path:str
        ):
    with open(file_pair_list_path, 'r') as f:
        file_pairs = [
            l[:-1].split('\t') for l in f 
            if not l.startswith('itemnum') and not l.startswith('/*')
        ]
    poly_absense = []

    for fname1, fname2 in file_pairs:
        aln_name = fname1.split('.')[0]

        out_path = os.path.join(out_aln_dir, f'{aln_name}{out_suffix}')
        if os.path.isfile(out_path):
            raise FileExistsError(out_path)

        aln_path1 = os.path.join(aln_dir1, fname1)

        if not os.path.isfile(aln_path1):
            poly_absense.append(aln_path1)
            continue

        aln_path2 = os.path.join(aln_dir2, fname2)
        if not os.path.isfile(aln_path2):
            poly_absense.append(aln_path2)
            continue
        
        # Read alignment file
        aln1 = al.fasta_file_to_alignment(aln_path1, 'msye alignment')

        # Find key sequence position
        ref_pos1_list = [
            i for i, sample in enumerate(aln1.sample_ids) 
            if sample.startswith(key_seq_prefix1)
        ]
        assert len(ref_pos1_list) == 1
        ref_pos1 = ref_pos1_list[0]

        # Read alignment file
        aln2 = al.fasta_file_to_alignment(aln_path2, 'msye alignment')

        # Find key sequence position
        ref_pos2_list = [
            i for i, sample in enumerate(aln2.sample_ids) 
            if sample.startswith(key_seq_prefix2)
        ]
        assert len(ref_pos2_list) == 1
        ref_pos2 = ref_pos2_list[0]
        
        refjoin_func(aln_path1, aln_path2, ref_pos1, ref_pos2, out_path)

    flist = to_filelist(out_aln_dir)
    print(len(flist), 'files saved.')

    with open(not_found_item_path, 'w') as f:
        print('itemnum: {}'.format(len(poly_absense)), file=f)
        print('\n'.join(poly_absense), file=f)
    
def to_filelist(dir_path):
    '''Return a file list in a given directory'''
    l1 = os.listdir(dir_path)
    l2= [a for a in l1 if not a.startswith('.')]
    flist = [a for a in l2 if not a.startswith('0')]

    with open(os.path.join(dir_path, '0.filelist'), 'w') as f:
        print('itemnum: '+str(len(flist)), file=f)
        print('\n'.join(flist), file=f)
    return flist

def to_control_file(output_file, start_date, **kwargs):
    lines = [f'{k} = {v}' for k, v in kwargs.items()]

    with open(output_file, 'w') as f:
        print(f'Run date: {start_date}.', file=f)
        print(f'refjoin package version: {__version__}\n', file=f)
        print('Listed parameters are used for the run.\n', file=f)
        print('\n'.join(lines), file=f)

def read_control_file(control_file):
    # Initialize ConfigParser object
    config = configparser.ConfigParser(
        strict=True,
        comment_prefixes=('#', ';'),
        inline_comment_prefixes=('#', ';')
    )

    # Parse control file
    paths = config.read(control_file)

    # Check number of read control files.
    if len(paths) == 0:
        raise FileNotFoundError(
        f'Specified control file, {control_file}, is not found.')
    elif len(paths) > 1:
        raise TypeError(f'Iterable {type(control_file)} is given as a control '\
            'file. Only one control file is supported.')

    # Check sections. Only 'REQUIRED' and 'OPTIONAL' sections will be used.
    assert 'REQUIRED' in config.sections(), \
        f'REQUIRED section is not found in {control_file}.'
    expected_sections = ['REQUIRED']
    not_expected_sections = [
        s for s in config.sections() if s not in expected_sections]
    if len(not_expected_sections) > 1:
        msg = f'Unexpected sections, {", ".join(not_expected_sections)}, '\
              'were found. These are not used in '\
              'the analysis. If you wish to include in the analysis, please '\
              'specify in "REQUIRED" or "OPTIONAL" sections'
        warn(msg)

    flattened  = {
        opt: v
        for s in expected_sections for opt, v in config[s].items()
    }

    return flattened

if __name__ == '__main__':
    desc = f'refjoin (version {__version__}): Join two alignments by referring '\
        'to the identical sequences as a template.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "ctl_file", 
        help="A path to control file where input arguments are listed."
    )
    args = parser.parse_args()
    params = read_control_file(args.ctl_file)

    print('\nINPUT')
    for k, v in params.items():
        print(f'{k}: {v}')

    main(
        params['mode'],
        params['aln_dir1'], 
        params['aln_dir2'],
        params['file_pair_list_path'],
        params['key_seq_prefix1'],
        params['key_seq_prefix2'],
        params['out_aln_dir'],
        params['out_suffix'],
        params['out_info_dir']
    )
