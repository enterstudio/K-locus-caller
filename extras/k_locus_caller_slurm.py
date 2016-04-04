#!/usr/bin/env python
'''
K locus caller - SLURM job generator

This tool generates SLURM jobs so the K locus caller program can be run in parallel. All instances
of the script will append to the same output table.

Author: Ryan Wick
email: rrwick@gmail.com
'''

from __future__ import print_function
from __future__ import division
import sys
import argparse
import os

def main():
    '''
    Script execution starts here.
    '''
    args = get_arguments()
    check_files_exist(args.assembly + [args.k_refs])
    fix_paths(args)
    script_path = find_script(args)
    for assembly in args.assembly:
        assembly_name = os.path.splitext(os.path.basename(assembly))[0]
        cmd = '#!/bin/bash'
        cmd += '\n#SBATCH -p sysgen'
        cmd += '\n#SBATCH --job-name=k_locus_caller_' + assembly_name
        cmd += '\n#SBATCH --ntasks=1'
        cmd += '\n#SBATCH --mem-per-cpu=' + args.memory
        cmd += '\n#SBATCH --time=' + args.walltime
        cmd += '\nmodule load Python/2.7.10-vlsci_intel-2015.08.25'
        cmd += '\nBLAST+/2.2.30-vlsci_intel-2015.08.25-Python-2.7.10'
        cmd += '\n' + script_path
        cmd += ' --assembly ' + assembly
        cmd += ' --k_refs ' + args.k_refs
        cmd += ' --out ' + args.out
        if args.verbose:
            cmd += ' --verbose'
        if args.no_seq_out:
            cmd += ' --no_seq_out'
        if args.start_end_margin != 10:
            cmd += ' --start_end_margin ' + str(args.start_end_margin)
        if args.min_gene_cov != 90.0:
            cmd += ' --min_gene_cov ' + str(args.min_gene_cov)
        if args.min_gene_id != 50.0:
            cmd += ' --min_gene_id ' + str(args.min_gene_id)
        if args.low_gene_id != 95.0:
            cmd += ' --low_gene_id ' + str(args.low_gene_id)
        if args.min_assembly_piece != 100:
            cmd += ' --min_assembly_piece ' + str(args.min_assembly_piece)
        if args.gap_fill_size != 100:
            cmd += ' --gap_fill_size ' + str(args.gap_fill_size)
        print(cmd)
        print()
        # os.system('echo "' + cmd + '" | sbatch')
        print()

def get_arguments():
    '''
    Specifies the command line arguments required by the script.
    '''
    parser = argparse.ArgumentParser(description='K locus caller - SLURM job generator',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--walltime', type=str, required=False,
                        help='wall time (default 0-0:30 = 30 min)', default='0-0:30')
    parser.add_argument('--memory', type=str, required=False,
                        help='mem (default 4096 = 4gb)', default='4096')
    parser.add_argument('--script', type=str, required=False, default='find automatically',
                        help='path to k_locus_caller.py')
    parser.add_argument('-a', '--assembly', nargs='+', type=str, required=True,
                        help='Fasta file(s) for assemblies')
    parser.add_argument('-k', '--k_refs', type=str, required=True,
                        help='Genbank file with reference K loci')
    parser.add_argument('-o', '--out', type=str, required=False, default='./k_locus_results',
                        help='Output directory/prefix')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Display detailed information about each assembly in stdout')
    parser.add_argument('--no_seq_out', action='store_true',
                        help='Suppress output files of sequences matching K locus')
    parser.add_argument('--start_end_margin', type=int, required=False, default=10,
                        help='Missing bases at the ends of K locus allowed in a perfect match.')
    parser.add_argument('--min_gene_cov', type=float, required=False, default=90.0,
                        help='minimum required %% coverage for genes')
    parser.add_argument('--min_gene_id', type=float, required=False, default=50.0,
                        help='minimum required %% identity for genes')
    parser.add_argument('--low_gene_id', type=float, required=False, default=95.0,
                        help='genes with a %% identity below this value will be flagged as low '
                             'identity')
    parser.add_argument('--min_assembly_piece', type=int, required=False, default=100,
                        help='minimum K locus matching assembly piece to return')
    parser.add_argument('--gap_fill_size', type=int, required=False, default=100,
                        help='when separate parts of the assembly are found within this distance, '
                             'they will be merged')
    return parser.parse_args()

def quit_with_error(message): # type: (str) -> None
    '''
    Displays the given message and ends the program's execution.
    '''
    print('Error:', message, file=sys.stderr)
    sys.exit(1)

def check_files_exist(filenames): # type: (list[str]) -> bool
    '''
    Checks to make sure each file in the given list exists.
    '''
    for filename in filenames:
        check_file_exists(filename)

def check_file_exists(filename): # type: (str) -> bool
    '''
    Checks to make sure the single given file exists.
    '''
    if not os.path.isfile(filename):
        quit_with_error('could not find ' + filename)

def find_script(args):
    '''
    Returns the location of k_locus_caller.py. If the user gave a location, this script just checks
    to make sure it's valid. If the user didn't give a location, it checks in this script's
    directory, the current working directory and all parent directories of those two locations.
    '''
    if args.script != 'find automatically':
        args.script = os.path.abspath(args.script)
        check_file_exists(args.script)
        return args.script
    script_dir = os.path.dirname(os.path.realpath(__file__))
    while not os.path.isfile(os.path.join(script_dir, 'k_locus_caller.py')):
        script_dir = os.path.dirname(script_dir)
    if not os.path.isfile(os.path.join(script_dir, 'k_locus_caller.py')):
        script_dir = os.getcwd()
        while not os.path.isfile(os.path.join(script_dir, 'k_locus_caller.py')):
            script_dir = os.path.dirname(script_dir)
    if not os.path.isfile(os.path.join(script_dir, 'k_locus_caller.py')):
        quit_with_error('Could not find k_locus_caller.py')
    else:
        return os.path.join(script_dir, 'k_locus_caller.py')

def fix_paths(args):
    '''
    Changes the paths given by the user to absolute paths, which are easier to work with later.
    Also creates the output directory, if necessary.
    '''
    args.assembly = [os.path.abspath(x) for x in args.assembly]
    args.k_refs = os.path.abspath(args.k_refs)
    if args.out[-1] == '/':
        args.out += 'k_locus_results' 
    args.out = os.path.abspath(args.out)
    out_dir = os.path.dirname(args.out)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

if __name__ == '__main__':
    main()
