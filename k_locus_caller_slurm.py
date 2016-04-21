#!/usr/bin/env python
'''
K locus caller - SLURM job generator

This is a tool for parallelising the work of K locus caller by making one SLURM job per assembly.

Author: Ryan Wick
email: rrwick@gmail.com
'''

from __future__ import print_function
from __future__ import division
import sys
import argparse
import os
import k_locus_caller



def main():
    '''Script execution starts here.'''
    args = get_arguments()
    k_locus_caller.check_for_blast()
    k_locus_caller.check_files_exist(args.assembly + [args.k_refs])

    if args.script:
        if args.script.endswith('k_locus_caller.py'):
            script_path = args.script
        else:
            script_path = os.path.join(args.script, 'k_locus_caller.py')
        k_locus_caller.check_file_exists(script_path)
    else:
        script_path_cwd = os.path.join(os.getcwd(), 'k_locus_caller.py')
        if os.path.isfile(script_path_cwd):
            script_path = script_path_cwd
        else:
            script_path_file_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                      'k_locus_caller.py')
            if os.path.isfile(script_path_file_directory):
                script_path = script_path_file_directory
            else:
                k_locus_caller.quit_with_error('could not find k_locus_caller.py')

    main_parser = k_locus_caller.get_argument_parser()
    using_start_end_margin = (args.start_end_margin != main_parser.get_default('start_end_margin'))
    using_min_gene_cov = (args.min_gene_cov != main_parser.get_default('min_gene_cov'))
    using_min_gene_id = (args.min_gene_id != main_parser.get_default('min_gene_id'))
    using_low_gene_id = (args.low_gene_id != main_parser.get_default('low_gene_id'))
    using_min_assembly_piece = (args.min_assembly_piece != main_parser.get_default('min_assembly_piece'))
    using_gap_fill_size = (args.gap_fill_size != main_parser.get_default('gap_fill_size'))
    using_low_gene_id = (args.low_gene_id != main_parser.get_default('low_gene_id'))

    for assembly_filename in args.assembly:
        assembly_name = simple_assembly_name(assembly_filename)
        assembly_filename = os.path.abspath(assembly_filename)

        # Build the k_locus_caller.py command.
        k_locus_caller_command = script_path
        k_locus_caller_command += ' --assembly ' + os.path.abspath(assembly_filename)
        k_locus_caller_command += ' --k_refs ' + os.path.abspath(args.k_refs)
        if args.out.endswith('/'):
            k_locus_caller_command += ' --out ' + os.path.abspath(args.out) + '/'
        else:
            k_locus_caller_command += ' --out ' + os.path.abspath(args.out)
        if args.verbose:
            k_locus_caller_command += ' --verbose'
        if args.no_seq_out:
            k_locus_caller_command += ' --no_seq_out'
        if using_start_end_margin:
            k_locus_caller_command += ' --start_end_margin ' + str(args.start_end_margin)
        if using_min_gene_cov:
            k_locus_caller_command += ' --min_gene_cov ' + str(args.min_gene_cov)
        if using_min_gene_id:
            k_locus_caller_command += ' --min_gene_id ' + str(args.min_gene_id)
        if using_low_gene_id:
            k_locus_caller_command += ' --low_gene_id ' + str(args.low_gene_id)
        if using_min_assembly_piece:
            k_locus_caller_command += ' --min_assembly_piece ' + str(args.min_assembly_piece)
        if using_gap_fill_size:
            k_locus_caller_command += ' --gap_fill_size ' + str(args.gap_fill_size)
        if using_low_gene_id:
            k_locus_caller_command += ' --low_gene_id ' + str(args.low_gene_id)

        # Build the SLURM script.
        slurm = []
        slurm.append('#!/bin/bash')
        slurm.append('#SBATCH -p sysgen')
        slurm.append('#SBATCH --job-name=k_locus_caller_' + assembly_name)
        slurm.append('#SBATCH --ntasks=1')
        slurm.append('#SBATCH --mem-per-cpu=' + args.memory)
        slurm.append('#SBATCH --time=' + args.walltime)
        slurm.append('module load Python/2.7.10-vlsci_intel-2015.08.25')
        slurm.append('module load BLAST+/2.2.30-vlsci_intel-2015.08.25-Python-2.7.10')
        slurm.append(k_locus_caller_command)

        slurm_string = '\n'.join(slurm)

        print(slurm_string)
        print()
        os.system('echo "' + slurm_string + '" | sbatch')
        print()

def get_arguments():
    '''Specifies the command line arguments required by the script.'''
    main_parser = k_locus_caller.get_argument_parser()

    parser = argparse.ArgumentParser(description='K locus caller - SLURM job generator',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--walltime', type=str, required=False,
                        help='wall time', default='0-0:10')
    parser.add_argument('--memory', type=str, required=False,
                        help='memory (Megabytes)', default='2048')
    parser.add_argument('--script', type=str, required=False,
                        help="path to k_locus_caller.py, if not in current directory or this "
                             "script's directory")
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
    parser.add_argument('--start_end_margin', type=int, required=False,
                        default=main_parser.get_default('start_end_margin'),
                        help='Missing bases at the ends of K locus allowed in a perfect match.')
    parser.add_argument('--min_gene_cov', type=float, required=False,
                        default=main_parser.get_default('min_gene_cov'),
                        help='minimum required %% coverage for genes')
    parser.add_argument('--min_gene_id', type=float, required=False,
                        default=main_parser.get_default('min_gene_id'),
                        help='minimum required %% identity for genes')
    parser.add_argument('--low_gene_id', type=float, required=False,
                        default=main_parser.get_default('low_gene_id'),
                        help='genes with a %% identity below this value will be flagged as low '
                             'identity')
    parser.add_argument('--min_assembly_piece', type=int, required=False,
                        default=main_parser.get_default('min_assembly_piece'),
                        help='minimum K locus matching assembly piece to return')
    parser.add_argument('--gap_fill_size', type=int, required=False,
                        default=main_parser.get_default('gap_fill_size'),
                        help='when separate parts of the assembly are found within this distance, '
                             'they will be merged')
    return parser.parse_args()


def simple_assembly_name(assembly_filename):
    '''
    Returns just the assembly name, without a path or file extension.
    '''
    assembly_filename = os.path.basename(assembly_filename)
    assembly_filename = rchop(assembly_filename, '.fasta')
    assembly_filename = rchop(assembly_filename, '.FASTA')
    assembly_filename = rchop(assembly_filename, '.fa')
    return assembly_filename

def rchop(thestring, ending):
    '''
    Returns thestring with ending removed from the end.
    http://stackoverflow.com/questions/3663450/python-remove-substring-only-at-the-end-of-string
    '''
    if thestring.endswith(ending):
        return thestring[:-len(ending)]
    return thestring

if __name__ == '__main__':
    main()
