#!/usr/bin/env python
'''
K-locus caller

This is a tool which reports information about the K-type for Klebsiella genome assemblies.

As input, this script takes:
  * One or more Klebsiella genome assemblies
  * Reference nucleotide sequences for known K-types
  * A multi-fasta of protein sequences for genes in the known K-types

For each sample, it will report:
  * Which K-type best matches the sample
  * Information about the best K-type match (e.g. % identity, whether it was in one piece, etc.)
  * The DNA sequence(s) in the sample which align to the best K-type match
  * The location and identity of K-locus proteins in the sample

This information will help a user to decide whether their sample has a known K-type or a novel
one, and if novel, how similar it is to a known type.

Author: Ryan Wick
email: rrwick@gmail.com
'''

from __future__ import print_function
from __future__ import division
import argparse
import sys
import os
import subprocess
from distutils import spawn


def main():
    '''
    Script execution starts here.
    '''
    args = get_arguments()
    check_for_blast()
    check_files_exist(args.assembly)
    check_file_exists(args.ref_types)
    check_file_exists(args.ref_genes)
    k_refs_dict = load_k_locus_references(args.ref_types)
    for assembly in args.assembly:
        best_k = get_best_k_type_match(assembly, args.ref_types, k_refs_dict)
        print('Assembly: ' + assembly + ', best K-type match: ' + best_k.name) # TEMP
        assembly_pieces, ideal = get_assembly_pieces(assembly, best_k)



def get_arguments():
    '''
    Specifies the command line arguments required by the script.
    '''
    parser = argparse.ArgumentParser(description='K-locus caller')
    parser.add_argument('-a', '--assembly', nargs='+', type=str, required=True,
                        help='Fasta file(s) for Klebsiella assemblies')
    parser.add_argument('-r', '--ref_types', type=str, required=True,
                        help='Fasta file with reference K-locus nucleotide sequences')
    parser.add_argument('-g', '--ref_genes', type=str, required=True,
                        help='Fasta file reference gene protein sequences')
    return parser.parse_args()

def check_for_blast():
    '''
    Checks to make sure the required BLAST+ tools are available.
    '''
    if not spawn.find_executable('makeblastdb'):
        quit_with_error('could not find makeblastdb tool (part of BLAST+)')
    if not spawn.find_executable('blastn'):
        quit_with_error('could not find blastn tool (part of BLAST+)')
    if not spawn.find_executable('tblastn'):
        quit_with_error('could not find tblastn tool (part of BLAST+)')

def check_files_exist(filenames):
    '''
    Checks to make sure each file in the given list exists.
    '''
    for filename in filenames:
        check_file_exists(filename)

def check_file_exists(filename):
    '''
    Checks to make sure the single given file exists.
    '''
    if not os.path.isfile(filename):
        quit_with_error('Error: could not find ' + filename)

def quit_with_error(message):
    '''
    Displays the given message and ends the program's execution.
    '''
    print('Error:', message, file=sys.stderr)
    quit()

def get_best_k_type_match(assembly, k_refs_fasta, k_refs_dict):
    '''
    Searches for all known K-types in the given assembly and returns the best match.
    Best match is defined as the K-type for which the largest fraction of the K-type has a BLAST
    hit to the assembly.
    '''
    assembly_name = os.path.splitext(assembly)[0]
    for k_ref in k_refs_dict.itervalues():
        k_ref.blast_hits[assembly_name] = []
    if not blast_database_exists(assembly):
        make_blast_database(assembly)
    blast_hits = get_blast_hits(assembly, k_refs_fasta)
    for hit in blast_hits:
        if hit.qseqid not in k_refs_dict:
            quit_with_error('BLAST hit (' + hit.qseqid + ') not found in K-locus references')
        k_refs_dict[hit.qseqid].blast_hits[assembly_name].append(hit)
    best_k_ref = None
    best_fraction_hit = 0.0
    for k_ref in k_refs_dict.itervalues():
        fraction_hit = k_ref.get_fraction_hit_length(assembly_name)
        if fraction_hit > best_fraction_hit:
            best_fraction_hit = fraction_hit
            best_k_ref = k_ref
    return best_k_ref

def blast_database_exists(filename):
    '''
    Returns whether or not a BLAST database exists for the given file.
    '''
    return os.path.isfile(filename + '.nin') and \
           os.path.isfile(filename + '.nhr') and \
           os.path.isfile(filename + '.nsq')

def make_blast_database(filename):
    '''
    Runs makeblastdb on the given assembly.
    '''
    makeblastdb_command = ['makeblastdb', '-dbtype', 'nucl', '-in', filename]
    process = subprocess.Popen(makeblastdb_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, err = process.communicate()
    if err:
        quit_with_error('makeblastdb encountered an error:\n' + err)

def get_blast_hits(database, query):
    '''
    Returns a list BlastHit objects for a search of the given query in the given database.
    '''
    blastn_command = ['blastn', '-task', 'blastn', '-db', database, '-query', query, '-outfmt',
                      '6 qseqid sseqid qstart qend sstart send sseq evalue bitscore length pident']
    process = subprocess.Popen(blastn_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    if err:
        quit_with_error('blastn encountered an error:\n' + err)
    blast_hits = [BlastHit(line) for line in line_iterator(out)]
    return blast_hits

def line_iterator(string_with_line_breaks):
    '''
    Iterates over a string containing line breaks, one line at a time.
    '''
    prev_newline = -1
    while True:
        next_newline = string_with_line_breaks.find('\n', prev_newline + 1)
        if next_newline < 0:
            break
        yield string_with_line_breaks[prev_newline + 1:next_newline]
        prev_newline = next_newline

def load_k_locus_references(fasta):
    '''
    Returns a dictionary of:
      key = K-locus name
      value = KLocusReference object
    '''
    return {seq[0]: KLocusReference(seq[0], seq[1]) for seq in load_fasta(fasta)}

def load_fasta(filename):
    '''
    Returns a list of tuples (name, seq) for the sequences in the fasta.
    '''
    fasta_seqs = []
    fasta_file = open(filename, 'r')
    name = ''
    sequence = ''
    for line in fasta_file:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>': # Header line = start of new contig
            if name:
                fasta_seqs.append((name, sequence))
                name = ''
                sequence = ''
            name = line[1:]
        else:
            sequence += line
    if name:
        fasta_seqs.append((name, sequence))
    return fasta_seqs

def get_assembly_pieces(assembly, k_type):
    '''
    This function uses the BLAST hits in the given K-type to find the corresponding pieces of the
    given assembly.
    It uses the following logic:
      * If the start and end (with a bit of wiggle room) of the K-locus both hit to the same
        contig, and the distance between matches the K-locus size (with a bit of wiggle room),
        then we just return that one piece of that one contig (ideal scenario).
      * If the first case doesn't apply (either because the start and end are on different contigs
        or because the length doesn't match up), then we gather the assembly pieces as follows:
          - Sort the BLAST hits by alignment length from largest to smallest
          - For each BLAST hit, we keep it if it offers new parts of the K-locus. If, on the other
            hand it lies entirely within an existing hit (in K-locus positions), we ignore it.
          - For all of the BLAST hits we've kept, I merge any overlapping ones (in assembly
            positions) and return those.
    The first (better) case will result in a '+' for the confidence call.  The second (worse) case
    will result in a '?'.
    '''
    TO DO
    TO DO
    TO DO
    TO DO
    TO DO
    TO DO
    TO DO
    TO DO
    TO DO
    TO DO


class BlastHit:
    def __init__(self, hit_string):
        parts = hit_string.split('\t')
        self.qseqid = parts[0]
        self.sseqid = parts[1]
        self.qstart = int(parts[2])
        self.qend = int(parts[3])
        self.sstart = int(parts[4])
        self.send = int(parts[5])
        self.sseq = parts[6]
        self.evalue = float(parts[7])
        self.bitscore = float(parts[8])
        self.length = int(parts[9])
        self.pident = float(parts[10])

    def __repr__(self):
        return 'Query: ' + self.qseqid + ' (' + str(self.qstart) + '-' + str(self.qend) + ') ' + \
               'Subject: ' + self.sseqid + ' (' + str(self.sstart) + '-' + str(self.send) + ')'

    def get_assembly_piece(self, assembly):
        '''
        Returns the piece of the assembly which corresponds to this BLAST hit.
        '''
        TO DO
        TO DO
        TO DO
        TO DO
        TO DO
        TO DO
        TO DO



class KLocusReference:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.blast_hits = {} # key = assembly, value = BlastHit

    def __repr__(self):
        return 'K-locus ' + self.name

    def get_total_hit_length(self, assembly_name):
        '''
        Returns the number of bases in this K-locus which are covered by BLAST hits in the given
        assembly.
        '''
        if assembly_name not in self.blast_hits:
            return 0
        hit_ranges = [(x.qstart, x.qend) for x in self.blast_hits[assembly_name]]
        hit_locations = [0] * len(self.seq)
        for hit_range in hit_ranges:
            fixed_range = (min(hit_range) - 1, max(hit_range))
            for i in xrange(fixed_range[0], fixed_range[1]):
                hit_locations[i] = 1
        return sum(hit_locations)

    def get_fraction_hit_length(self, assembly_name):
        '''
        Returns the fraction of this K-locus which is covered by BLAST hits in the given assembly.
        '''
        return self.get_total_hit_length(assembly_name) / len(self.seq)


class AssemblyPiece:
    def __init__(self, assembly, contig, seq):
        self.assembly_name = assembly_name
        self.contig = contig
        self.seq = seq
        TO DO
        TO DO
        TO DO
        TO DO
        TO DO
        TO DO
        TO DO
        TO DO
        TO DO
        TO DO
        TO DO
        TO DO
        TO DO



class IntRange:
    def __init__(self, ranges):
        self.ranges = []
        self.add_ranges(ranges)
        self.simplify()

    def add_range(self, one_range):
        self.add_ranges([one_range])

    def add_ranges(self, ranges):
        self.ranges += ranges
        self.simplify()

    def simplify(self):
        '''
        Collapses overlapping ranges together.
        '''
        TO DO: double check 0-based vs 1-based, inclusive vs exclusive
        TO DO: double check 0-based vs 1-based, inclusive vs exclusive
        TO DO: double check 0-based vs 1-based, inclusive vs exclusive
        TO DO: double check 0-based vs 1-based, inclusive vs exclusive
        TO DO: double check 0-based vs 1-based, inclusive vs exclusive
        TO DO: double check 0-based vs 1-based, inclusive vs exclusive
        TO DO: double check 0-based vs 1-based, inclusive vs exclusive
        TO DO: double check 0-based vs 1-based, inclusive vs exclusive

        fixed_ranges = []
        for int_range in ranges:
            if int_range[1] < int_range[0]:
                fixed_ranges.append((int_range[1], int_range[0]))
            else:
                fixed_ranges.append(int_range)
        starts_ends = [(x[0], 1) for x in fixed_ranges]
        starts_ends += [(x[1], -1) for x in fixed_ranges]
        starts_ends.sort(key=lambda x: x[0])
        current_sum = 0
        cumulative_sum = []
        for start_end in starts_ends:
            current_sum += start_end[1]
            cumulative_sum.append((start_end[0], current_sum))
        prev_depth = 0
        start = 0
        combined = []
        for pos, depth in cumulative_sum:
            if prev_depth == 0:
                start = pos
            elif depth == 0:
                combined.append((start, pos))
            prev_depth = depth
        return combined
        


if __name__ == '__main__':
    main()
