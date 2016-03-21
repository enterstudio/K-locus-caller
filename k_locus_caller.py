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


def main():
    '''
    Script execution starts here.
    '''
    args = get_arguments()
    check_for_blast()
    check_files_exist(args.assembly + [args.k_ref_seqs] + [args.k_ref_genes] + [args.gene_seqs])
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    k_refs = load_k_locus_references(args.k_ref_seqs, args.k_ref_genes) # type: dict[str, KLocusReference]
    for fasta_file in args.assembly:
        assembly = Assembly(fasta_file)
        best_k = get_best_k_type_match(assembly, args.k_ref_seqs, k_refs)
        pieces_matching_k_locus, ideal = get_assembly_pieces(assembly, best_k, args)
        save_assembly_pieces_to_file(pieces_matching_k_locus, args.outdir, best_k.name)
        gene_results = protein_blast(assembly, pieces_matching_k_locus, best_k, args)
        output_to_table(assembly, best_k, pieces_matching_k_locus, ideal, gene_results)

        print() # TEMP
        

    sys.exit(0)



def get_arguments():
    '''
    Specifies the command line arguments required by the script.
    '''
    parser = argparse.ArgumentParser(description='K-locus caller')
    parser.add_argument('-a', '--assembly', nargs='+', type=str, required=True,
                        help='Fasta file(s) for Klebsiella assemblies')
    parser.add_argument('-k', '--k_ref_seqs', type=str, required=True,
                        help='Fasta file with reference K-locus nucleotide sequences')
    parser.add_argument('-g', '--k_ref_genes', type=str, required=True,
                        help='Table file specifying which genes occur in which K-locus')
    parser.add_argument('-s', '--gene_seqs', type=str, required=True,
                        help='Fasta file with protein sequences for each gene')
    parser.add_argument('-o', '--outdir', type=str, required=True,
                        help='Output directory')
    parser.add_argument('--start_end_margin', type=int, required=False, default=10,
                        help='Missing bases at the ends of K-locus allowed in a perfect match.')
    parser.add_argument('--allowed_length_error', type=int, required=False, default=5,
                        help='% error in K-locus length allowed in a perfect match')
    parser.add_argument('--min_gene_cov', type=float, required=False, default=90.0,
                        help='minimum required % coverage for genes')
    parser.add_argument('--min_gene_id', type=float, required=False, default=90.0,
                        help='minimum required % identity for genes')
    return parser.parse_args()

def check_for_blast(): # type: () -> bool
    '''
    Checks to make sure the required BLAST+ tools are available.
    '''
    if not find_program('makeblastdb'):
        quit_with_error('could not find makeblastdb tool (part of BLAST+)')
    if not find_program('blastn'):
        quit_with_error('could not find blastn tool (part of BLAST+)')
    if not find_program('tblastn'):
        quit_with_error('could not find tblastn tool (part of BLAST+)')

def find_program(name): # type: (str) -> bool
    '''
    Checks to see if a program exists.
    '''
    process = subprocess.Popen(['which', name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    return bool(out) and not bool(err)

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

def quit_with_error(message): # type: (str) -> None
    '''
    Displays the given message and ends the program's execution.
    '''
    print('Error:', message, file=sys.stderr)
    sys.exit(1)

def get_best_k_type_match(assembly, k_refs_fasta, k_refs):
    # type: (Assembly, str, dict[str, KLocusReference]) -> KLocusReference
    '''
    Searches for all known K-types in the given assembly and returns the best match.
    Best match is defined as the K-type for which the largest fraction of the K-type has a BLAST
    hit to the assembly.
    '''
    for k_ref in k_refs.itervalues():
        k_ref.clear_hits()
    blast_hits = get_blast_hits(assembly.fasta, k_refs_fasta)
    for hit in blast_hits:
        if hit.qseqid not in k_refs:
            quit_with_error('BLAST hit (' + hit.qseqid + ') not found in K-locus references')
        k_refs[hit.qseqid].add_blast_hit(hit)
    for k_ref in k_refs.itervalues():
        k_ref.sort_hits()
    best_k_ref = None
    best_fraction_hit = 0.0
    for k_ref in k_refs.itervalues():
        fraction_hit = k_ref.get_fraction_hit_length()
        if fraction_hit > best_fraction_hit:
            best_fraction_hit = fraction_hit
            best_k_ref = k_ref
    return best_k_ref

def get_assembly_pieces(assembly, k_type, args):
    '''
    This function uses the BLAST hits in the given K-type to find the corresponding pieces of the
    given assembly.  It returns a list of the pieces and a boolean value to indicate whether an
    ideal match was found.

    It uses the following logic:
      * If the start and end (with a bit of wiggle room) of the K-locus both hit to the same
        contig, and the distance between matches the K-locus size (with a bit of wiggle room),
        then we just return that one piece of that one contig (ideal scenario).
      * If the first case doesn't apply (either because the start and end are on different contigs
        or because the length doesn't match up), then we gather the assembly pieces as follows:
          - For each BLAST hit, we keep it if it offers new parts of the K-locus. If, on the other
            hand it lies entirely within an existing hit (in K-locus positions), we ignore it.
            Since the BLAST hits are sorted longest to shortest, this strategy will prioritise long
            hits over short ones.
          - For all of the BLAST hits we've kept, I merge any overlapping ones (in assembly
            positions) and return those.
    The first (ideal) case will result in a '+' for the confidence call.  The second case will
    result in a '?'.
    '''
    if not k_type.blast_hits:
        return [], False

    # Check for the ideal case.
    earliest_hit = k_type.get_earliest_hit()
    latest_hit = k_type.get_latest_hit()
    start = earliest_hit.qstart
    end = latest_hit.qend
    k_len = k_type.get_length()
    good_start = start <= args.start_end_margin
    good_end = end >= k_len - args.start_end_margin
    same_contig = earliest_hit.sseqid == latest_hit.sseqid
    same_strand = earliest_hit.strand == latest_hit.strand
    proper_length = 100.0 * abs(k_len - (end - start)) / k_len <= args.allowed_length_error
    if good_start and good_end and same_contig and same_strand and proper_length:
        if earliest_hit.strand == '+':
            contig_start, contig_end = earliest_hit.sstart, latest_hit.send
        else:
            contig_start, contig_end = latest_hit.sstart, earliest_hit.send
        one_piece = AssemblyPiece(assembly, earliest_hit.sseqid, contig_start, contig_end,
                                  earliest_hit.strand)
        return [one_piece], True

    # If we got here, then it's the non-ideal case.  Keep BLAST hits that give us new parts of the
    # K-locus and return the corresponding pieces of the assembly.
    used_hits = []
    k_range_so_far = IntRange()
    for hit in k_type.blast_hits:
        hit_range = hit.get_query_range()
        if not k_range_so_far.contains(hit_range):
            k_range_so_far.merge_in_range(hit_range)
            used_hits.append(hit)
    assembly_pieces = [x.get_assembly_piece(assembly) for x in used_hits]
    return merge_assembly_pieces(assembly_pieces), False

def protein_blast(assembly, pieces_matching_k_locus, k_locus_ref, args):
    '''
    Conducts a BLAST search of all known K-locus proteins.  Returns the BLAST hits for 
    '''
    blast_hits = get_blast_hits(assembly.fasta, args.gene_seqs, genes=True)
    filtered_hits = [x for x in blast_hits if x.query_cov >= args.min_gene_cov and x.pident >= args.min_gene_id]
    expected_gene_hits = []
    missing_expected_genes = []
    for expected_gene in k_locus_ref.gene_names:
        best_hit = get_best_hit_for_query(filtered_hits, expected_gene)
        if not best_hit:
            missing_expected_genes.append(expected_gene)
        else:
            expected_gene_hits.append(best_hit)
            filtered_hits = [x for x in filtered_hits if x is not best_hit]
            filtered_hits = cull_conflicting_hits(best_hit, filtered_hits)
    other_hits = cull_all_conflicting_hits(filtered_hits)
    expected_genes_inside_locus = [x for x in expected_gene_hits if x.in_assembly_pieces(pieces_matching_k_locus)]
    expected_genes_outside_locus = [x for x in expected_gene_hits if not x.in_assembly_pieces(pieces_matching_k_locus)]
    other_genes_inside_locus = [x for x in other_hits if x.in_assembly_pieces(pieces_matching_k_locus)]
    other_genes_outside_locus = [x for x in other_hits if not x.in_assembly_pieces(pieces_matching_k_locus)]
    return missing_expected_genes, expected_genes_inside_locus, expected_genes_outside_locus, \
           other_genes_inside_locus, other_genes_outside_locus

def output_to_table(assembly, best_k, pieces, ideal, gene_results):

    missing, exp_in, exp_out, other_in, other_out = gene_results

    # TEMP
    print('\n')
    print('Assembly: ' + str(assembly) + ', best K-type match: ' + best_k.name + ', ideal: ' + str(ideal))
    for piece in pieces:
        print(str(piece) + ': ' + piece.get_sequence_short())

    # TEMP
    print('\n')
    print('missing_expected_genes')
    for gene in missing:
        print(gene)
    print('\n')
    print('expected_genes_inside_locus')
    for hit in exp_in:
        print(hit)
    print('\n')
    print('expected_genes_outside_locus')
    for hit in exp_out:
        print(hit)
    print('\n')
    print('other_genes_inside_locus')
    for hit in other_in:
        print(hit)
    print('\n')
    print('other_genes_outside_locus')
    for hit in other_out:
        print(hit)





def get_blast_hits(database, query, genes=False):
    '''
    Returns a list BlastHit objects for a search of the given query in the given database.
    '''
    if genes:
        command = ['tblastn']
    else:
        command = ['blastn']
    command += ['-db', database, '-query', query, '-outfmt',
                '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    if err:
        quit_with_error('blastn encountered an error:\n' + err)
    if genes:
        blast_hits = [GeneBlastHit(line) for line in line_iterator(out)]
    else:
        blast_hits = [BlastHit(line) for line in line_iterator(out)]
    return blast_hits

def get_best_hit_for_query(blast_hits, query_name):
    '''
    Given a list of BlastHits, this function returns the best hit for the given query, based on
    identity.  It returns None if no BLAST hits match that query.
    '''
    matching_hits = [x for x in blast_hits if x.qseqid == query_name]
    if matching_hits:
        return sorted(matching_hits, key=lambda x: x.pident, reverse=True)[0]
    else:
        return None

def cull_conflicting_hits(hit_to_keep, blast_hits):
    '''
    This function returns a (potentially) reduced set of BLAST hits which excludes BLAST hits that
    overlap with the hit to keep (same part of assembly) and are in the same cluster as the hit to
    keep.
    '''
    return [x for x in blast_hits if not x.conflicts(hit_to_keep)]


def cull_all_conflicting_hits(blast_hits):
    '''
    This function returns a (potentially) reduced set of BLAST hits where none of the remaining
    hits conflict.
    '''
    blast_hits.sort(key=lambda x: x.pident, reverse=True)
    kept_hits = []
    while blast_hits:
        kept_hits.append(blast_hits.pop(0))
        blast_hits = cull_conflicting_hits(kept_hits[-1], blast_hits)
    return kept_hits

def merge_assembly_pieces(pieces):
    '''
    Takes a list of AssemblyPiece objects and returns another list of AssemblyPiece objects where
    the overlapping pieces have been merged.
    '''
    merged_pieces = []
    while pieces:
        merged_piece = pieces[0]
        unmerged = []
        for other_piece in pieces[1:]:
            combined = merged_piece.combine(other_piece)
            if not combined:
                unmerged.append(other_piece)
            else:
                merged_piece = combined
        merged_pieces.append(merged_piece)
        pieces = unmerged
    return merged_pieces

def reverse_complement(seq):
    '''
    Given a DNA sequences, this function returns the reverse complement sequence.
    '''
    rev_comp = ''
    for i in reversed(range(len(seq))):
        rev_comp += complement_base(seq[i])
    return rev_comp

def complement_base(base):
    '''
    Given a DNA base, this returns the complement.
    '''
    forward = 'ATGCatgcRYSWKMryswkmBDHVbdhvNn.-?'
    reverse = 'TACGtacgYRSWMKyrswmkVHDBvhdbNn.-?N'
    return reverse[forward.find(base)]

def save_assembly_pieces_to_file(assembly_pieces, outdir, k_locus_name):
    '''
    Creates a single FASTA file for all of the assembly pieces.
    Assumes all assembly pieces are from the same assembly.
    '''
    if not assembly_pieces:
        return
    assembly_name = assembly_pieces[0].assembly.name
    fasta_file_name = os.path.join(outdir, assembly_name + '_' + k_locus_name + '.fasta')
    fasta_file = open(fasta_file_name, 'w')
    for piece in assembly_pieces:
        fasta_file.write('>' + piece.get_header() + '\n')
        fasta_file.write(add_line_breaks_to_sequence(piece.get_sequence(), 60))

def add_line_breaks_to_sequence(sequence, length):
    '''
    Wraps sequences to the defined length.  All resulting sequences end in a line break.
    '''
    seq_with_breaks = ''
    while len(sequence) > length:
        seq_with_breaks += sequence[:length] + '\n'
        sequence = sequence[length:]
    if sequence:
        seq_with_breaks += sequence
        seq_with_breaks += '\n'
    return seq_with_breaks

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

def load_k_locus_references(fasta, table): # type: (str, str) -> dict[str, KLocusReference]
    '''
    Returns a dictionary of:
      key = K-locus name
      value = KLocusReference object
    '''
    return {seq[0]: KLocusReference(seq[0], seq[1], table) for seq in load_fasta(fasta)}

def load_fasta(filename): # type: (str) -> list[tuple[str, str]]
    '''
    Returns the names and sequences for the given fasta file.
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



class BlastHit(object):
    '''
    Stores the BLAST hit output mostly verbatim.  However, it does convert the BLAST ranges
    (1-based, inclusive end) to Python ranges (0-based, exclusive end).
    '''
    def __init__(self, hit_string):
        parts = hit_string.split('\t')
        self.qseqid = parts[0]
        self.sseqid = parts[1]
        self.qstart = int(parts[2]) - 1
        self.qend = int(parts[3])
        self.sstart = int(parts[4])
        self.send = int(parts[5])
        if self.sstart <= self.send:
            self.strand = '+'
        else:
            self.sstart, self.send = self.send, self.sstart
            self.strand = '-'
        self.sstart -= 1
        self.evalue = float(parts[6])
        self.bitscore = float(parts[7])
        self.length = int(parts[8])
        self.pident = float(parts[9])
        self.query_cov = 100.0 * len(parts[11]) / float(parts[10])

    def __repr__(self):
        return 'Query: ' + self.qseqid + ' (' + str(self.qstart) + '-' + str(self.qend) + ') ' + \
               'Subject: ' + self.sseqid + ' (' + str(self.sstart) + '-' + str(self.send) + \
               ', ' + self.strand + ' strand) ' + \
               'Cov: ' + str(self.query_cov) + '%, ID: ' + str(self.pident) + '%'

    def get_assembly_piece(self, assembly):
        '''
        Returns the piece of the assembly which corresponds to this BLAST hit.
        '''
        return AssemblyPiece(assembly, self.sseqid, self.sstart, self.send, self.strand)

    def get_query_range(self):
        '''
        Produces an IntRange object for the hit query.
        '''
        return IntRange([(self.qstart, self.qend)])

    def in_assembly_pieces(self, assembly_pieces):
        '''
        Returns True if the hit is in (or at least overlaps with) any of the given assembly pieces.
        '''
        for piece in assembly_pieces:
            if piece.overlaps(self.sseqid, self.sstart, self.send):
                return True
        return False



class GeneBlastHit(BlastHit):
    '''
    This class adds a few gene-specific things to the BlastHit class.
    '''
    def __init__(self, hit_string):
        BlastHit.__init__(self, hit_string)
        gene_name_parts = self.qseqid.split('__')
        self.cluster = int(gene_name_parts[0])
        self.gene_name = gene_name_parts[2]

    def conflicts(self, other):
        '''
        Returns whether or not this hit conflicts with the other hit.  Conflicting hits overlap on
        the assembly and are in the same cluster.
        A hit is not considered to conflict with itself.
        '''
        if self is other:
            return False
        if self.sseqid != other.sseqid:
            return False
        if self.strand != other.strand:
            return False
        if self.cluster != other.cluster:
            return False
        return self.sstart <= other.send and other.sstart <= self.send



class KLocusReference(object):
    def __init__(self, name, seq, table):
        self.name = name
        self.seq = seq
        self.gene_names = []
        self.blast_hits = []
        self.hit_ranges = IntRange()
        self.load_genes(table)

    def __repr__(self):
        return 'K-locus ' + self.name

    def get_length(self):
        '''
        Returns the K-locus sequence length.
        '''
        return len(self.seq)

    def add_blast_hit(self, hit):
        '''
        Adds a BLAST hit and updates the hit ranges.
        '''
        self.blast_hits.append(hit)
        self.hit_ranges.add_range(hit.qstart, hit.qend)

    def clear_hits(self):
        '''
        Clears all BLAST hits and the corresponding hit ranges.
        '''
        self.blast_hits = []
        self.hit_ranges = IntRange()

    def sort_hits(self):
        '''
        Sorts the BLAST hits from longest to shortest.
        '''
        self.blast_hits.sort(key=lambda x: x.length, reverse=True)

    def get_fraction_hit_length(self):
        '''
        Returns the fraction of this K-locus which is covered by BLAST hits in the given assembly.
        '''
        return self.hit_ranges.get_total_length() / len(self.seq)

    def get_earliest_hit(self):
        '''
        Returns the BLAST hit that is earliest in the K-locus sequence (lowest start coordinate).
        '''
        earliest_hit = self.blast_hits[0]
        for hit in self.blast_hits[1:]:
            if hit.qstart < earliest_hit.qstart:
                earliest_hit = hit
        return earliest_hit

    def get_latest_hit(self):
        '''
        Returns the BLAST hit that is latest in the K-locus sequence (highest end coordinate).
        '''
        latest_hit = self.blast_hits[0]
        for hit in self.blast_hits[1:]:
            if hit.qend > latest_hit.qend:
                latest_hit = hit
        return latest_hit

    def load_genes(self, table_filename): # type: (str) -> None
        '''
        Reads the table file which gives the genes in each K-locus and remembers which genes belong
        in this K-locus.
        '''
        table_file = open(table_filename, 'r')
        for line in table_file:
            line_parts = line.strip().split('\t')
            if len(line_parts) != 2:
                continue
            k_locus_name = line_parts[0]
            gene_name = line_parts[1]
            if k_locus_name == self.name:
                self.gene_names.append(gene_name)



class Assembly(object):
    def __init__(self, fasta_file):
        '''
        Loads in an assembly and builds a BLAST database for it (if necessary).
        '''
        self.fasta = fasta_file
        self.name = os.path.splitext(fasta_file)[0]
        self.contigs = {x[0]: x[1] for x in load_fasta(fasta_file)} # key = name, value = sequence
        if not self.blast_database_exists():
            self.make_blast_database()

    def __repr__(self):
        return self.name

    def blast_database_exists(self):
        '''
        Returns whether or not a BLAST database already exists for this assembly.
        '''
        return os.path.isfile(self.fasta + '.nin') and \
               os.path.isfile(self.fasta + '.nhr') and \
               os.path.isfile(self.fasta + '.nsq')

    def make_blast_database(self):
        '''
        Runs makeblastdb on the assembly.
        '''
        command = ['makeblastdb', '-dbtype', 'nucl', '-in', self.fasta]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, err = process.communicate()
        if err:
            quit_with_error('makeblastdb encountered an error:\n' + err)




class AssemblyPiece(object):
    '''
    This class describes a piece of an assembly: which contig the piece is on and what the range is.
    '''
    def __init__(self, assembly, contig_name, contig_start, contig_end, strand):
        self.assembly = assembly
        self.contig_name = contig_name
        self.start = contig_start
        self.end = contig_end
        self.strand = strand

    def __repr__(self):
        return self.assembly.name + '_' + self.get_header()

    def get_header(self):
        '''
        Returns a descriptive string for the FASTA header when saving this piece to file.
        '''
        contig_name_parts = self.contig_name.split('_')
        contig_number = contig_name_parts[1]
        return 'NODE_' + contig_number + '_' + str(self.start + 1) + '_to_' + str(self.end) + \
               '_' + self.strand + '_strand'

    def get_sequence(self):
        '''
        Returns the DNA sequence for this piece of the assembly.
        '''
        seq = self.assembly.contigs[self.contig_name][self.start:self.end]
        if self.strand == '+':
            return seq
        else:
            return reverse_complement(seq)

    def get_sequence_short(self):
        '''
        Returns a shortened format of the sequence
        '''
        seq = self.get_sequence()
        length = len(seq)
        if len(seq) > 9:
            seq = seq[0:6] + '...' + seq[-6:]
        return seq + ' (' + str(length) + ' bp)'

    def combine(self, other):
        '''
        If this assembly piece and the other can be combined, this function returns the combined
        piece.  If they can't, it returns None.
        '''
        if self.contig_name != other.contig_name or self.strand != other.strand:
            return False
        combined = IntRange([(self.start, self.end)])
        combined.add_range(other.start, other.end)
        if len(combined.ranges) == 1:
            new_start, new_end = combined.ranges[0]
            return AssemblyPiece(self.assembly, self.contig_name, new_start, new_end, self.strand)
        else:
            return None

    def overlaps(self, contig_name, start, end):
        '''
        Returns whether this assembly piece overlaps with the given parameters.
        '''
        return self.contig_name == contig_name and self.start <= end and start <= self.end



class IntRange(object):
    '''
    This class contains one or more integer ranges.  Overlapping ranges will be merged together.
    It stores its ranges in a Python-like fashion where the last value in each range is
    exclusive.
    '''
    def __init__(self, ranges = []):
        self.ranges = []
        self.add_ranges(ranges)
        self.simplify()

    def __repr__(self):
        return str(self.ranges)

    def add_range(self, start, end):
        '''
        Adds a single range.
        '''
        self.add_ranges([(start, end)])

    def add_ranges(self, ranges):
        '''
        Adds multiple ranges (list of tuples)
        '''
        self.ranges += ranges
        self.simplify()

    def merge_in_range(self, other):
        '''
        Merges the other IntRange object into this one.
        '''
        self.add_ranges(other.ranges)

    def get_total_length(self):
        '''
        Returns the number of integers in the ranges.
        '''
        return sum([x[1] - x[0] for x in self.ranges])

    def simplify(self):
        '''
        Collapses overlapping ranges together.
        '''
        fixed_ranges = []
        for int_range in self.ranges:
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
        self.ranges = combined

    def contains(self, other):
        '''
        Returns True if the other IntRange is entirely contained within this IntRange.
        '''
        for other_range in other.ranges:
            other_start, other_end = other_range
            contained = False
            for this_range in self.ranges:
                this_start, this_end = this_range
                if other_start >= this_start and other_end <= this_end:
                    contained = True
                    break
            if not contained:
                return False
        return True



if __name__ == '__main__':
    main()
