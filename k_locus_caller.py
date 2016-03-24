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
    make_paths_absolute(args)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    check_inputs(args)
    k_refs = load_k_locus_references(args.k_ref_seqs, args.k_ref_genes) # type: dict[str, KLocus]
    table_file = create_table_file(args.outdir)
    for fasta_file in args.assembly:
        assembly = Assembly(fasta_file)
        best_k = get_best_k_type_match(assembly, args.k_ref_seqs, k_refs)
        find_assembly_pieces(assembly, best_k, args)
        protein_blast(assembly, best_k, args)
        output(table_file, assembly, best_k, args)
        save_assembly_pieces_to_file(best_k, assembly, args.outdir)
    sys.exit(0)



def get_arguments():
    '''
    Specifies the command line arguments required by the script.
    '''
    parser = argparse.ArgumentParser(description='K-locus caller',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    parser.add_argument('--min_gene_cov', type=float, required=False, default=90.0,
                        help='minimum required %% coverage for genes')
    parser.add_argument('--min_gene_id', type=float, required=False, default=70.0,
                        help='minimum required %% identity for genes')
    parser.add_argument('--low_gene_id', type=float, required=False, default=95.0,
                        help='genes with a %% identity below this value will be flagged as low '
                             'identity')
    parser.add_argument('--min_assembly_piece', type=int, required=False, default=100,
                        help='minimum K-locus matching assembly piece to return')
    parser.add_argument('--gap_fill_size', type=int, required=False, default=100,
                        help='when separate parts of the assembly are found within this distance, '
                             'they will be merged')
    parser.add_argument('--verbose', action='store_true',
                        help='Display detailed information about each assembly in stdout')

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

def make_paths_absolute(args):
    '''
    Changes the paths given by the user to absolute paths, which are easier to work with later.
    '''
    args.assembly = [os.path.abspath(x) for x in args.assembly]
    args.k_ref_seqs = os.path.abspath(args.k_ref_seqs)
    args.k_ref_genes = os.path.abspath(args.k_ref_genes)
    args.gene_seqs = os.path.abspath(args.gene_seqs)
    args.outdir = os.path.abspath(args.outdir)

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

def check_inputs(args):
    '''
    Makes sure that:
      1) the gene names are in the correct format
      2) the gene names in the table are present in the FASTA file
      3) each K locus sequence has at least one gene in the table
      4) the K loci in the table are present in the FASTA file
    If any of these are not true, it quits with an error message.
    '''
    k_names_from_fasta = set([x[0] for x in load_fasta(args.k_ref_seqs)])
    if not k_names_from_fasta:
        quit_with_error('No K-locus reference sequences found in ' +
                        os.path.basename(args.k_ref_seqs))
    gene_names_from_fasta = set([x[0] for x in load_fasta(args.gene_seqs)])
    if not gene_names_from_fasta:
        quit_with_error('No gene sequences found in ' + os.path.basename(args.gene_seqs))
    k_names_from_table = set()
    gene_names_from_table = set()
    table_file = open(args.k_ref_genes, 'r')
    for line in table_file:
        line_parts = line.strip().split('\t')
        if len(line_parts) == 2:
            k_names_from_table.add(line_parts[0])
            gene_names_from_table.add(line_parts[1])
    if not k_names_from_table or not gene_names_from_table:
        quit_with_error('No K-locus or gene names found in ' + 
                        os.path.basename(args.k_ref_genes))
    missing_genes = gene_names_from_table.difference(gene_names_from_fasta)
    if missing_genes:
        quit_with_error('These genes are present in the table but not in the gene FASTA '
                        'file: ' + ', '.join(list(missing_genes)))
    missing_k_loci = k_names_from_table.difference(k_names_from_fasta)
    if missing_k_loci:
        quit_with_error('These K loci are present in the table but not in the K locus '
                        'FASTA file: ' + ', '.join(list(missing_k_loci)))
    missing_k_loci = k_names_from_fasta.difference(k_names_from_table)
    if missing_k_loci:
        quit_with_error('These K loci are present in the FASTA file but not in the '
                        'table: ' + ', '.join(list(missing_k_loci)))
    for gene_name in gene_names_from_table:
        gene_name_parts = gene_name.split('__')
        if len(gene_name_parts) < 3:
            quit_with_error('This gene name is not in the correct format: ' + gene_name)
        try:
            int(gene_name_parts[0])
        except ValueError:
            quit_with_error('This gene\'s cluster number is not a number: ' + gene_name)
        if not gene_name_parts[2]:
            quit_with_error('This gene is missing an allele name: ' + gene_name)

def get_best_k_type_match(assembly, k_refs_fasta, k_refs):
    # type: (Assembly, str, dict[str, KLocus]) -> KLocus
    '''
    Searches for all known K-types in the given assembly and returns the best match.
    Best match is defined as the K-type for which the largest fraction of the K-type has a BLAST
    hit to the assembly.
    '''
    for k_ref in k_refs.itervalues():
        k_ref.clear()
    blast_hits = get_blast_hits(assembly.fasta, k_refs_fasta)
    for hit in blast_hits:
        if hit.qseqid not in k_refs:
            quit_with_error('BLAST hit (' + hit.qseqid + ') not found in K-locus references')
        k_refs[hit.qseqid].add_blast_hit(hit)
    best_k_ref = None
    best_cov = 0.0
    for k_ref in k_refs.itervalues():
        cov = k_ref.get_coverage()
        if cov > best_cov:
            best_cov = cov
            best_k_ref = k_ref
    best_k_ref.clean_up_blast_hits()
    return best_k_ref

def find_assembly_pieces(assembly, k_locus, args):
    '''
    This function uses the BLAST hits in the given K-type to find the corresponding pieces of the
    given assembly.  It saves its results in the KLocus
    '''
    if not k_locus.blast_hits:
        return
    assembly_pieces = [x.get_assembly_piece(assembly) for x in k_locus.blast_hits]
    merged_pieces = merge_assembly_pieces(assembly_pieces)
    length_filtered_pieces = [x for x in merged_pieces if x.get_length() >= args.min_assembly_piece]
    if not length_filtered_pieces:
        return
    k_locus.assembly_pieces = fill_assembly_piece_gaps(length_filtered_pieces, args.gap_fill_size)

    # Now check to see if the biggest assembly piece seems to capture the whole locus.  If so, this
    # is an ideal match.
    biggest_piece = sorted(k_locus.assembly_pieces, key=lambda x: x.get_length(), reverse=True)[0]
    start = biggest_piece.earliest_hit_coordinate()
    end = biggest_piece.latest_hit_coordinate()
    if good_start_and_end(start, end, k_locus.get_length(), args.start_end_margin):
        k_locus.assembly_pieces = [biggest_piece]

    # If it isn't the ideal case, we still want to check if the start and end of the K-locus were
    # found in the same contig.  If so, fill all gaps in between so we include the entire
    # intervening sequence.
    else:
        earliest, latest, same_contig_and_strand = k_locus.get_earliest_and_latest_pieces()
        start = earliest.earliest_hit_coordinate()
        end = latest.latest_hit_coordinate()
        if good_start_and_end(start, end, k_locus.get_length(), args.start_end_margin) and \
        same_contig_and_strand:
            gap_filling_piece = AssemblyPiece(assembly, earliest.contig_name, start, end,
                                              earliest.strand)
            k_locus.assembly_pieces = merge_assembly_pieces(k_locus.assembly_pieces + \
                                                            [gap_filling_piece])
    k_locus.identity = get_mean_identity(k_locus.assembly_pieces)

def protein_blast(assembly, k_locus, args):
    '''
    Conducts a BLAST search of all known K-locus proteins.  Stores the results in the KLocus
    object.
    '''
    hits = get_blast_hits(assembly.fasta, args.gene_seqs, genes=True)
    hits = [x for x in hits if x.query_cov >= args.min_gene_cov and x.pident >= args.min_gene_id]
    expected_hits = []
    for expected_gene in k_locus.gene_names:
        best_hit = get_best_hit_for_query(hits, expected_gene)
        if not best_hit:
            k_locus.missing_expected_genes.append(expected_gene)
        else:
            best_hit.over_identity_threshold = best_hit.pident >= args.low_gene_id
            expected_hits.append(best_hit)
            hits = [x for x in hits if x is not best_hit]
            hits = cull_conflicting_hits(best_hit, hits)
    other_hits = cull_all_conflicting_hits(hits)
    k_locus.expected_hits_inside_locus = [x for x in expected_hits if x.in_assembly_pieces(k_locus.assembly_pieces)]
    k_locus.expected_hits_outside_locus = [x for x in expected_hits if not x.in_assembly_pieces(k_locus.assembly_pieces)]
    k_locus.other_hits_inside_locus = [x for x in other_hits if x.in_assembly_pieces(k_locus.assembly_pieces)]
    k_locus.other_hits_outside_locus = [x for x in other_hits if not x.in_assembly_pieces(k_locus.assembly_pieces)]

def create_table_file(outdir):
    '''
    Creates the table file and writes a header line.
    '''
    table_path = os.path.join(outdir, 'k-locus_results.txt')
    table = open(table_path, 'w')
    headers = []
    headers.append('Assembly')
    headers.append('Best matching K-locus')
    headers.append('K-locus match problems')
    headers.append('K-locus match coverage')
    headers.append('K-locus match identity')
    headers.append('K-locus match length discrepancy')
    headers.append('Expected genes found in K-locus')
    headers.append('Expected genes found in K-locus, details')
    headers.append('Expected genes not found in K-locus')
    headers.append('Other genes found in K-locus')
    headers.append('Other genes found in K-locus, details')
    headers.append('Expected genes found outside K-locus')
    headers.append('Expected genes found outside K-locus, details')
    headers.append('Other genes found outside K-locus')
    headers.append('Other genes found outside K-locus, details')
    table.write('\t'.join(headers))
    table.write('\n')
    table.flush()
    return table

def output(table, assembly, k_locus, args):
    '''
    Writes a line to the output table describing all that we've learned about the given K-locus and
    writes to stdout as well.
    '''
    uncertainty_chars = k_locus.get_match_uncertainty_chars()
    expected_genes_str = str(len(k_locus.expected_hits_inside_locus)) + ' / ' + \
                         str(len(k_locus.gene_names))
    missing_genes_str = str(len(k_locus.missing_expected_genes)) + ' / ' + \
                        str(len(k_locus.gene_names))
    coverage_str = '%.2f' % k_locus.get_coverage() + '%'
    identity_str = '%.2f' % k_locus.identity + '%'

    line = []
    line.append(assembly.name)
    line.append(k_locus.name)
    line.append(uncertainty_chars)
    line.append(coverage_str)
    line.append(identity_str)
    line.append(k_locus.get_length_discrepancy_string())
    line.append(expected_genes_str)
    line.append(get_gene_info_string(k_locus.expected_hits_inside_locus))
    line.append(k_locus.get_missing_gene_string())
    line.append(str(len(k_locus.other_hits_inside_locus)))
    line.append(get_gene_info_string(k_locus.other_hits_inside_locus))
    line.append(str(len(k_locus.expected_hits_outside_locus)))
    line.append(get_gene_info_string(k_locus.expected_hits_outside_locus))
    line.append(str(len(k_locus.other_hits_outside_locus)))
    line.append(get_gene_info_string(k_locus.other_hits_outside_locus))

    table.write('\t'.join(line))
    table.write('\n')
    table.flush()

    if not args.verbose:
        print(assembly.name + ': ' + k_locus.name + uncertainty_chars)
    if args.verbose:
        print('Assembly: ' + assembly.name)
        print('    Best K-type match: ' + k_locus.name)
        print('    Uncertainties: ' + uncertainty_chars)
        print('    Coverage: ' + coverage_str)
        print('    Identity: ' + identity_str)
        print('    Length discrepancy: ' + k_locus.get_length_discrepancy_string())
        print('    K-locus assembly pieces:')
        for piece in k_locus.assembly_pieces:
            print('        ' + piece.get_header() + ', ' + piece.get_sequence_short())
        print('    Expected genes found in K-locus: ' + expected_genes_str)
        for hit in k_locus.expected_hits_inside_locus:
            print('        ' + str(hit))
        print('    Expected genes missing in K-locus: ' + missing_genes_str)
        for gene in k_locus.missing_expected_genes:
            print('        ' + str(gene))
        print('    Other genes found in K-locus: ' + str(len(k_locus.other_hits_inside_locus)))
        for hit in k_locus.other_hits_inside_locus:
            print('        ' + str(hit))
        print('    Expected genes found outside K-locus: ' + \
              str(len(k_locus.expected_hits_outside_locus)))
        for hit in k_locus.expected_hits_outside_locus:
            print('        ' + str(hit))
        print('    Other genes found outside K-locus: ' + \
              str(len(k_locus.other_hits_outside_locus)))
        for hit in k_locus.other_hits_outside_locus:
            print('        ' + str(hit))
        print()

def get_blast_hits(database, query, genes=False):
    '''
    Returns a list BlastHit objects for a search of the given query in the given database.
    '''
    if genes:
        command = ['tblastn']
    else:
        command = ['blastn', '-task', 'blastn']
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
    bit score.  It returns None if no BLAST hits match that query.
    '''
    matching_hits = [x for x in blast_hits if x.qseqid == query_name]
    if matching_hits:
        return sorted(matching_hits, key=lambda x: x.bitscore, reverse=True)[0]
    else:
        return None

def cull_conflicting_hits(hit_to_keep, blast_hits):
    '''
    This function returns a (potentially) reduced set of BLAST hits which excludes BLAST hits that
    overlap with the hit to keep (same part of assembly) and are in the same cluster as the hit to
    keep.
    If any of the other blast hits have the exact same target range as the hit to keep, they'll be
    culled, regardless of their cluster.
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
    while True:
        merged_pieces = []
        merge_count = 0
        while pieces:
            merged_piece = pieces[0]
            unmerged = []
            for other_piece in pieces[1:]:
                combined = merged_piece.combine(other_piece)
                if not combined:
                    unmerged.append(other_piece)
                else:
                    merged_piece = combined
                    merge_count += 1
            merged_pieces.append(merged_piece)
            pieces = unmerged
        if merge_count == 0:
            break
        else:
            pieces = merged_pieces
    return merged_pieces

def fill_assembly_piece_gaps(pieces, max_gap_fill_size):
    '''
    This function takes a list of assembly pieces, and if any of them are close enough to each
    other, the gap will be merged in.
    It assumes that all given pieces are from the same assembly.
    '''
    pieces_by_contig_and_strand = {}
    fixed_pieces = []
    for piece in pieces:
        contig = piece.contig_name
        strand = piece.strand
        if (contig, strand) not in pieces_by_contig_and_strand:
            pieces_by_contig_and_strand[(contig, strand)] = []
        pieces_by_contig_and_strand[(contig, strand)].append(piece)
    for (contig, strand), pieces_in_contig_and_strand in pieces_by_contig_and_strand.iteritems():
        gap_filling_pieces = []
        sorted_pieces = sorted(pieces_in_contig_and_strand, key=lambda x: x.start)
        max_end = sorted_pieces[0].end
        gaps = []
        for piece in sorted_pieces[1:]:
            if piece.start > max_end and piece.start - max_end <= max_gap_fill_size:
                gaps.append((max_end, piece.start))
            max_end = max(max_end, piece.end)
        assembly = sorted_pieces[0].assembly
        for gap in gaps:
            gap_filling_pieces.append(AssemblyPiece(assembly, contig, gap[0], gap[1], strand))
        before_merge = pieces_in_contig_and_strand + gap_filling_pieces
        filled_pieces = merge_assembly_pieces(before_merge)
        fixed_pieces += filled_pieces
    return fixed_pieces

def get_mean_identity(pieces):
    '''
    Returns the mean identity (weighted by sequence length) for a list of assembly pieces.
    '''
    identity_sum = 0.0
    length_sum = 0
    for piece in pieces:
        for hit in piece.blast_hits:
            length_sum += hit.length
            identity_sum += hit.length * hit.pident
    if identity_sum == 0.0:
        return 0.0
    else:
        return identity_sum / length_sum

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

def save_assembly_pieces_to_file(k_locus, assembly, outdir):
    '''
    Creates a single FASTA file for all of the assembly pieces.
    Assumes all assembly pieces are from the same assembly.
    '''
    if not k_locus.assembly_pieces:
        return
    assembly_and_locus = assembly.name + '_' + k_locus.name
    fasta_file_name = os.path.join(outdir, assembly_and_locus + '.fasta')
    fasta_file = open(fasta_file_name, 'w')
    uncertainties = k_locus.get_match_uncertainty_words()
    for piece in k_locus.assembly_pieces:
        fasta_file.write('>' + assembly_and_locus + '_' + piece.get_header() + ' ' + \
                         uncertainties + '\n')
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

def load_k_locus_references(fasta, table): # type: (str, str) -> dict[str, KLocus]
    '''
    Returns a dictionary of:
      key = K-locus name
      value = KLocus object
    '''
    return {seq[0]: KLocus(seq[0], seq[1], table) for seq in load_fasta(fasta)}

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
                fasta_seqs.append((name.split()[0], sequence))
                name = ''
                sequence = ''
            name = line[1:]
        else:
            sequence += line
    if name:
        fasta_seqs.append((name.split()[0], sequence))
    return fasta_seqs

def good_start_and_end(start, end, k_length, allowed_margin):
    '''
    Checks whether the given start and end coordinates are within the accepted margin of error.
    '''
    good_start = start <= allowed_margin
    good_end = end >= k_length - allowed_margin
    start_before_end = start < end
    return good_start and good_end and start_before_end

def get_gene_info_string(gene_hit_list):
    '''
    Returns a single comma-delimited string summarising the gene hits in the given list.
    '''
    return ';'.join([x.allele_name + ',' + str(x.pident) + '%' for x in gene_hit_list])


def is_contig_name_spades_format(contig_name):
    '''
    Returns whether or not the contig name appears to be in the SPAdes/Velvet format.
    Example: NODE_5_length_150905_cov_4.42519
    '''
    contig_name_parts = contig_name.split('_')
    return len(contig_name_parts) > 5 and contig_name_parts[0] == 'NODE' and \
           contig_name_parts[2] == 'length' and contig_name_parts[4] == 'cov'


def get_nice_contig_name(contig_name):
    '''
    For a contig with a SPAdes/Velvet format, this function returns a simplified string that is
    just NODE_XX where XX is the contig number.
    For any other format, this function trims off everything following the first whitespace.
    '''
    if is_contig_name_spades_format(contig_name):
        return 'NODE_' + contig_name.split('_')[1]
    else:
        return contig_name.split()[0]




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
        return self.qseqid + ' (' + str(self.qstart) + '-' + str(self.qend) + ')   ' + \
               'Contig: ' + self.sseqid + ' (' + str(self.sstart) + '-' + str(self.send) + \
               ', ' + self.strand + ' strand)   ' + \
               'Cov: ' + '%.2f' % self.query_cov + '%, ID: ' + '%.2f' % self.pident + '%'

    def get_assembly_piece(self, assembly):
        '''
        Returns the piece of the assembly which corresponds to this BLAST hit.
        '''
        return AssemblyPiece(assembly, self.sseqid, self.sstart, self.send, self.strand, [self])

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
        self.allele_name = gene_name_parts[2]
        self.over_identity_threshold = False

    def conflicts(self, other):
        '''
        Returns whether or not this hit conflicts with the other hit.
        If the hits are in the same cluster, then any overlap between them counts as a conflict.
        If the hits are in different clusters, then it takes 90% overlap to conflict.
        A hit is not considered to conflict with itself.
        '''
        if self is other:
            return False
        if self.sseqid != other.sseqid:
            return False

        max_start = max(self.sstart, other.sstart)
        min_end = min(self.send, other.send)
        if max_start < min_end:
            overlap = min_end - max_start
        else:
            overlap = 0

        if self.cluster != other.cluster:
            min_length = min(self.send - self.sstart, other.send - other.sstart)
            frac_overlap = overlap / min_length
            return frac_overlap > 0.9
        else:
            return overlap > 0




class KLocus(object):
    def __init__(self, name, seq, table):
        self.name = name
        self.seq = seq
        self.gene_names = []
        self.load_genes(table)
        self.blast_hits = []
        self.hit_ranges = IntRange()
        self.assembly_pieces = []
        self.identity = 0.0
        self.expected_hits_inside_locus = []
        self.missing_expected_genes = []
        self.expected_hits_outside_locus = []
        self.other_hits_inside_locus = []
        self.other_hits_outside_locus = []

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

    def clear(self):
        '''
        Clears everything in the KLocus object relevant to a particular assembly - gets it ready
        for the next assembly.
        '''
        self.blast_hits = []
        self.hit_ranges = IntRange()
        self.assembly_pieces = []
        self.identity = 0.0
        self.expected_hits_inside_locus = []
        self.missing_expected_genes = []
        self.expected_hits_outside_locus = []
        self.other_hits_inside_locus = []
        self.other_hits_outside_locus = []

    def get_coverage(self):
        '''
        Returns the % of this K-locus which is covered by BLAST hits in the given assembly.
        '''
        return 100.0 * self.hit_ranges.get_total_length() / len(self.seq)

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

    def clean_up_blast_hits(self):
        '''
        This function removes unnecessary BLAST hits from self.blast_hits.
        For each BLAST hit, we keep it if it offers new parts of the K-locus. If, on the other
        hand, it lies entirely within an existing hit (in K-locus positions), we ignore it. Since
        we first sort the BLAST hits longest to shortest, this strategy will prioritise long hits
        over short ones.
        '''
        self.blast_hits.sort(key=lambda x: x.length, reverse=True)
        kept_hits = []
        k_range_so_far = IntRange()
        for hit in self.blast_hits:
            hit_range = hit.get_query_range()
            if not k_range_so_far.contains(hit_range):
                k_range_so_far.merge_in_range(hit_range)
                kept_hits.append(hit)
        self.blast_hits = kept_hits

    def get_match_uncertainty_chars(self):
        '''
        Returns the character code which indicates uncertainty with how this K-locus was found in
        the current assembly.
        '?' means the K-locus was found in multiple discontinuous assembly pieces.
        '-' means that one or more expected genes were missing.
        '+' means that one or more additional genes were found in the K-locus assembly parts.
        '*' means that at least one of the expected genes in the K-locus is low identity.
        '''
        uncertainty_chars = ''
        if len(self.assembly_pieces) > 1:
            uncertainty_chars += '?'
        if self.missing_expected_genes:
            uncertainty_chars += '-'
        if self.other_hits_inside_locus:
            uncertainty_chars += '+'
        if not all([x.over_identity_threshold for x in self.expected_hits_inside_locus]):
            uncertainty_chars += '*'
        return uncertainty_chars

    def get_match_uncertainty_words(self):
        '''
        Returns a string with similar info as get_match_uncertainty_chars, but uses words instead
        of single characters.
        '''
        uncertainty_words = []
        if len(self.assembly_pieces) > 1:
            uncertainty_words.append('broken')
        if self.missing_expected_genes:
            uncertainty_words.append('incomplete')
        if self.other_hits_inside_locus:
            uncertainty_words.append('extra_genes')
        if not all([x.over_identity_threshold for x in self.expected_hits_inside_locus]):
            uncertainty_words.append('low_gene_identity')
        return ' '.join(uncertainty_words)

    def get_length_discrepancy(self):
        '''
        Returns an integer of the base discrepancy between the K-locus in the assembly and the
        reference K-locus sequence.
        E.g. if the assembly match was 5 bases shorter than the reference, this returns -5.
        This function only applies to cases where the K-locus was found in a single contig.  In
        other cases, it returns None.
        '''
        earliest_piece, latest_piece, same_contig_and_strand = self.get_earliest_and_latest_pieces()
        if not same_contig_and_strand:
            return None
        a_start = earliest_piece.start
        a_end = earliest_piece.end
        k_start = earliest_piece.earliest_hit_coordinate()
        k_end = latest_piece.latest_hit_coordinate()
        expected_length = k_end - k_start
        actual_length = a_end - a_start
        return actual_length - expected_length

    def get_length_discrepancy_string(self):
        '''
        Returns the length discrepancy, not as an integer but as a string with a sign and units.
        '''
        length_discrepancy = self.get_length_discrepancy()
        if length_discrepancy is None:
            return 'n/a'
        length_discrepancy_string = str(length_discrepancy) + ' bp'
        if length_discrepancy > 0:
            length_discrepancy_string = '+' + length_discrepancy_string
        return length_discrepancy_string


    def get_earliest_and_latest_pieces(self):
        '''
        Returns the AssemblyPiece with the earliest coordinate (closest to the K-locus start) and
        the AssemblyPiece with the latest coordinate (closest to the K-locus end)
        '''
        earliest_piece = sorted(self.assembly_pieces, key=lambda x: x.earliest_hit_coordinate())[0]
        latest_piece = sorted(self.assembly_pieces, key=lambda x: x.latest_hit_coordinate())[-1]
        same_contig_and_strand = earliest_piece.contig_name == latest_piece.contig_name and \
                                 earliest_piece.strand == latest_piece.strand
        return earliest_piece, latest_piece, same_contig_and_strand

    def get_missing_gene_string(self):
        return ';'.join([x.split('__')[2] for x in self.missing_expected_genes])



class Assembly(object):
    def __init__(self, fasta_file):
        '''
        Loads in an assembly and builds a BLAST database for it (if necessary).
        '''
        self.fasta = fasta_file
        self.name = os.path.splitext(os.path.basename(fasta_file))[0]
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
    def __init__(self, assembly, contig_name, contig_start, contig_end, strand, blast_hits=[]):
        self.assembly = assembly
        self.contig_name = contig_name
        self.start = contig_start
        self.end = contig_end
        self.strand = strand
        self.blast_hits = blast_hits

    def __repr__(self):
        return self.assembly.name + '_' + self.get_header()

    def get_header(self):
        '''
        Returns a descriptive string for the FASTA header when saving this piece to file.
        '''
        nice_contig_name = get_nice_contig_name(self.contig_name)
        return nice_contig_name + '_' + str(self.start + 1) + '_to_' + str(self.end) + \
               '_' + self.strand + '_strand'

    def get_bandage_range(self):
        '''
        Returns the assembly piece in a Bandage path format.
        '''
        if is_contig_name_spades_format(self.contig_name):
            name = self.contig_name.split('_')[1]
        else:
            name = self.contig_name.split()[0]
        return '(' + str(self.start + 1) + ') ' + name + '+ (' + str(self.end) + ')'

    def get_sequence(self):
        '''
        Returns the DNA sequence for this piece of the assembly.
        '''
        seq = self.assembly.contigs[self.contig_name][self.start:self.end]
        if self.strand == '+':
            return seq
        else:
            return reverse_complement(seq)

    def get_length(self):
        '''
        Returns the sequence length for this piece.
        '''
        return self.end - self.start

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
        To be able to combine, pieces must be overlapping or directly adjacent and on the same
        strand.
        '''
        if self.contig_name != other.contig_name or self.strand != other.strand:
            return None
        combined = IntRange([(self.start, self.end)])
        combined.add_range(other.start, other.end)
        if len(combined.ranges) == 1:
            new_start, new_end = combined.ranges[0]
            combined_hits = self.blast_hits + other.blast_hits
            return AssemblyPiece(self.assembly, self.contig_name, new_start, new_end, self.strand,
                                 combined_hits)
        else:
            return None

    def overlaps(self, contig_name, start, end):
        '''
        Returns whether this assembly piece overlaps with the given parameters.
        '''
        return self.contig_name == contig_name and self.start < end and start < self.end

    def earliest_hit_coordinate(self):
        '''
        Returns the lowest query start coordinate in the BLAST hits.
        '''
        if not self.blast_hits:
            return None
        return sorted([x.qstart for x in self.blast_hits])[0]

    def latest_hit_coordinate(self):
        '''
        Returns the highest query end coordinate in the BLAST hits.
        '''
        if not self.blast_hits:
            return None
        return sorted([x.qend for x in self.blast_hits])[-1]



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
            if int_range[0] > int_range[1]:
                fixed_ranges.append((int_range[1], int_range[0]))
            elif int_range[0] < int_range[1]:
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
