#!/usr/bin/env python
'''
K locus caller

This is a tool which reports information about the K type for Klebsiella genome assemblies. It
will help a user to decide whether their Klebsiella sample has a known or novel K type, and if
novel, how similar it is to a known type.

This script needs the following input files to run:
* A FASTA file with nucleotide sequences for your known K types
* A FASTA file with proteins sequences for the genes in your known K types
* A tab-delimited file which specifies which genes go in which K types
* One or more Klebsiella assemblies in FASTA format

Example command:
k_locus_caller.py -a path/to/assemblies/*.fasta -k k_loci_refs.fasta -g k_loci_gene_list.txt
    -s genes.fasta -o output_directory

For each input assembly file, the script will identify the closest known K type and report
information about the corresponding K locus genes.

It generates the following output files:
* A FASTA file for each input assembly with the nucleotide sequences matching the closest K type
* A table summarising the results for all input assemblies

Character codes indicate problems with the K locus match:
* `?` indicates that the match was not in a single piece, possible due to a poor match or
      discontiguous assembly
* `-` indicates that genes expected in the K locus were not found
* `+` indicates that extra genes were found in the K locus
* `*` indicates that one or more expected genes was found but with low identity

Author: Ryan Wick
email: rrwick@gmail.com
'''

from __future__ import print_function
from __future__ import division
import argparse
import sys
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation


def main():
    '''
    Script execution starts here.
    '''
    args = get_arguments()
    check_for_blast()
    check_files_exist(args.assembly + [args.k_refs])
    fix_paths(args)
    temp_dir = make_temp_dir(args)
    k_ref_seqs, gene_seqs, k_ref_genes = parse_genbank(args.k_refs, temp_dir)
    all_genes = []
    for genes in k_ref_genes.itervalues():
        all_genes += genes
    k_refs = load_k_locus_references(k_ref_seqs, k_ref_genes) # type: dict[str, KLocus]
    create_table_file(args.out)
    for fasta_file in args.assembly:
        assembly = Assembly(fasta_file)
        best_k = get_best_k_type_match(assembly, k_ref_seqs, k_refs)
        find_assembly_pieces(assembly, best_k, args)
        protein_blast(assembly, best_k, gene_seqs, all_genes, args)
        output(args.out, assembly, best_k, args)
        if not args.no_seq_out:
            save_assembly_pieces_to_file(best_k, assembly, args.out)
    clean_up(k_ref_seqs, gene_seqs, temp_dir)
    sys.exit(0)



def get_arguments():
    '''
    Specifies the command line arguments required by the script.
    '''
    parser = argparse.ArgumentParser(description='K locus caller',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--assembly', nargs='+', type=str, required=True,
                        help='Fasta file(s) for Klebsiella assemblies')
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

def make_temp_dir(args):
    '''
    Makes the temporary directory, if necessary. Returns the temp directory path.
    '''
    temp_dir = os.path.join(os.path.dirname(args.out), 'temp')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    return temp_dir

def clean_up(k_ref_seqs, gene_seqs, temp_dir):
    '''
    Deletes the temporary FASTA files. If the temp directory is then empty, it is deleted too.
    '''
    os.remove(k_ref_seqs)
    os.remove(gene_seqs)
    if not os.listdir(temp_dir):
        os.rmdir(temp_dir)

def parse_genbank(genbank, temp_dir):
    '''
    This function reads the input Genbank file and produces two temporary FASTA files: one with the
    K loci nucleotide sequences and one with the gene sequences.
    It returns the file paths for these two FASTA files along with a dictionary that links genes to
    K loci.
    '''
    k_ref_genes = {}
    temp_prefix = 'temp_' + str(os.getpid()) + '_'
    k_ref_seqs_filename = os.path.join(temp_dir, temp_prefix + 'k_ref_seqs.fasta')
    gene_seqs_filename = os.path.join(temp_dir, temp_prefix + 'gene_seqs.fasta')
    k_ref_seqs = open(k_ref_seqs_filename, 'w')
    gene_seqs = open(gene_seqs_filename, 'w')
    for record in SeqIO.parse(genbank, 'genbank'):
        k_locus_name = ''
        for feature in record.features:
            if feature.type == 'source' and 'note' in feature.qualifiers:
                for note in feature.qualifiers['note']:
                    if note.startswith('K locus:'):
                        k_locus_name = note[8:].strip()
        if not k_locus_name:
            quit_with_error('Genbank record missing K locus name')
        if k_locus_name in k_ref_genes:
            quit_with_error('Duplicate reference K locus name: ' + k_locus_name)
        k_ref_genes[k_locus_name] = []
        k_ref_seqs.write('>' + k_locus_name + '\n')
        k_ref_seqs.write(add_line_breaks_to_sequence(str(record.seq), 60))
        
        gene_num = 1
        for feature in record.features:
            if feature.type == 'CDS':
                gene = Gene(k_locus_name, gene_num, feature, record.seq)
                k_ref_genes[k_locus_name].append(gene)
                gene_num += 1
                gene_seqs.write(gene.get_fasta())
    return k_ref_seqs_filename, gene_seqs_filename, k_ref_genes

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
    # type: (Assembly, str, dict[str, KLocus]) -> KLocus
    '''
    Searches for all known K types in the given assembly and returns the best match.
    Best match is defined as the K type for which the largest fraction of the K type has a BLAST
    hit to the assembly.
    '''
    for k_ref in k_refs.itervalues():
        k_ref.clear()
    blast_hits = get_blast_hits(assembly, k_refs_fasta)
    for hit in blast_hits:
        if hit.qseqid not in k_refs:
            quit_with_error('BLAST hit (' + hit.qseqid + ') not found in K locus references')
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
    This function uses the BLAST hits in the given K type to find the corresponding pieces of the
    given assembly. It saves its results in the KLocus
    '''
    if not k_locus.blast_hits:
        return
    assembly_pieces = [x.assembly_piece for x in k_locus.blast_hits]
    merged_pieces = merge_assembly_pieces(assembly_pieces)
    length_filtered_pieces = [x for x in merged_pieces if x.get_length() >= args.min_assembly_piece]
    if not length_filtered_pieces:
        return
    k_locus.assembly_pieces = fill_assembly_piece_gaps(length_filtered_pieces, args.gap_fill_size)

    # Now check to see if the biggest assembly piece seems to capture the whole locus. If so, this
    # is an ideal match.
    biggest_piece = sorted(k_locus.assembly_pieces, key=lambda x: x.get_length(), reverse=True)[0]
    start = biggest_piece.earliest_hit_coordinate()
    end = biggest_piece.latest_hit_coordinate()
    if good_start_and_end(start, end, k_locus.get_length(), args.start_end_margin):
        k_locus.assembly_pieces = [biggest_piece]

    # If it isn't the ideal case, we still want to check if the start and end of the K locus were
    # found in the same contig. If so, fill all gaps in between so we include the entire
    # intervening sequence.
    else:
        earliest, latest, same_contig_and_strand = k_locus.get_earliest_and_latest_pieces()
        start = earliest.earliest_hit_coordinate()
        end = latest.latest_hit_coordinate()
        if good_start_and_end(start, end, k_locus.get_length(), args.start_end_margin) and \
        same_contig_and_strand:
            gap_filling_piece = AssemblyPiece(assembly, earliest.contig_name, earliest.start,
                                              latest.end, earliest.strand)
            k_locus.assembly_pieces = merge_assembly_pieces(k_locus.assembly_pieces + \
                                                            [gap_filling_piece])
    k_locus.identity = get_mean_identity(k_locus.assembly_pieces)

def protein_blast(assembly, k_locus, gene_seqs, all_genes, args):
    '''
    Conducts a BLAST search of all known K locus proteins. Stores the results in the KLocus
    object.
    '''
    hits = get_blast_hits(assembly, gene_seqs, genes=all_genes)
    hits = [x for x in hits if x.query_cov >= args.min_gene_cov and x.pident >= args.min_gene_id]
    expected_hits = []
    for expected_gene in k_locus.genes:
        best_hit = get_best_hit_for_query(hits, expected_gene)
        if not best_hit:
            k_locus.missing_expected_genes.append(expected_gene.full_name)
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

def create_table_file(output_prefix):
    '''
    Creates the table file and writes a header line if necessary.
    If the file already exists and the header line is correct, then it does nothing (to allow
    multiple independent processes to append to the file).
    '''
    table_path = output_prefix + '_table.txt'
    if os.path.isfile(table_path):
        with open(table_path, 'r') as existing_table:
            first_line = existing_table.readline().strip()
            if first_line == ('Assembly\tBest match locus\tProblems\tCoverage\tIdentity\t'
                              'Length discrepancy\tExpected genes in locus\t'
                              'Expected genes in locus, details\tMissing expected genes\t'
                              'Other genes in locus\tOther genes in locus, details\t'
                              'Expected genes outside locus\t'
                              'Expected genes outside locus, details\tOther genes outside locus\t'
                              'Other genes outside locus, details'):
                return




    table = open(table_path, 'w')
    headers = []
    headers.append('Assembly')
    headers.append('Best match locus')
    headers.append('Problems')
    headers.append('Coverage')
    headers.append('Identity')
    headers.append('Length discrepancy')
    headers.append('Expected genes in locus')
    headers.append('Expected genes in locus, details')
    headers.append('Missing expected genes')
    headers.append('Other genes in locus')
    headers.append('Other genes in locus, details')
    headers.append('Expected genes outside locus')
    headers.append('Expected genes outside locus, details')
    headers.append('Other genes outside locus')
    headers.append('Other genes outside locus, details')
    table.write('\t'.join(headers))
    table.write('\n')
    table.close()

def output(output_prefix, assembly, k_locus, args):
    '''
    Writes a line to the output table describing all that we've learned about the given K locus and
    writes to stdout as well.
    '''
    uncertainty_chars = k_locus.get_match_uncertainty_chars()
    expected_genes_per = 100.0 * len(k_locus.expected_hits_inside_locus) / len(k_locus.genes)
    expected_genes_str = str(len(k_locus.expected_hits_inside_locus)) + ' / ' + \
                         str(len(k_locus.genes)) + \
                         ' (' + float_to_str(expected_genes_per) + '%)'
    missing_genes_str = str(len(k_locus.missing_expected_genes)) + ' / ' + \
                        str(len(k_locus.genes))
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

    table_path = output_prefix + '_table.txt'
    table = open(table_path, 'a')
    table.write('\t'.join(line))
    table.write('\n')
    table.close()

    if not args.verbose:
        print(assembly.name + ': ' + k_locus.name + uncertainty_chars)
    if args.verbose:
        print('Assembly: ' + assembly.name)
        print('    Best match locus: ' + k_locus.name)
        print('    Problems: ' + uncertainty_chars)
        print('    Coverage: ' + coverage_str)
        print('    Identity: ' + identity_str)
        print('    Length discrepancy: ' + k_locus.get_length_discrepancy_string())
        print('    Locus assembly pieces:')
        for piece in k_locus.assembly_pieces:
            print('        ' + piece.get_header() + ', ' + piece.get_sequence_short())
        print('    Expected genes in locus: ' + expected_genes_str)
        for hit in k_locus.expected_hits_inside_locus:
            print('        ' + str(hit))
        print('    Missing expected genes: ' + missing_genes_str)
        for gene in k_locus.missing_expected_genes:
            print('        ' + str(gene))
        print('    Other genes in locus: ' + str(len(k_locus.other_hits_inside_locus)))
        for hit in k_locus.other_hits_inside_locus:
            print('        ' + str(hit))
        print('    Expected genes outside locus: ' + \
              str(len(k_locus.expected_hits_outside_locus)))
        for hit in k_locus.expected_hits_outside_locus:
            print('        ' + str(hit))
        print('    Other genes outside locus: ' + \
              str(len(k_locus.other_hits_outside_locus)))
        for hit in k_locus.other_hits_outside_locus:
            print('        ' + str(hit))
        print()

def float_to_str(float_in):
    '''
    This function converts a float to a string in a special manner: if the float is an integer,
    the resulting string has no decimal point. Otherwise, one decimal point is used.
    '''
    if float_in == int(float_in):
        return str(int(float_in))
    else:
        return '%.1f' % float_in

def get_blast_hits(assembly, query, genes=None):
    '''
    Returns a list BlastHit objects for a search of the given query in the given assembly.
    '''
    if genes:
        command = ['tblastn']
    else:
        command = ['blastn', '-task', 'blastn']
    command += ['-db', assembly.fasta, '-query', query, '-outfmt',
                '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    if err:
        quit_with_error('blastn encountered an error:\n' + err)
    if genes:
        blast_hits = [GeneBlastHit(line, assembly, genes) for line in line_iterator(out)]
    else:
        blast_hits = [BlastHit(line, assembly) for line in line_iterator(out)]
    return blast_hits

def get_best_hit_for_query(blast_hits, query):
    '''
    Given a list of BlastHits, this function returns the best hit for the given query, based on
    bit score. It returns None if no BLAST hits match that query.
    '''
    matching_hits = [x for x in blast_hits if x.qseqid == query.full_name]
    if matching_hits:
        return sorted(matching_hits, key=lambda x: x.bitscore, reverse=True)[0]
    else:
        return None

def cull_conflicting_hits(hit_to_keep, blast_hits):
    '''
    This function returns a (potentially) reduced set of BLAST hits which excludes BLAST hits that
    overlap too much (same part of assembly) with the hit to keep.
    '''
    return [x for x in blast_hits if not x.conflicts(hit_to_keep)]


def cull_all_conflicting_hits(blast_hits):
    '''
    This function returns a (potentially) reduced set of BLAST hits where none of the remaining
    hits conflict.
    '''
    blast_hits.sort(key=lambda x: x.bitscore, reverse=True)
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

def save_assembly_pieces_to_file(k_locus, assembly, output_prefix):
    '''
    Creates a single Genbank file for all of the assembly pieces.
    Assumes all assembly pieces are from the same assembly.
    '''
    if not k_locus.assembly_pieces:
        return
    records = []
    for piece in k_locus.assembly_pieces:
        seq = Seq(piece.get_sequence(), generic_dna)
        record = SeqRecord(seq)
        record.features.append(SeqFeature(FeatureLocation(0, len(seq)), type='source'))
        hits_to_include = k_locus.expected_hits_inside_locus + k_locus.other_hits_inside_locus
        hits_to_include = [x for x in hits_to_include if piece.overlaps_other(x.assembly_piece)]
        for hit in hits_to_include:
            start = hit.sstart - piece.start
            end = hit.send - piece.start
            hit_strand = 1
            if hit.strand == '-':
                hit_strand = -1
            hit_feature = SeqFeature(FeatureLocation(start, end, strand=hit_strand), type='CDS')
            hit_feature.qualifiers['codon_start'] = '1'
            if 'gene' in hit.gene.feature.qualifiers:
                hit_feature.qualifiers['gene'] = hit.gene.feature.qualifiers['gene']
            if 'product' in hit.gene.feature.qualifiers:
                hit_feature.qualifiers['product'] = hit.gene.feature.qualifiers['product']
            if 'transl_table' in hit.gene.feature.qualifiers:
                hit_feature.qualifiers['transl_table'] = hit.gene.feature.qualifiers['transl_table']
            hit_feature.qualifiers['note'] = ['Created by K locus caller']
            hit_feature.qualifiers['translation'] = hit.prot_seq
            record.features.append(hit_feature)
        records.append(record)
    filename = output_prefix + '_' + assembly.name + '.gbk'
    SeqIO.write(records, filename, 'genbank')




    # out_file = open(filename, 'w')
    # for piece in k_locus.assembly_pieces:
    #     fasta_file.write('>' + assembly.name + '_' + piece.get_header() + '\n')
    #     fasta_file.write(add_line_breaks_to_sequence(piece.get_sequence(), 60))

def add_line_breaks_to_sequence(sequence, length):
    '''
    Wraps sequences to the defined length. All resulting sequences end in a line break.
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

def load_k_locus_references(fasta, k_ref_genes): # type: (str, str) -> dict[str, KLocus]
    '''
    Returns a dictionary of:
      key = K locus name
      value = KLocus object
    '''
    return {seq[0]: KLocus(seq[0], seq[1], k_ref_genes[seq[0]]) for seq in load_fasta(fasta)}

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
    return ';'.join([x.qseqid + ',' + str(x.pident) + '%' for x in gene_hit_list])


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
    Stores the BLAST hit output mostly verbatim. However, it does convert the BLAST ranges
    (1-based, inclusive end) to Python ranges (0-based, exclusive end).
    '''
    def __init__(self, hit_string, assembly):
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
        self.assembly_piece = AssemblyPiece(assembly, self.sseqid, self.sstart, self.send, self.strand, [self])

    def __repr__(self):
        return self.qseqid + ', Contig: ' + self.sseqid + ' (' + str(self.sstart) + '-' + \
               str(self.send) + ', ' + self.strand + ' strand), ' + \
               'Cov: ' + '%.2f' % self.query_cov + '%, ID: ' + '%.2f' % self.pident + '%'

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
    def __init__(self, hit_string, assembly, genes):
        BlastHit.__init__(self, hit_string, assembly)
        self.over_identity_threshold = False
        self.prot_seq = Seq(self.assembly_piece.get_sequence(), generic_dna).translate(table=11)
        for gene in genes:
            if self.qseqid == gene.full_name:
                self.gene = gene

    def conflicts(self, other):
        '''
        Returns whether or not this hit conflicts with the other hit.
        A conflict is defined as the hits overlapping by 50% or more of the shortest hit's length.
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
        min_length = min(self.send - self.sstart, other.send - other.sstart)
        frac_overlap = overlap / min_length
        return frac_overlap > 0.5



class KLocus(object):
    def __init__(self, name, seq, genes):
        self.name = name
        self.seq = seq
        self.genes = genes
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
        return 'K locus ' + self.name

    def get_length(self):
        '''
        Returns the K locus sequence length.
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
        Returns the % of this K locus which is covered by BLAST hits in the given assembly.
        '''
        return 100.0 * self.hit_ranges.get_total_length() / len(self.seq)

    def clean_up_blast_hits(self):
        '''
        This function removes unnecessary BLAST hits from self.blast_hits.
        For each BLAST hit, we keep it if it offers new parts of the K locus. If, on the other
        hand, it lies entirely within an existing hit (in K locus positions), we ignore it. Since
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
        Returns the character code which indicates uncertainty with how this K locus was found in
        the current assembly.
        '?' means the K locus was found in multiple discontinuous assembly pieces.
        '-' means that one or more expected genes were missing.
        '+' means that one or more additional genes were found in the K locus assembly parts.
        '*' means that at least one of the expected genes in the K locus is low identity.
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

    def get_length_discrepancy(self):
        '''
        Returns an integer of the base discrepancy between the K locus in the assembly and the
        reference K locus sequence.
        E.g. if the assembly match was 5 bases shorter than the reference, this returns -5.
        This function only applies to cases where the K locus was found in a single contig. In
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
        Returns the AssemblyPiece with the earliest coordinate (closest to the K locus start) and
        the AssemblyPiece with the latest coordinate (closest to the K locus end)
        '''
        earliest_piece = sorted(self.assembly_pieces, key=lambda x: x.earliest_hit_coordinate())[0]
        latest_piece = sorted(self.assembly_pieces, key=lambda x: x.latest_hit_coordinate())[-1]
        same_contig_and_strand = earliest_piece.contig_name == latest_piece.contig_name and \
                                 earliest_piece.strand == latest_piece.strand
        return earliest_piece, latest_piece, same_contig_and_strand

    def get_missing_gene_string(self):
        return ';'.join(self.missing_expected_genes)



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
        piece. If they can't, it returns None.
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

    def overlaps_other(self, other):
        '''
        Conducts an overlap test against another assembly piece.
        '''
        return self.overlaps(other.contig_name, other.start, other.end)

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
    This class contains one or more integer ranges. Overlapping ranges will be merged together.
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



class Gene(object):
    '''
    This class prepares and stores a gene taken from the input Genbank file.
    '''
    def __init__(self, k_locus_name, num, feature, k_locus_seq):
        self.k_locus_name = k_locus_name
        self.feature = feature
        gene_num_string = str(num).zfill(2)
        self.full_name = k_locus_name + '_' + gene_num_string
        if 'gene' in feature.qualifiers:
            self.full_name += '_' + feature.qualifiers['gene'][0]
        self.nuc_seq = feature.extract(k_locus_seq)
        self.prot_seq = self.nuc_seq.translate(table=11)

    def get_fasta(self):
        '''
        Returns the FASTA version of this gene: a header line followed by sequence lines (of
        protein sequence) ending in a line break.
        '''
        return '>' + self.full_name + '\n' + \
               add_line_breaks_to_sequence(str(self.prot_seq), 60)




if __name__ == '__main__':
    main()
