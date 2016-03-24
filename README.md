# K-locus caller

This is a tool which reports information about the K-type for Klebsiella genome assemblies.  It will help a user to decide whether their Klebsiella sample has a known or novel K-type, and if novel, how similar it is to a known type.


## Quick version (for the impatient)

This script needs the following input files to run (included in this repository):
* A FASTA file with nucleotide sequences for your known K-types
* A FASTA file with proteins sequences for the genes in your known K-types
* A tab-delimited file which specifies which genes go in which K-types
* One or more Klebsiella assemblies in FASTA format

Example command:

`k_locus_caller.py -a path/to/assemblies/*.fasta -k k_loci_refs.fasta -g k_loci_gene_list.txt -s genes.fasta -o output_directory`

For each input assembly file, the script will identify the closest known K-type and report information about the corresponding K-locus genes.

It generates the following output files:
* A FASTA file for each input assembly with the nucleotide sequences matching the closest K-type
* A table summarising the results for all input assemblies

Character codes indicate problems with the K-locus match:
* `?` indicates that the match was not in a single piece, possible due to a poor match or discontiguous assembly
* `-` indicates that genes expected in the K-locus were not found
* `+` indicates that extra genes were found in the K-locus
* `*` indicates that one or more expected genes was found but with low identity


## Installation

No explicit installation is required - simply clone (or download) from GitHub and run the script.  The script uses BLAST+, and so that tool is required (specifically the commands `makeblastdb`, `blastn` and `tblastn`).


## Input files

### K-type sequences

### Gene protein sequences

### K-type/gene table


## Standard output

### Basic

This tool will write a simple line to stdout for each assembly:
* the assembly name
* the best K-type match
* character codes for any match problems

Example:
```
assembly_1: K2*
assembly_2: K4
assembly_3: K17?-*
```

### Verbose

If run without the `--verbose` option, this script will give detailed information about each assembly including:
* Which K-type reference best matched the assembly
* Information about the nucleotide sequence match between the assembly and the best K-type reference:
  * % Coverage and % identity
  * Length discrepancy (only available if assembled K-type match is in one piece)
  * Contig names and coordinates for matching sequences
* Details about found genes:
  * Whether they were expected or unexpected
  * Whether they were found inside or outside the K-locus matching sequence
  * % Coverage and % identity
  * Contig names and coordinates for matching sequences


## Output files

### Summary table

### K-locus matching sequences


## Advanced options

Each of these options has a default and is not required on the command line, but can be adjusted if desired:

* `--start_end_margin`
* `--min_gene_cov`
* `--min_gene_id`
* `--low_gene_id`
* `--min_assembly_piece`
* `--gap_fill_size`



## License

GNU General Public License, version 3
