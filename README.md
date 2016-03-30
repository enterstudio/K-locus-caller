# K locus caller

This is a tool which reports information about the K locus for genome assemblies of capsule locus-containing bacterial species. It will help a user to decide whether their sample has a known or novel K locus.

It carries out the following for each input assembly:
* BLAST for all known K locus nucleotide sequences (using `blastn`) to identify the best match ('best' defined as having the highest coverage).
* Extract the region(s) of the assembly which correspond to the BLAST hits (i.e. the K locus sequence in the assembly) and save it to a FASTA file.
* BLAST for all known K locus genes (using `tblastn`) to identify which expected genes (genes in the best matching K locus) are present/missing and whether any unexpected genes (genes from other K loci) are present.
* Output a summary to a table file.

In cases where your input assembly closely matches a known K locus, this tool should make that obvious. When your assembly has a novel type, that too should be clear. However, this tool cannot reliably extract or annotate K locus sequences for totally novel types – if this tool indicates a novel K locus is present then extracting an annotating the sequence is up to you! Poor assemblies can also confound the results, so be sure to closely examine any case where the K locus sequence in your assembly is broken into multiple pieces.


## Table of Contents

* [Quick version (for the impatient)](https://github.com/rrwick/K locus-caller#quick-version-for-the-impatient)
* [Installation](https://github.com/rrwick/K locus-caller#installation)
* [Input files](https://github.com/rrwick/K locus-caller#input-files)
  * [K locus sequences](https://github.com/rrwick/K locus-caller#K locus-sequences)
  * [Gene protein sequences](https://github.com/rrwick/K locus-caller#gene-protein-sequences)
  * [K locus/gene table](https://github.com/rrwick/K locus-caller#K locusgene-table)
* [Standard output](https://github.com/rrwick/K locus-caller#standard-output)
  * [Basic](https://github.com/rrwick/K locus-caller#basic)
  * [Verbose](https://github.com/rrwick/K locus-caller#verbose)
* [Output files](https://github.com/rrwick/K locus-caller#output-files)
  * [Summary table](https://github.com/rrwick/K locus-caller#summary-table)
  * [K locus matching sequences](https://github.com/rrwick/K locus-caller#K locus-matching-sequences)
* [Example results and interpretation](https://github.com/rrwick/K locus-caller#example-results-and-interpretation)
  * [Very close match](https://github.com/rrwick/K locus-caller#very-close-match)
  * [More distant match](https://github.com/rrwick/K locus-caller#more-distant-match)
  * [Broken assembly](https://github.com/rrwick/K locus-caller#broken-assembly)
  * [Poor match](https://github.com/rrwick/K locus-caller#poor-match)
* [Advanced options](https://github.com/rrwick/K locus-caller#advanced-options)
* [License](https://github.com/rrwick/K locus-caller#license)


## Quick version (for the impatient)

This tool needs the following input files to run (included in this repository):
* A FASTA file with nucleotide sequences for your known K loci
* A FASTA file with proteins sequences for the genes in your known K loci
* A tab-delimited file specifying which genes go in which K loci
* One or more Klebsiella assemblies in FASTA format

Example command:

`k_locus_caller.py -a path/to/assemblies/*.fasta -k k_loci_refs.fasta -g k_loci_gene_list.txt -s genes.fasta -o output_directory`

For each input assembly file, this tool will identify the closest known K locus type and report information about the corresponding locus genes.

It generates the following output files:
* A FASTA file for each input assembly with the nucleotide sequences matching the closest K locus
* A table summarising the results for all input assemblies

Character codes in the output indicate problems with the K locus match:
* `?` = the match was not in a single piece, possible due to a poor match or discontiguous assembly.
* `-` = genes expected in the K locus were not found.
* `+` = extra genes were found in the K locus.
* `*` = one or more expected genes was found but with low identity.


## Installation

No explicit installation is required – simply clone (or download) from GitHub and run the script.

It uses [BLAST+](http://www.ncbi.nlm.nih.gov/books/NBK279690/), so that tool must also be installed (specifically the commands `makeblastdb`, `blastn` and `tblastn`). BLAST+ can usually be easily installed using a package manager such as [Homebrew](http://brew.sh/) (on Mac) or [apt-get](https://help.ubuntu.com/community/AptGet/Howto) (on Ubuntu and related Linux distributions).


## Input files

#### K locus sequences

This is a FASTA file containing the nucleotide sequences of each known K locus. The header for each sequence is simply the K locus name, e.g. K1, K2, etc.

Example:
```
>K1
ATGAATATGGCGAATTTGAAAGCGGTTATTCCGGTCGCAGGACTAGGCATGCATATGCTG
CCGGCCACAAAGGCAATTCCAAAGGAGATGCTGCCGATCGTTGATAAGCCAATGATTCAG
...
CCGCGATCTGTTTGGTAACGATTAA
>K2
ATGAATATGGCGAATTTGAAAGCGGTTATTCCGGTCGCAGGACTAGGCATGCATATGCTG
CCGGCCACAAAAGCAATTCCAAAGGAGATGCTGCCGATCGTTGATAAGCCAATGATTCAG
...
```

#### Gene protein sequences

This is a FASTA file containing the protein sequences for each protein in all known K loci. The headers must be in an [SRST2](https://github.com/katholt/srst2)-style format, which assigns each sequence to a cluster. The header contains four parts separated by double underscores. The important parts for this program are the first (cluster ID number) and the third (allele name).

Example:
```
>36__K1__K1-CDS1-galF__00064
MANLKAVIPVAGLGMHMLPATKAIPKEMLPIVDKPMIQYIVDEIVAAGIKEIVLVTHSSK
NAVENHFDTSYELEALLEQRVKRQLLAEVQAICPPGVTIMNVRQAQPLGLGHSILCARPV
...
AIAELAKKQSVDAMLMTGESYDCGKKMGYMQAFVTYGMRNLKEGAKFRESIKKLLA*
>43__K1__K1-CDS2-cpsACP__00003
MNWQLISFFGDSTVLLPSAAALFIVLMLRKTSRLLAWQWSLLFGITGAIVCASKLAFMGW
GLGIRELDYTGFSGHSALSAAFWPIFLWLLSARFSAGLQKAAVATGYILAAVVGYSRLVI
...
```

#### K locus/gene table

This is a tab-delimited table specifying which alleles are in which known K locus types. It has no header line and only two columns: the K locus name (which matches the headers in the K loci FASTA input) and the full allele name (which matches the headers in the gene FASTA input).

Example:
```
K1	36__K1__K1-CDS1-galF__00064
K1	43__K1__K1-CDS2-cpsACP__00003
...
K1	21__K1__K1-CDS20-ugd__00015
K2	36__K2__K2-CDS1-galF__00065
K2	43__K2__K2-CDS2-__00004
...
```


## Standard output

#### Basic

This tool will write a simple line to stdout for each assembly:
* the assembly name
* the best K locus match
* character codes for any match problems

Example:
```
assembly_1: K2*
assembly_2: K4
assembly_3: K17?-*
```

#### Verbose

If run without the `-v` or `--verbose` option, this tool will give detailed information about each assembly including:
* Which K locus reference best matched the assembly
* Information about the nucleotide sequence match between the assembly and the best K locus reference:
  * % Coverage and % identity
  * Length discrepancy (only available if assembled K locus match is in one piece)
  * Contig names and coordinates for matching sequences
* Details about found genes:
  * Whether they were expected or unexpected
  * Whether they were found inside or outside the K locus matching sequence
  * % Coverage and % identity
  * Contig names and coordinates for matching sequences


## Output files

#### Summary table

This tool produces a single tab-delimited table summarising the results of all input assemblies. It has the following columns:
* **Assembly**: the name of the input assembly, taken from the assembly filename.
* **Best match K locus**: the K locus type which most closely matches the assembly, based on BLAST coverage.
* **Problems**: characters indicating issues with the K locus match. An absence of any such characters indicates a very good match.
  * `?` = the match was not in a single piece, possible due to a poor match or discontiguous assembly.
  * `-` = genes expected in the K locus were not found.
  * `+` = extra genes were found in the K locus.
  * `*` = one or more expected genes was found but with low identity.
* **Coverage**: the percent of the K locus reference which BLAST found in the assembly.
* **Identity**: the nucleotide identity of the BLAST hits between K locus reference and assembly.
* **Length discrepancy**: the difference in length between the K locus match and the corresponding part of the assembly. Only available if the K locus was found in a single piece (i.e. the `?` problem character is not used).
* **Expected genes in K locus**: a fraction indicating how many of the genes in the best matching K locus were found in the K locus part of the assembly.
* **Expected genes in K locus, details**: gene names and percent identity (from the BLAST hits) for the expected genes found in the K locus part of the assembly.
* **Missing expected genes**: a string listing the gene names of expected genes that were not found.
* **Other genes in K locus**: the number of unexpected genes (genes from K loci other than the best match) which were found in the K locus part of the assembly.
* **Other genes in K locus, details**: gene names and percent identity (from the BLAST hits) for the other genes found in the K locus part of the assembly.
* **Expected genes outside K locus**: the number of expected genes which were found in the assembly but not in the K locus part of the assembly (usually zero)
* **Expected genes outside K locus, details**: gene names and percent identity (from the BLAST hits) for the expected genes found outside the K locus part of the assembly.
* **Other genes outside K locus**: the number of unexpected genes (genes from K loci other than the best match) which were found outside the K locus part of the assembly.
* **Other genes outside K locus, details**: gene names and percent identity (from the BLAST hits) for the other genes found outside the K locus part of the assembly.

#### K locus matching sequences

For each input assembly, this tool produces a FASTA file of the region(s) of the assembly which correspond to the best K locus match. This may be a single piece (in cases of a good assembly and a strong match) or it may be in multiple pieces (in cases of poor assembly and/or a novel K locus). The file is named using the assembly name and then the best matching K locus name.


## Example results and interpretation

These examples show what the tool's results might look like in the output table. The gene details columns of the table have been excluded for brevity, as they can be quite long.

#### Very close match

Assembly | Best match K locus | Problems | Coverage | Identity | Length discrepancy | Expected genes in K locus | Expected genes in K locus, details | Missing expected genes | Other genes in K locus | Other genes in K locus, details | Expected genes outside K locus | Expected genes outside K locus, details | Other genes outside K locus | Other genes outside K locus, details
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
assembly_1 | K1 |  | 99.94% | 99.81% | -22 bp | 20 / 20 | ... |  | 0 |  | 0 |  | 2 | ...

This is a case where our assembly very closely matches a known K locus type. There are no characters in the 'Problems' column, the coverage and identity are both high, the length discrepency is low, and all expected genes were found with high identity. A couple of other low-identity K locus genes hits were elsewhere in the assembly, but that's not abnormal and no cause for concern.

Overall, this is a nice, solid match for K1.

#### More distant match

Assembly | Best match K locus | Problems | Coverage | Identity | Length discrepancy | Expected genes in K locus | Expected genes in K locus, details | Missing expected genes | Other genes in K locus | Other genes in K locus, details | Expected genes outside K locus | Expected genes outside K locus, details | Other genes outside K locus | Other genes outside K locus, details
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
assembly_2 | K1 | * | 99.84% | 95.32% | +97 bp | 20 / 20 | ... |  | 0 |  | 0 |  | 2 | ...

This case shows an assembly that also matches the K1 locus sequence, but not as closely as our previous case. The `*` character indicates that one or more of the expected genes falls below the identity threshold (default 95%). The 'Expected genes in K locus, details' columns, excluded here for brevity, would show the identity for each gene.

Our sample still almost certainly has a K locus type of K1, but it has diverged a bit more from our K1 reference, possibly due to mutation and/or recombination.

#### Broken assembly

Assembly | Best match K locus | Problems | Coverage | Identity | Length discrepancy | Expected genes in K locus | Expected genes in K locus, details | Missing expected genes | Other genes in K locus | Other genes in K locus, details | Expected genes outside K locus | Expected genes outside K locus, details | Other genes outside K locus | Other genes outside K locus, details
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
assembly_3 | K2 | ?- | 99.95% | 98.38% | n/a | 17 / 18 | ... | K2-CDS17-manB | 0 |  | 0 |  | 1 | ...

Here is a case where our assembly matched a known K locus type well (high coverage and identity) but with a couple of problems. First, the `?` character indicates that the K locus sequence was not found in one piece in the assembly. Second, one of the expected genes (K2-CDS17-manB) was not found in the gene BLAST search.

In cases like this, it is worth examining the case in more detail outside of this tool. For this example, such an examination revealed that the assembly was poor (broken into many small pieces) and the manB gene happened to be split between two contigs. So the manB gene isn't really missing, it's just broken in two. Our sample most likely is a very good match for K2, but the poor assembly quality made it difficult for this tool to determine that automatically.

#### Poor match

Assembly | Best match K locus | Problems | Coverage | Identity | Length discrepancy | Expected genes in K locus | Expected genes in K locus, details | Missing expected genes | Other genes in K locus | Other genes in K locus, details | Expected genes outside K locus | Expected genes outside K locus, details | Other genes outside K locus | Other genes outside K locus, details
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
assembly_4 | K3 | ?-* | 77.94% | 83.60% | n/a | 15 / 20 | ... | ... | 0 |  | 0 |  | 5 | ...

In this case, the tool did not find a close match to any known K locus sequence. The best match was to K3, but BLAST only found alignments for 78% of the K3 sequence, and only at 84% nucleotide identity. Five of the twenty K3 genes were not found, and the 15 which were found had low identity. The assembly sequences matching K3 did not come in one piece (indicated by `?`), possibly due to assembly problems, but more likely due to the fact that our sample is not in fact K3 but rather has some novel K locus that was not in our reference inputs.

A case such as this demands a closer examination outside of this tool. It is likely a novel K locus type, and you may wish to extract and annotate the K locus sequence from the assembly.


## Advanced options

Each of these options has a default and is not required on the command line, but can be adjusted if desired:

* `--start_end_margin`: this tool tries to identify whether the start and end of a K locus are present in an assembly and in the same contig. This option allows for a bit of wiggle room in this determination. For example, if this value is 10 (the default), a K locus match that is missing the first 8 base pairs will still count as capturing the start of the locus. If set to zero, then the BLAST hit(s) must extend to the very start and end of the K locus for this tool to consider the match complete.
* `--min_gene_cov`: the minimum required percent coverage for the gene BLAST search. For example if this value is 90 (the default), then a gene BLAST hit which only covers 85% of the gene will be ignored. Using a lower value will allow smaller pieces of genes to be included in the results.
* `--min_gene_id`: the mimimum required percent identity for the gene BLAST search. For example if this value is 70 (the default), then a gene BLAST hit which has only 65% amino acid identity will be ignored. A lower value will allow for more distant gene hits to be included in the results (possibly resulting in more genes in the 'Other genes outside K locus' category). A higher value will make the tool only accept very close gene hits (possibly resulting in low-identity K locus genes not being found and included in 'expected genes not found in K locus').
* `--low_gene_id`: the percent identity threshold for what counts as a low identity match in the gene BLAST search. This only affects whether or not the `*` character is included in the 'Problems'. Default is 95.
* `--min_assembly_piece`: the smallest piece of the assembly (measured in bases) that will be included in the output FASTA files. For example, if this value is 100 (the default), then a 50 bp match between the assembly and the best matching K locus reference will be ignored.
* `--gap_fill_size`: the size of assembly gaps to be filled in when producing the output FASTA files. For example, if this value is 100 (the default) and an assembly has two separate K locus BLAST hits which are only 50 bp apart in a contig, they will be merged together into one sequence for the output FASTA. But if the two BLAST hits were 150 bp apart, they will be included in the output FASTA as two separate sequences. A lower value will possibly result in more fragmented output FASTA sequences. A higher value will possibly result in more sequences being included in the K locus output.


## License

GNU General Public License, version 3
