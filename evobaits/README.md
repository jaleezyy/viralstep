# EvoBaits extension for ViralSTEP

Below are details on the individual scripts used as part of EvoBaits. Future work will harmonize the code to be executable via `viralstep` as `viralstep evobaits`.

## Other software

Ancestral reconstruction was done using HyPhy's Ancestral Reconstruction standalone (https://github.com/veg/hyphy-analyses/tree/master/AncestralSequences).

From HyPhy standalone, the two files used in EvoBaits is:
* `AncestralSequences.bf`
* `FitMG94.bf`

## Imported scripts

The following scripts provide functions that can be used across core and auxiliary scripts:
* `convert.py`
* `random_seqeunces.py`

## Core scripts

All shell scripts ending with `.sh` are used to execute ViralSTEP for coverage assessment.

Python: 

`dot_cluster_framework.py`

```
usage: dot_cluster_framework.py [-h] -i INPUT_FASTA_FILE [-o OUTPUT_DIRECTORY]
                                [-d INPUT_FAMILY_FILE] [-c {nucl,prot}]
                                [-p TRAINING_FILE] [-r] [-k KMER]
                                [--overwrite] [--cleanup] [--test] [--log]
                                [-t THREADS]

Perform alignment based clustering based on predicted open reading frames
(ORFs). If taxonomy metadata provided in ViralSTEP format, the clustering will
be broken down per virus family.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTA_FILE, --input INPUT_FASTA_FILE
                        Input FASTA file
  -o OUTPUT_DIRECTORY, --outdir OUTPUT_DIRECTORY
                        Output directory containing FASTA files per determined
                        ORF.
  -d INPUT_FAMILY_FILE, --database INPUT_FAMILY_FILE
                        ViralSTEP classificaton database file. Will be used to
                        isolate sequences between virus families to produce a
                        framework.
  -c {nucl,prot}, --content {nucl,prot}
                        Seqeunce content for design. Choose between nucleotide
                        (nucl) or amino acid/protein (prot). Default is
                        'nucl'.
  -p TRAINING_FILE, --prodigaltf TRAINING_FILE
                        Use existing Prodigal training file. Useful for
                        replicating clusters.
  -r, --realign-unclustered
                        Update mash sketches and realign unclustered
                        sequences. Will loop until no further reductions are
                        possible. If any unclustered sequences remain, they
                        will be output as a FASTA file.
  -k KMER, --kmer KMER  K-mer size for MASH alignment. Default is 21, matching
                        default parameters for MASH
  --overwrite           Overwrite any existing intermediate or temporary
                        files/directories
  --cleanup             Delete temporary directories or files produced during
                        the run.
  --test                Run test function to ensure basic Cluster and
                        Framework objects are operational. Will exit after
                        running!
  --log                 Produce log file(s) with any Cluster frameworks
                        produced
  -t THREADS, --threads THREADS
                        Number of threads. Default is 1.
```

`ancestral_targets.py`

```
usage: ancestral_targets.py [-h] -i INPUT_FILES [-o OUTPUT_DIR]
                            [--keep-intermediate] [--overwrite]

Take FASTA clusters and generate phylogenetic alignments/trees and perform
ancestral reconstruction using HyPhy developmental Ancestral Reconstuction
software (v0.2).

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILES, --input INPUT_FILES
                        Input FASTA file or directory of FASTAs. Will auto-
                        detect content.
  -o OUTPUT_DIR, --output OUTPUT_DIR
                        Output directory containing ancestral nodes
  --keep-intermediate   Keep intermediate temporary files produced in the run
  --overwrite           Overwrite any prior existing intermediate files. If
                        not provided, and intermediates found, they will be
                        used instead of generating new file(s)
```

`node_selection.py`
```
usage: node_selection.py [-h] -i INPUT_FASTA_FILE -t INPUT_TREE
                         [-o OUTPUT_DIR] [--mode {seq,median,node_median}]
                         [--rev] [--overwrite]

Run Newick_utils on a given tree alongside node sequences to parse and extract
seqeunces for bait design.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTA_FILE, --input INPUT_FASTA_FILE
                        Input FASTA file containing node sequences
  -t INPUT_TREE, --tree INPUT_TREE
                        Input Newick tree file
  -o OUTPUT_DIR, --outdir OUTPUT_DIR
                        Output directory. Default='node_parse'
  --mode {seq,median,node_median}
                        Parsing mode. 'Seq' = Sequenctial mode. Output FASTA
                        files will contain a single sequence that is closest
                        to the tree root. Additional FASTA output will include
                        one additional sequence until end of tree. If the '--
                        rev' flag is used, output will begin from the furthest
                        node from the root. 'Median' = Median mode. The median
                        distance from the root will be determined from the
                        tips of the tree. Output FASTA will contain all node
                        sequences closest (or furthest if '--rev') to the
                        median. Node Median = Same principle as Median mode,
                        but median will instead be calculated using the
                        distances of the nodes rather than tips
  --rev                 Reverse order of sequences pulled
  --overwrite           Overwrite existing output directory!
```

## Auxiliary scripts

`orf_correction.py`

```
usage: orf_correction.py [-h] -i INPUT_FASTA_FILE [-o OUTPUT_PREFIX]
                         [-c {nucl,prot}] [--offset OFFSET] [--all-stops]
                         [--keep-intermediate]

Given nucleotide or protein ORF sequences, examine for consistent characters
corresponding to the content and remove stop codon traces, if desired.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTA_FILE, --input INPUT_FASTA_FILE
                        Input FASTA file. Will auto-detect content.
  -o OUTPUT_PREFIX, --output OUTPUT_PREFIX
                        Output file prefix. Will be appended to the beginning
                        of input FASTA file (i.e., 'corrected_'input.fasta)
  -c {nucl,prot}, --content {nucl,prot}
                        Sequence content for input sequences. Choose between
                        'nucl' or 'prot'. Leave blank for auto-detection of
                        content per sequence.
  --offset OFFSET       Determine how many characters/codons will be examined
                        per sequence at the 5' and 3' ends to validate
                        start/stop codons. Default = 9
  --all-stops           Remove all traces of stop codons
  --keep-intermediate   Keep intermediate temporary files
```

`orf_frame_count.py`

```
usage: orf_frame_count.py [-h] -i INPUT_FASTA_FILE [-o OUTPUT_PREFIX]
                          [--count] [--keep-intermediate] [-f]

Count length of each nucleotide sequence to determine if overall size of ORF
sequences are correct i.e., are divisible by 3; count and correction options
provided.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTA_FILE, --input INPUT_FASTA_FILE
                        Input FASTA file.
  -o OUTPUT_PREFIX, --output OUTPUT_PREFIX
                        Output file prefix. Will be appended to the beginning
                        of input FASTA file (i.e., 'corrected_'input.fasta)
  --count               Only output sequences that do not have correct frame
                        i.e., are not divisible by 3. No output will be
                        generated!
  --keep-intermediate   Keep intermediate temporary files
  -f, --force           By default, no output is generated if no trimming is
                        done. This flag will force output to be created
                        otherwise, producing a copy of the input file
```

`determine_taxonomy.py`
```
usage: determine_taxonomy.py [-h] -i INPUT_FASTA_FILE [-o OUTPUT_FASTA_FILE]
                             [-d INPUT_DB_FILE] -r MANUAL_REF [-k KMER]
                             [--overwrite] [--retain-headers]
                             [--keep-intermediates] [-t THREADS]

Perform alignment per sequence in a given cluster to reference sequences and
determine likely taxonomy. Headers will be updated through this script.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTA_FILE, --input INPUT_FASTA_FILE
                        Input FASTA file with sequences needing taxonomic
                        identification
  -o OUTPUT_FASTA_FILE, --output OUTPUT_FASTA_FILE
                        Output FASTA file. By default an output filename will
                        be generated based on the input file
  -d INPUT_DB_FILE, --database INPUT_DB_FILE
                        ViralSTEP classificaton database file.
  -r MANUAL_REF, --reference MANUAL_REF
                        Single FASTA file containing all desired reference
                        sequences.
  -k KMER, --kmer KMER  Optional argument for KMA aligner. Will index and
                        align accordingly. Default k-mer size of 16 will be
                        used if not specified.
  --overwrite           Overwrite any existing intermediate or temporary
                        files/directories
  --retain-headers      If toggled, then old header names will be maintained
                        in the final corrected output. Underscore characters
                        (_) will be replaced with hyphens (-).
  --keep-intermediates  For debugging purposes, keep all intermediate files.
  -t THREADS, --threads THREADS
                        Number of threads. Default is 1.
```

`determine_ancestral_labels.py`
```
usage: determine_ancestral_labels.py [-h] -i1 INPUT_FASTA_ORIG -i2
                                     INPUT_FASTA_UPDATED [-it INPUT_TREE]
                                     [-ij INPUT_JSON] [-od OUT_DIR]
                                     [-op OUT_PRE] [--overwrite]

Translation script to update labels in HyPhy AncestralSequences JSON and TREE
output as generated via ViralSTEP clustering.

optional arguments:
  -h, --help            show this help message and exit
  -i1 INPUT_FASTA_ORIG, --input1 INPUT_FASTA_ORIG
                        Input FASTA file with unclear headers. This should be
                        the FASTA from the early ViralSTEP clustering.
  -i2 INPUT_FASTA_UPDATED, --input2 INPUT_FASTA_UPDATED
                        Input FASTA file with updated or taxonomically
                        identified headers. This should be the output from
                        'determine_taxonomy.py'.
  -it INPUT_TREE, --input-tree INPUT_TREE
                        Input tree file to update labels.
  -ij INPUT_JSON, --input-json INPUT_JSON
                        NOT SUPPORTED YET! Input JSON file to update labels.
  -od OUT_DIR, --output-dir OUT_DIR
                        Output directory for corrected files. Default in the
                        current working directory.
  -op OUT_PRE, --output-prefix OUT_PRE
                        Output prefix for resulting translated files. Default
                        will retain input 2 filename with appended
                        'corrected'.
  --overwrite           If output file(s) already exist, overwrite them.
```

