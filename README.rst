CoVa
====

Variant analysis pipeline for Coronavirus genomes. Specifically, COVID19.
This pipeline performs following analyses on a collection of coronavirus genomes.
 
* Whole-genome Multiple Sequence Alignment
* Variant calling
* Diversity computation
* Phylogenetic inference
* Positive selection analysis

Following external programs are required:

- MAFFT        : Build Multiple Sequence Alignment
- FastTreeMP   : Build phylogeny
- Hyphy        : Perform selection analysis to identify sites under positive selection

Installation
-----------------

1. Create a python virtual environment, at a location of your choice. 

   a. Use the following command, or something equivalent.
   
   .. code-block:: bash

      python3 -m venv <path/name>

   Replace <path/name> by a path of your choice. 
   
   b. Activate this environment by running
   
   .. code-block:: bash

      source <path/name>/bin/activate

2. Clone this repository.

   .. code-block:: bash

      git clone https://github.com/A-Farhan/cova.git
 
3. Install the package from the binary distribution. 

   .. code-block:: bash
   
      pip install <repo-path>/dist/*.whl
      
   Replace <repo-path> by the full path of your cloned cova repository. 
   
   This step may throw an error related to pyqt5 instllation. In that case, try upgrading pip.
   
   .. code-block:: bash
   
      pip install --upgrade pip

4. Install external programs.

   a. ``cd`` to cova repository

   b. run ``./INSTALL.sh``. This is an interactive script. It will ask for the permission to proceed with each installation *viz*. MAFFT, FastTree and Hyphy. Appropriate response should be *y* (yes) OR *n* (no). If you have any of these programs already installed, then you can choose to skip its installation here.

   c. CoVa expects external programs to be available system-wide, by default. In other words, programs are invoked without providing the full path. If you installed external programs using the previous step, a new directory *local* will be created inside the repository. All 3 programs shall then be accessible from the path <repo-path>/local/bin/.

      You'll need to add the following line to your ~/.bashrc file.

      .. code-block:: bash

         export COVA_BIN_PATH=<repo-path>/cova/local/bin  

   d. You can also access your pre-installed external programs similarly without using full paths. Simply, add their respective paths to your PATH variable.

      .. code-block:: bash

         export PATH=$PATH:<full-path-to-external-program>
      
   e. You can also find instructions to download and install these programs from the following urls:

      	https://mafft.cbrc.jp/alignment/software/installation_without_root.html

      	http://www.microbesonline.org/fasttree/#Install

      	https://hyphy.org/download/

Run CoVa from the command-line
------------------------------

This pipeline is built to be run as a command-line tool **CoVa**

To check if the command is available, run the following command inside the virtual environment

.. code-block:: bash

   CoVa --help
   
CoVa runs in the curent directory by default. You can provide any directory-path through ``--indr`` option. 

CoVa expects a minimum input of a whole-genome multi-FASTA file under this path, named "genomes.fna" by default. Several commands require a reference accession [default: ``NC_045512``]. Make sure, your FASTA file has this accession, or any reference of your choice, included. If you wish to use a different accession, you'll need to change several other arguments from within the source, as these are not available directly from the command-line. Also, you'll need to replace reference data files in the ``<package-path>/cova/data/`` directory. These files follow NCBI genome-assembly file formats. 

To get familiar with CoVa, and its outputs, you can run CoVa on ``<repo-path>/example`` directory. You can also copy the input file from this directory into an empty directory of your choce, for a fresh run. 

CoVa has multiple subcommands, and these commands have their own arguments. To see these arguments, you can run, for example 

.. code-block:: bash

   CoVa msabuild --help
   
You can run these commands individually or in combination, with or without arguments, as all arguments have defaults. Combination works like this.

.. code-block:: bash

   CoVa --indr <full-path-to-input-directory> msabuild msaref msaunq
   
Not all combinations would work, of course. As many commands depend on the input from specific preceding command(s), it is a largely rigid chain. To run the entire pipeline, use the sub-command ``full``.

.. code-block:: bash

   CoVa full
   
which is equivalent to 

.. code-block:: bash
   
   CoVa msabuild msaref msaunq msap vcalpd annpv vcali div tree sel tabvs   

Sub-commands
------------

1. **MSABUILD**

   This command is a MAFFT wrapper to build whole-genome Multiple Sequence Alignments (MSA).
   To speed up the process, MSABUILD performs no more than 5 refinement iterations and to reduce the 
   memory requirement, particularly a problem with wide whole-genome alignments, it uses the ``--nomemsave`` option.

   Input:

   A multi-FASTA file of *unaligned* whole-genomes present in the working directory, named "genomes.fna" by default.

   Output:

   A multi-FASTA file of *aligned* whole-genomes present in the same directory, named "genome_aln.fna" by default.

2. **MSAD**

   MAFFT allows for addition of new sequences to pre-existing MSAs. CoVa makes use of this feature to simplify incorporation of incoming genomic data and update analysis results. To do so, the main command can be run with the flag ``--addseq``. To facilitate rest of the analysis without changing any arguments, the MSA is changed in place and a copy is kept for back up. All the other analysis files are updated without a backup. If you wish to retain previous analysis, you can separately copy these files to a directory. 

   Input:

   MSA file generated by ``msabuild``.

   A FASTA file of possibly multiple genome sequences to be incorporated in the above MSA.

   Output:

   Updated input MSA file. 

3. **MSAREF**

   Before we can call variants ( point mutations and deletions) relative to a reference, our MSA must be restricted to the sites present in this reference. That's the job of this command.

   Input:

   MSA file generated by ``msabuild``, and a reference accession included in this MSA.

   Output:

   A multi-FASTA file of the above MSA limited to sites present in the given reference.

4. **MSAUNQ**

   Since it is possible that the original set of unaligned sequences, or the reference-limited MSA has duplicate sequences, it may be of interest to remove these duplicate sequences before further analysis. It would serve to speed up certain downstream analysis and also, to present polytomies in the phylogeny. 

   Input:

   MSA file generated by ``msaref``.

   Output:

   A multi-FASTA file of the above MSA excluding any duplicate sequences.

   A tab-delimited table of duplicate genomes:
   
   column 1 - genome id included in the output MSA
   column 2 - ids of genomes identical to the one selected 

5. **MSAP**

   This command extracts nucleotide MSAs for all proteins/peptides-encoding regions from the whole-genome reference limited MSA. These MSAs are later used for selection analysis on individual proteins.

   Input:

   MSA file generated by ``msaref``.

   A directory path to store output MSAs.

   Output:

   Nucleotide MSA files of individual protein/peptide-encoding regions.

6. **VCALPD**

   Variant CALling ( Point mutations / Deletions).  

   Input:

   MSA file generated by ``msaunq``.

   Output:

   Point mutation table with 1 row per variant and 1 column per genome, except the first 2 columns are for 1-indexed genomic coordinate and reference allele respectively.

   Deletion table with 1 row per deletion, and following columns:

      a. pos - 1-indexed genomic coordinate of the first base of deletion
      b. ref - deleted reference sequence
      c. len - length of deletion
      d. id  - Bits for absence(0) OR presence(1) of deletion in the respective genome

      There is one id column for every genome in the MSA. 

7. **ANNPV**

   ANNotation of Point mutation Variants.

   Input:

   Point mutation table generated by ``vcalpd``.

   Output:

   A tab-delimited table with following columns:

   1) protein_id - protein's accession in the reference genome
   2) name 	     - common name or abbreviation for the protein
   3) position 	 - 1-indexed genomic position
   4) ref_base   - nucleotide at the above position in the reference
   5) var_base   - a different allele at this position in some genome
   6) old_codon  - codon at this position in the protein-coding sequence of reference
   7) new_codon  - modified codon due to nucleotide substitution in some genome
   8) aa_change  - amino acid change due to this substitution
   9) genomes 	 - comma-separated list of genome ids with this variant

8. **TABVS**

  It may be of interest to characterize each isolate in terms of its unique variants and the variants 
  that it shares with the others, for further analyses. These results are summarized by this command. 
  Also, only non-synonymous variants are considered, in the interest of readability of the output table. 

  Input:

  Point mutation table generated by ``vcalpd``.
  
  Annotated point mutation table generated by ``annpv``.

  Output:

  A tab-delimited table with 1 row per genome and with following columns:
  
  a. genome    - genome id 
  b. #variants - total number of variants in the genome
  c. #shared   - number of shared non-synonymous variants
  d. #unique   - number of unique non-synonymous variants
  e. shared    - comma-separated list of shared variants
  f. unique    - comma-separated list of unique variants

9. **VCALI**

   Variant Calling for Insertions relative to a reference.

   Input:

   MSA file generated by ``msabuild``.

   Output:

   A tab-delimited table with 1 row per insertion and following columns:

   a. pos - 1-indexed genomic position of the reference base in the immediate left of the insertion
   b. ref - the reference base at the above position
   c. id  - either the reference base, if no insertion is present, OR an insertion sequence in the
      respective genome

   There is 1 id column for every genome.

10. **DIV**

   This command calculates Nucleotide Diversity for the whole-genome, as well as for all proteins/peptide-encoding regions. Nucleotide Diversity is the average pairwise-difference per base.

   Input:

   whole-genome MSA generated by ``msaref``.

   MSAs of protein/peptide-encoding regions generated by ``msap``.

   Output:

   A comma-delimited table. First row is for the whole-genome and following rows are for other regions.
   First column is the region's name and second column is for its nucleotide diversity.

11. **TREE**

    This command builds whole-genome based phylogeny using FastTree and plots a tree using python ETE3 module. The date and location information, if available, can be used to annotate the tree.

    Input:

    whole-genome MSA generated by ``msaref``.

    Output:

    Output tree generated by FastTree in NEWICK format.

    PNG image file for the above tree.

12. **SEL**

    This command runs HYPHY FUBAR which perform selection analysis on protein-encoding regions by estimating synonymous and non-synonymous rates. It also identifies putative sites under positive selection. 

    Input:

    MSAs generated by ``msap``.

    Phylogeny tree generated by ``tree``.

    Output:

    Output files generated by FUBAR.

    A comma-delimited table of *rates* with 1 row per protein and following columns:

    a. protein 	- common name or abbreviation for the protein
    b. exp_subs - expected substitution rate
    c. syn 	- synonymous rate
    d. nonsyn 	- non-synonymous rate
    e. dnds 	- (nonsyn-syn) 

    A comma-delimited table of *sites* with 1 row per site and following columns:

    a. protein 	 - common name or abbreviation for the protein
    b. site 	 - 1-indexed position in the protein
    c. syn 	 - site-specific synonymous rate
    d. nonsyn 	 - site-specific non-synonymous rate
    e. post_prob - posterior probability (nonsyn > syn)

External commands
---------------------------------

**MAFFT**

Cova runs the following MAFFT command

.. code-block:: bash

   mafft --quiet --nomemsave --maxiterate 5 --thread <ncpu> <infile>

**FastTree**

FastTree in cova was built from the source with

* Double-precision: improves branch length precision for highly similar sequences, AND
* OpenMP: allows multi-threading for faster computations 

using the following command

.. code-block:: bash

   gcc -DUSE_DOUBLE -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
   
Cova runs the following FastTree command. To reduce memory requirement, split supports are not calculated as these are not very informative for wide alignments. 

.. code-block:: bash

   FastTree -quiet -nt -mlnni 4 -nosupport

**Hyphy**

Cova makes use of Hyphy's **FUBAR** program to do selection analysis and identify sites under positive selection. It runs FUBAR as below

.. code-block:: bash

   hyphy fubar --alignment <msafile> --tree <treefile> --cache <cachefile>
