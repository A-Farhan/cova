CoVa
====

Variant analysis pipeline for Coronavirus genomes. Specifically, COVID19.
This pipeline performs following analyses on a collection of coronavirus genomes.
 
* Whole-genome Multiple Sequence Alignment
* Variant calling
* Diversity computation
* Phylogenetic inference
* Positive selection analysis

The following external programs must be present to run CoVa

* MAFFT 		   : Build Multiple Sequence Alignment
* FastTreeMP : Build phylogeny
* Hyphy 		   : Selection analysis to identify sites under positive selection

Installation of external programs
---------------------------------

**MAFFT**

You can follow the instructions from https://mafft.cbrc.jp/alignment/software/installation_without_root.html
to download and install MAFFT. 

Cova runs the following MAFFT command

.. code-block:: bash

   mafft --quiet --nomemsave --maxiterate 5 --thread <ncpu> <infile>
   
**FastTree**

Instructions for downloading and installing are available at 
http://www.microbesonline.org/fasttree/#Install

You can either download the executables OR build from the source following the instructions
from the above link. 

The following command builds FastTree from the source file (FastTree.c) with

* Double-precision: improves branch length precision for highly similar sequences, AND
* OpenMP: allows multi-threading for faster computations 

.. code-block:: bash

   gcc -DUSE_DOUBLE -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
   
Cova runs the following FastTree command. To reduce memory requirement, split supports are not calculated as these are not very informative for wide alignments. 

.. code-block:: bash

   FastTree -quiet -nt -mlnni 4 -nosupport
   
**Hyphy**

To setup Hyphy, see https://hyphy.org/download/

Cova makes use of Hyphy's **FUBAR** program to do selection analysis and identify sites under positive selection. It runs FUBAR as below

.. code-block:: bash

   hyphy fubar --alignment <msafile> --tree <treefile> --cache <cachefile>
   
CoVa installation
-----------------

1. Create a python virtual environment, at a location of your choice. 

   a. Use the following command, or something equivalent.
   
   .. code-block:: bash

      python3 -m venv <path/name>
   
   b. Activate this environment by running
   
   .. code-block:: bash

      source <path/name>/bin/activate

2. Clone this repository.

   .. code-block:: bash

      git clone https://github.com/A-Farhan/cova.git
 
   Copy wheel to a location of your choice.
   
   .. code-block:: bash
   
      cp <path-to-cova-repo>/dist/cova-0.0.1-py3-none-any.whl <location-of-choice>
      
      cd <location-of-choice>
 
3. Install the package from the wheel

   .. code-block:: bash
   
      pip install <wheel-file>
      
4. Open python and test your installation by importing cova

   .. code-block:: python
   
      import cova
      
Run CoVa from the command-line
------------------------------

This pipeline is built to be run as a command-line tool **CoVa**

To check if the command is available, run the following command inside the virtual environment

.. code-block:: bash

   CoVa --help
   
CoVa runs in the curent directory by default. You can provide any directory-path through ``--indr`` option. 

CoVa expects a minimum input of a whole-genome multi-FASTA file under this path, named "genomes.fna" by default. Several commands require a reference accession [default: ``NC_045512``]. Make sure, your FASTA file has this accession, or any reference of your choice, included. If you wish to use a different accession, you'll need to change several other arguments from within the source, as these are not available directly from the command-line. Also, you'll need to replace reference data files in the ``<package-path>/cova/data/`` directory. These files follow NCBI genome-assembly file formats. 

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

To get familiar with these commands, and their outputs, you can run CoVa on ``<path-to-cova-repo>/example`` directory. You can also copy the input file from this directory into an empty directory of your choce, for a fresh run. 
