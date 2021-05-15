# GRIT

Genome Recoding Informatics Toolbox (GRIT)

Stephen Casper, scasper@college.harvard.edu

TODO other authors

## About

This repository accompanies the paper [TODO]()

The GRIT software is a pythonic platform for working with human genome data, specifically [GRCh38.p13](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/). GRIT is meant for human genome recoding of TAA-->TAG using NG base editors, but it can be easily adapted to find a variety of statistics and sequences on chromosomes, genes, and coding DNA sequences. 

## Installation and Configuration

GRIT requires python 3 and has only a few main dependencies, so running it should be straightforward in an environment of your own. After cloning the repo and activating a virtual environment, install requirements.

  ```bash
  pip install -r requirements.txt
  ```

Then data will need to be gathered. To fetch CRCh38.p13 data, the UNIX command line [utilities for Entrez](https://www.ncbi.nlm.nih.gov/books/NBK179288/) must be installed. Run the following commands. 

  ```bash
  cd data
  curl -O https://www.ncbi.nlm.nih.gov/books/NBK179288/bin/install-edirect.sh
  chmod +x install-edirect.sh
  ./install-edirect.sh
  y
  sudo apt install acedb-other
  sudo apt install ncbi-entrez-direct
  y
  chmod +x get_data.sh
  ./get_data.sh
  ```

## Using GRIT

GRIT is primarily meant to be a platform for writing your own functions and using them to gather bioinformatic and recoding data on the human genome. However, four functions are build in which demo GRIT features, count editing sites, find genes to recode, count editable sites, and getting genome-wide data on sites. To execute these, run the following.

  ```bash
  mkdir results
  python GRIT.py --function=demo --write_file=results/demo_out.txt
  python GRIT.py --function=count_total_sites --write_file=results/count_total_sites_out.txt
  python GRIT.py --function=count_editing_sites --write_file=results/count_editing_sites_out.txt
  python GRIT.py --function=find_genes_to_recode --write_file=results/find_genes_to_recode_out.txt
  python GRIT.py --function=get_all_site_data ----write_file=results/get_all_site_data.txt
  ```

## Working with GRIT

The code for the repository is both long and fairly well commented, so provided here is a high-level overview of its functionality. See the actual code as well as the outputs of the `demo()` function for more specific info. 

GRIT consists of several major sections in two `.py` files. The main file, `GRIT.py` contains key imports, sample code, and functions to replicate results. The second file, `GRIT_utill.py` contains imports, initializations, helper functions, the `Chromosome` class, and the `Gene` class. 

The helper functions at the top of `GRIT_util.py` are used for parsing data, reverse complementing a string, and ordering indices by centrality in a base editor activity window for base editor site selection. 

A `Chromosome` object is instantiated by passing in a number or letter giving the chromosome as well as arguments for the data directory, what to consider the editing windows for A and C base editors, and what size range of primer products to find. When one is instantiated, the `__init__()` function calls each method inside the class to instantiate a variety of sequence, site, and gene attributes. Sites are found that can be directly edited with a C base editor or edited with a daisy chain of A and C editors. See the output of `demo()` to see how these are formatted. `Gene` objects are generally meant to be instantiated from within the `Chromosome` class. When one is instantiated, the `__init__()` function will sequences and recoding sites. 

To extend GRIT to specific purposes beyond the ones in the functions provided, it is recommended that you write your own functions on top of the provided code to customize your analysis. 
