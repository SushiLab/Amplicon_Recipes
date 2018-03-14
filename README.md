# Amplicon recipes

Recipes for data analysis of amplicon sequencing data

**Available recipes:**  
 1. 16S rRNA gene amplicon pipeline (`16S_pipeline.sh`)

## 1. 16S rRNA gene amplicon pipeline (*16S_pipeline.sh*)

Data processing pipeline for 16S rRNA gene amplicon raw sequencing data

### Dependencies

This pipeline depends on [USEARCH](https://www.drive5.com/usearch/) and [CUTADAPT](https://github.com/marcelm/cutadapt).  
The pipeline has been tested with **USEARCH v.10.0.240** and **CUTADAPT v.1.9.1**.

### Installation

 1. Clone this repository
 2. Add the `16S_pipeline.sh` script to your PATH by modifying your `.bashrc` file: `export PATH="$PATH:<path_to_16S_pipeline.sh>"`
 3. Add both to your PATH by modifying your `.bashrc` file: `export PATH="$PATH:<path_to_USEARCH>"` and `export PATH="$PATH:<path_to_CUTADAPT>"`. Then source it: `source ~/.bashrc`.
 4. The pipeline requires USEARCH and CUTADAPT to be executable with the exact commands `usearch`and `cutadapt`. You can test that by typing both commands. If not working add an alias to your `.bashrc` file: `alias usearch=<path_to_USEARCH>` and `alias usearch=<path_to_CUTADAPT>`. Then source it: `source ~/.bashrc`.

### Pipeline

The 16S rRNA gene amplicon pipeline processes demultiplexed pair-end `fastq` files and produces OTU/zOTU abundance tables through the following steps:

 1. Merging of pair-end reads.
 2. Quality filtering.
 3. Primer matching (optional but recommended).
 4. Dereplication.
 5. OTU clustering with the UPARSE algorithm.
 6. zOTU denoting with the UNOISE3 algorithm.
 7. Taxonomic annotation of OTUs and zOTUs against SILVA database with LCA approach (optional).
 8. Quantification of OTU and zOTU abundances.

### Usage

**Minimal usage** (not recommended): will skip the primer matching and taxonomic annotation steps and run on a single CPU.

`16S_pipeline.sh -input_f <input_folder> -output_f <output_folder>`

**Recomended usage:** will include the primer matching and taxonomic annotation steps.

`16S_pipeline.sh -input_f ./data/ -output_f ./out -db <path_to_SILVA_database.fasta> -primerF <forward_primer> -primerR <reverse_primer> -threads <num_threads>`

Example: using *SILVA_128_SSURef_Nr99* database and the *515F-Y / 806R* primers for the V4 region:

`16S_pipeline.sh -input_f ./data/ -output_f ./out -db SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta -primerF GTGYCAGCMGCCGCGGTAA -primerR ATTAGAWACCCNDGTAGTCC -threads 10`

### Parameters
The following parameters may be used:

Mandatory options:

**-input_f** Path to input folder with demultiplexed raw reads fastq files.
Files must start with the sample identifier. Sample identifier is taken from the FASTQ file name by truncating at the first underscore or period.
R1 and R2 files for each sample must have the exact same file name but 'SampleID_*R1*.fastq' and 'SampleID_*R2*.fastq'.  
**-output_f** Path to output folder.

General options:

**-threads[1-N]** Number of threads used (default=1)  
**-help** Show this help

Merging:

**-pctid[0-100]** Percentage of maximum number of mismatchesas an integer (default=10)  
**-minoverlap [0-N]** Discard pair if alignment is shorter than given value (default=16)

Quality filtering:

**-maxee[0-N]** Expected errors: discard reads with expected errors > maxee (default=0.1)  
**-minlength[0-N]** Discard sequences with length < minlength (default=100)

Primer match: (skipped if -primerF or primerR options are missing)

**-primerF** Forward primer with IUPAC wildcard characters.  
**-primerR** Reverse primer with IUPAC wildcard characters (reverse-complement needed).  
**-minprimfrac[0-1]** Minimum fraction of the primer searched by cutadapt (default=1, i.e. the entire primer).  
**-maxmismatch[1-N]** Number of mismatches allowed by cutadapt in each primer (default=0)

Clustering:

**-minsize[0-1]** Minimum sequence abundance to be considered (default=1, i.e. include singletons)

Taxonomical annotation: (skipped if -db option is missing)

**-db** Path to database for taxonomical annotation (SILVA db suggested: https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz)  
**-tax_id [0-1]** Minimum identity for taxonomic search (default=0.90)
