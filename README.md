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

 1. Clone this repository: `git clone https://github.com/SushiLab/Amplicon_Recipes.git`
 2. Add the `16S_pipeline.sh` script to your PATH by modifying your `.bashrc` file: `export PATH="$PATH:<path_to_16S_pipeline.sh>"`
 3. The pipeline requires USEARCH and CUTADAPT to be executable with the exact commands `usearch`and `cutadapt`. You can test that by typing both commands.  
 If this does not work, you can add aliases to your `.bashrc` file. Use a text editor to add the lines `alias usearch=<path_to_USEARCH>` and `alias usearch=<path_to_CUTADAPT>` to your `.bashrc` file and source it: `source ~/.bashrc`.

### Pipeline

The 16S rRNA gene amplicon pipeline processes demultiplexed pair-end `fastq` files and produces OTU/zOTU abundance tables through the following steps:

 1. Merging of pair-end reads.
 2. Quality filtering.
 3. Primer matching (optional but recommended).
 4. Dereplication.
 5. OTU clustering with the UPARSE algorithm (97% id).
 6. OTU clustering with the UCLUST algorithm (100% id).
 7. zOTU denoting with the UNOISE3 algorithm.
 8. Taxonomic annotation of OTUs and zOTUs against SILVA database with LCA approach (optional).
 9. Quantification of OTU and zOTU abundances.

### Usage

**No primer sequences:** will skip the primer matching. Only recommended if primer sequences have been already removed and reads trimmed.

`16S_pipeline.sh -input_f <input_folder> -output_f <output_folder> -db <path_to_SILVA_database.fasta>`

**Primer sequences present:** will include the primer matching. Recommended if primer sequences have not been yet removed.

`16S_pipeline.sh -input_f ./data/ -output_f ./out -db <path_to_SILVA_database.fasta> -primerF <forward_primer> -primerR <reverse_primer> -threads <num_threads>`

Example: using *SILVA_128_SSURef_Nr99* database and the *515F-Y / 806RB* primers for the V4 region:

`16S_pipeline.sh -input_f ./data/ -output_f ./out -db SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta -primerF GTGYCAGCMGCCGCGGTAA -primerR ATTAGAWACCCBNGTAGTCC -threads 10`


### Parameters

Mandatory options:

**-input_f** Path to input folder with demultiplexed raw reads fastq files.
Files must start with the sample identifier. Sample identifier is taken from the FASTQ file name by truncating at the first underscore or period.
R1 and R2 files for each sample must have the exact same file name but 'SampleID_*R1*.fastq' and 'SampleID_*R2*.fastq'.  
**-output_f** Path to output folder.

General options:

**-threads[1-N]** Number of threads used (default=1)  
**-help** Show this help

Merging:

**-pctid[0-100]** Minimum percentage id of alignment (default=90)  
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

### Results

The complete pipeline will produce the following output files:

Report:

`report.txt` Report of the whole pipeline with basic statistics

Merging:

`merged.fq` Fastq file with merged reads  
`merging.log` Logfile for the merging step

Filtering:

`filtered.fa` Fasta file with quality filetered reads  
`filtered_primermatch.fa` Fasta file with reads matchin the primers  
`filter.log` Logfile for the quality filtering step

De-replication:

`uniques.fa` Fasta file with de-replicated sequences  
`dereplication.log` Logfile for the quality de-replication step

Clustering / Denoising:

`otus_uparse.fa` Fasta file with OTU representative sequences  
`otutab_uparse.*` OTU table (3 available formats)  
`clustering.log` Logfile for the UPARSE clustering step  
`make_otutab_uparse.log` Logfile for the OTU table quantification step

`otus_unoise.fa` Fasta file with zOTU sequences  
`otutab_unoise.*` zOTU table (3 available formats)  
`denoising.log` Logfile for the denoising step  
`make_otutab_unoise.log` Logfile for the zOTU table quantification step
 
Taxonomic annotation:

`taxonomy_uparse_lca.txt` Taxonomic annotation of OTUs  
`taxonomy_unoise_lca.txt` Taxonomic annotation of zOTUs  
`taxsearch_uparse.tax` All hits to the taxonomic database for OTUs  
`taxsearch_unoise.tax` All hits to the taxonomic database for zOTUs  
`taxsearch_uparse.log` Logfile for the OTU taxonomic annotation step  
`taxsearch_unoise.log` Logfile for the zOTU taxonomic annotation step
