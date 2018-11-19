#!/bin/bash

comm_exec=$(echo $(basename $0) $@)

# If no argument is passed
if [[ ! $@ =~ ^\-.+ ]]
then
  printf "No options have been passed. Options -input_f and -output_f are mandatory.\nUse '$(basename $0) -help' for help.\nExiting\n"
  exit
fi

# Usage function
function usage
{
cat <<EOF

-----------------------------------------------------------
 Data processing pipeline for 16S amplicon sequencing data 
 Guillem Salazar - Sunagawa Lab (guillems@ethz.ch)          
 2018                                                      
 Modified by Chris Field (fieldc@ethz.ch) 2018
-----------------------------------------------------------

Usage: $(basename $0) -input_f <INPUT_FOLDER> -output_f <OUTPUT_FOLDER> [options]


  -input_f      Path to input folder with demultiplexed raw reads fastq files.
                Files must start with the sample identifier. Sample identifier is taken from the FASTQ file name by truncating at the first underscore (_) or period (.).
                R1 and R2 files for each sample must have the exact same file name but 'SampleID_*R1*.fastq' and 'SampleID_*R2*.fastq'.
  -output_f     Path to output folder.

OPTIONS:

General options:

  -threads      [1-N] Number of threads used (default=1)
  -help         Show this help

Merging:
  
  -pctid        [0-100] Minimum percentage id of alignment (default=90)
  -minoverlap   [0-N] Discard pair if alignment is shorter than given value (default=16)
  
Quality filtering:
  
  -maxee        [0-N] Expected errors: discard reads with expected errors > maxee (default=0.1)
  -minlength    [0-N] Discard sequences with length < minlength (default=100)
  
Primer match: (skipped if -primerF or primerR options are missing)

  -primerF      Forward primer with IUPAC wildcard characters.
  -primerR      Reverse primer with IUPAC wildcard characters (reverse-complement needed).
  -minprimfrac  [0-1] Minimum fraction of the primer searched by cutadapt (default=1, i.e. the entire primer).
  -maxmismatch  [1-N] Number of mismatches allowed by cutadapt in each primer (default=0)
  
Clustering:

  -minsize      [0-1] Minimum sequence abundance to be considered (default=1, i.e. include singletons)

Taxonomical annotation: (skipped if -db option is missing)

  -db           Path to database for taxonomical annotation (SILVA db suggested: 'https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz')
  -tax_id       [0-1] Minimum identity for taxonomic search (default=0.90)

Defined community:

  -ref          Path to a fasta file of reference sequences for a defined community. If this option is given, only unclassifiable sequences will be de novo clustered.
  -ref_id        [0-1] Identity threshold for alignment to reference sequences (default=0.97)

EOF
}


# Default values
threads=1
pctid=90
maxee=0.1
minoverlap=16
minprimfrac=1
maxmismatch=0
minlength=100
minsize=1
tax_id=0.90
ref_id=0.97

while [ $# -gt 0 ]
do
    case "$1" in
    -input_f) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else input_f="$2"; fi; shift;;
	-output_f) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else output_f="$2"; fi; shift;;
	-db) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else db="$2"; fi; shift;;
	-primerF) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else primerF="$2"; fi; shift;;
	-primerR) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else primerR="$2"; fi; shift;;
	-threads) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else threads="$2"; fi; shift;;
	-help) usage; exit;;
	-pctid) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else pctid="$2"; fi; shift;;
	-maxee) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else maxee="$2"; fi; shift;;
        -minoverlap) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else minoverlap="$2"; fi; shift;; 
        -minlength) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else minlength="$2"; fi; shift;;
	-minprimfrac) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else minprimfrac="$2"; fi; shift;;
	-maxmismatch) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else maxmismatch="$2"; fi; shift;;
	-minsize) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else minsize="$2"; fi; shift;;
	-tax_id) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else tax_id="$2"; fi; shift;;
	-ref) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else ref="$2"; fi; shift;;
    -ref_id) if [[ $2 == -* ]] || [ -z "$2" ]; then echo "Option $1 needs an argument"; exit; else ref_id="$2"; fi; shift;;
	--) shift; break;;
	-*) printf "Unrecognized option $1\n'$(basename $0) -help' for help.\nExiting\n"
	    exit 1;;
	*)  printf "Unrecognized option $1\nAll options should start by '-'\n '$(basename $0) -help' for help.\nExiting\n"; exit 1;;	# terminate while loop
    esac
    shift
done


# check all mandatory options have arguments
if [[ -z ${input_f+x} || -z ${output_f+x} ]]; then printf "Not all mandatory options have been set.\nOptions -input_f -output_f are mandatory.\nExiting\n"; exit; fi

# check whether input and output directory exist
if [ ! -d $input_f ]; then printf "Input directory ($input_f) does not exist.\nExiting.\n"; exit; fi
if [ ! -d $output_f ]; then mkdir $output_f; printf "\nOutput directory ($output_f) did not exist. It has been created\n"; fi

# check variables suposed to be integers
if [[ ! "$maxmismatch" =~ ^[0-9]+$ ]]; then printf " -maxmismatch only accepts integers\nExiting.\n"; exit; fi
if [[ ! "$threads" =~ ^[0-9]+$ ]]; then printf " -threads only accepts integers\nExiting.\n"; exit; fi
if [[ ! "$pctid" =~ ^[0-9]+$ ]]; then printf " -pctid only accepts integers\nExiting.\n"; exit; fi
if [[ ! "$minsize" =~ ^[0-9]+$ ]]; then printf " -minsize only accepts integers\nExiting.\n"; exit; fi
if [[ ! "$minlength" =~ ^[0-9]+$ ]]; then printf " -minlength only accepts integers\nExiting.\n"; exit; fi
if [[ ! "$minoverlap" =~ ^[0-9]+$ ]]; then printf " -minoverlap only accepts integers\nExiting.\n"; exit; fi

# check variables suposed to be between >0 
if [[ $(echo $maxee'<'0 | bc -l) -eq 1 ]]; then printf " -maxee only accepts numeric positive values\nExiting.\n"; exit; fi

# check variables suposed to be between 0 and 1
if [[ $(echo $minprimfrac'<'0 | bc -l) -eq 1 ]] || [[ $(echo $minprimfrac'>'1 | bc -l) -eq 1 ]]; then printf " -minprimfrac only accepts values within 0-1\nExiting.\n"; exit; fi
if [[ $(echo $tax_id'<'0 | bc -l) -eq 1 ]] || [[ $(echo $tax_id'>'1 | bc -l) -eq 1 ]]; then printf " -tax_id only accepts values within 0-1\nExiting.\n"; exit; fi
if [[ $(echo $ref_id'<'0 | bc -l) -eq 1 ]] || [[ $(echo $ref_id'>'1 | bc -l) -eq 1 ]]; then printf " -ref_id only accepts values between 0 and 1\nExiting.\n"; exit; fi

# check reference file exists if -ref is given
if [[ -v ref ]]; then if [ ! -f $ref ]; then printf " Reference sequences file ($ref) does not exist.\nExiting.\n";exit;fi;fi

# Define some extra variables
MIN_F=$(echo "scale=1;${#primerF}*$minprimfrac" | bc | awk '{print int($1)}'); # Minimum length of primer overlap allowed in primer search (F primer).
MIN_R=$(echo "scale=1;${#primerR}*$minprimfrac" | bc | awk '{print int($1)}'); # Minimum length of primer overlap allowed in primer search (R primer).
ERR_F=$(echo "scale=8;$maxmismatch.4/${#primerF}" | bc) # Error tolerance for primer search: mismatches/primer_length (primer F).
ERR_R=$(echo "scale=8;$maxmismatch.4/${#primerR}" | bc) # Error tolerance for primer search: mismatches/primer_length (primer R).
usearch=`which usearch`
cutadapt=`which cutadapt`

if [ -z "$usearch" ]; then echo -e "\nUSEARCH not found. Include USEARCH in your path or create an alias. It should be executable by using the exact command 'usearch'\n"; exit 1; fi
if [ -z "$cutadapt" ]; then echo -e "\nCUTADAPT not found. Include CUTADAPT in your path or create an alias. It should be executable by using the exact command 'cutadapt'\n"; exit 1; fi

# check USEARCH version
version=$(usearch -version | sed "-es/usearch v10.*/v10/")
if [[ $version != "v10" ]]; then echo "usearch version too old, need v10"; fi


# Print options executed
function show_params
{
cat <<EOF

-- Command executed --

$comm_exec

-- Software versions --

usearch: $usearch
cutadapt: $cutadapt

---- Parameters ----
input_f: $input_f
output_f: $output_f
$(if [[ -v db ]]; then printf "db: $db"; else printf "db: not used"; fi)
$(if [[ -v ref ]]; then printf "Defined community references: $ref"; fi)
$(if [[ -v ref ]]; then printf "Reference identity threshold: $ref_id"; fi)
$(if [[ -v primerF ]]; then printf "primerF: $primerF"; else printf "primerF: not used"; fi)
$(if [[ -v primerR ]]; then printf "primerR: $primerR"; else printf "primerR: not used"; fi)
threads: $threads
minoverlap: $minoverlap
pctid: $pctid
maxee: $maxee
minlength: $minlength
minprimfrac: $minprimfrac
maxmismatch: $maxmismatch
minsize: $minsize
tax_id: $tax_id
--------------------

EOF
}

# Show parameters used
show_params

#####################
## PIPELINE STARTS ##
#####################

# Pipefail
set -Eeou pipefail

## MERGING PAIRED READS
if [ -e $output_f/merged.fq ]
then
    echo -e "\nMerged paired reads already exist. Skip this step.\n"
else
    echo -e "\nMerging paired reads...\n"
    $usearch -fastq_mergepairs $(ls -v $input_f/*R1*.fastq) -fastqout $output_f/merged.fq -fastqout_notmerged_fwd $output_f/unmerged_fwd.fq -fastqout_notmerged_rev $output_f/unmerged_rev.fq -fastq_minovlen ${minoverlap} -relabel @ -fastq_pctid ${pctid} -fastq_maxdiffs 300 -threads ${threads} &> $output_f/merging.log
    echo -e "\n...done merging reads.\n"
fi

if [ ! -e $output_f/merged.fq ]; then echo -e "\n${output_f}/merged.fq was not created. Merging failed. Exiting...\n"; exit 1; fi

## FILTERING MERGED READS

# Quality filtering
if [[ -e $output_f/filtered.fa ]]
then
    echo -e "\nFiltered reads by quality already exist. Skip this step.\n"
else
    echo -e "\nFiltering merged reads...\n"
    $usearch -fastq_filter $output_f/merged.fq -fastq_maxee ${maxee} -fastq_minlen ${minlength} -fastaout $output_f/filtered.fa -threads ${threads} &> $output_f/filter.log
    echo -e "\nDone filtering merged reads.\n"
fi

if [ ! -e $output_f/filtered.fa ]; then echo -e "\n${output_f}/filtered.fa was not created. Quality filtering failed. Exiting...\n"; exit 1; fi

# Primer match filtering
if [ -e $output_f/filtered_primermatch.fa ]
then
    echo -e "\nFiltered reads after primer match already exist. Skip this step.\n"
else

    if [[ -z ${primerF+x} || -z ${primerR+x} ]]
    then
        echo -e "\nPair of primers not provided. Skipping primer match step.\n"
    else
    
        echo -e "\nSelecting reads with primer matches...\n"
        $cutadapt --discard-untrimmed -g ${primerF} -O ${MIN_F} -e ${ERR_F} -m ${minlength} -o $output_f/filtered.tmp.F.fa $output_f/filtered.fa &> $output_f/primermatch.log
        $cutadapt --discard-untrimmed -a ${primerR} -O ${MIN_R} -e ${ERR_R} -m ${minlength} -o $output_f/filtered_primermatch.fa $output_f/filtered.tmp.F.fa &>> $output_f/primermatch.log
        rm $output_f/filtered.tmp.F.fa
        echo -e "...done selecting reads.\n"
    
        echo -e "\nReverse-complementing reads...\n"
        $usearch -fastx_revcomp $output_f/filtered.fa -label_suffix _RC -fastaout $output_f/filtered.tmp.RC.fa &>> $output_f/primermatch.log
    
        echo -e "\nSelecting RC reads with primer matches...\n"
        $cutadapt --discard-untrimmed -g ${primerF} -O ${MIN_F} -e ${ERR_F} -m ${minlength} -o $output_f/filtered.tmp.RC.F.fa $output_f/filtered.tmp.RC.fa &>> $output_f/primermatch.log
        $cutadapt --discard-untrimmed -a ${primerR} -O ${MIN_R} -e ${ERR_R} -m ${minlength} -o $output_f/filtered_primermatch_RC.fa $output_f/filtered.tmp.RC.F.fa &>> $output_f/primermatch.log
        if [ -e $output_f/filtered.tmp.RC.F.fa ]; then rm $output_f/filtered.tmp.RC.F.fa; fi
        if [ -e $output_f/filtered.tmp.RC.fa ]; then rm $output_f/filtered.tmp.RC.fa; fi
    
        if [ -e $output_f/filtered_primermatch_RC.fa ]; then cat $output_f/filtered_primermatch_RC.fa >> $output_f/filtered_primermatch.fa; fi
        if [ -e $output_f/filtered_primermatch_RC.fa ]; then rm $output_f/filtered_primermatch_RC.fa; fi
        echo -e "\n...done selecting RC reads.\n"
        
        if [ ! -e $output_f/filtered_primermatch.fa ]; then echo -e "\n${output_f}/filtered_primermatch.fa was not created. Primer match failed. Exiting...\n"; exit 1; fi
        fi
fi

## DEREPLICATION
if [ -e $output_f/uniques.fa ]
then
    echo -e "\nDereplicated sequences already exist. Skip this step.\n"
else
    echo -e "\nDereplicating reads...\n"
    if [ -e $output_f/filtered_primermatch.fa ]
    then
        $usearch -fastx_uniques $output_f/filtered_primermatch.fa -minuniquesize ${minsize} -sizeout -relabel Uniq -fastaout $output_f/uniques.fa -threads ${threads} &> $output_f/dereplication.log
    else
        $usearch -fastx_uniques $output_f/filtered.fa -minuniquesize ${minsize} -sizeout -relabel Uniq -fastaout $output_f/uniques.fa -threads ${threads} &> $output_f/dereplication.log
    fi
    echo -e "\n...done dereplicating reads.\n"
fi

if [ ! -e $output_f/uniques.fa ]; then echo -e "\n${output_f}/uniques_uparse.fa was not created. Dereplication failed. Exiting...\n"; exit 1; fi


## START OF DEFINED COMMUNITY ROUTINE
if [[ -v ref ]]
then

## ALIGNMENT to reference sequences with USEARCH
if [ -e $output_f/initial_classification.txt ]
then
    echo -e "\nInitial classification file already exists. Skip this step.\n"
else
    echo -e "\nClassifying reads according to given reference sequences (USEARCH_GLOBAL $ref_id ID)...\n"
    $usearch -usearch_global $output_f/uniques.fa -db $ref -strand both -id $ref_id -top_hit_only -output_no_hits -blast6out $output_f/initial_classification.txt -threads ${threads} &> $output_f/initial_classification.log
    echo -e "\n...done classifying reads.\n"
fi

if [ ! -e $output_f/initial_classification.txt ]; then echo -e "\n${output_f}/initial_classification.txt was not created. Classification failed. Exiting...\n"; exit 1; fi

## QUANTIFICATION of reference sequences
if [ -e $output_f/otutab_initial_classified.txt ]
then
    echo -e "\nInitial classified OTU table already exists. Skip this step.\n"
else
    echo -e "\nQuantifying vs reference sequences...\n"
    $usearch -otutab $output_f/filtered.fa -otus $ref -strand both -id $ref_id -notmatched $output_f/unclassified.fa -otutabout $output_f/otutab_initial_classified.txt -biomout $output_f/otutab_initial_classified.json -mothur_shared_out $output_f/otutab_initial_classified.mothur -sample_delim . -threads ${threads} &> $output_f/make_otutab_initial_classified.log
    # APPEND unclassified counts CURRENTLY ONLY FOR .TXT FILE, OTHER FORMATS ARE NOT NICE
    echo -ne "Unclassified" >> $output_f/otutab_initial_classified.txt
    awk -F '.' '/^>/ {print $1}' $output_f/unclassified.fa | sort -V | uniq -c | awk '{printf "\t"$1}' >> $output_f/otutab_initial_classified.txt
    echo -e "\n...done quantifying reads.\n"
fi

if [ ! -e $output_f/otutab_initial_classified.txt ]; then echo -e "\n{$output_f}/otutab_initial_classified.txt was not created. Quantification failed. Exiting...\n"; exit 1; fi

## DEREPLICATION of unclassified reads
if [ -e $output_f/unclassified_uniques.fa ]
then
    echo -e "\nDereplicated unclassified sequences already exist. Skip this step.\n"
else
    echo -e "\nDereplicating unclassified reads...\n"
    $usearch -fastx_uniques $output_f/unclassified.fa -sizeout -relabel Uniq -fastaout $output_f/unclassified_uniques.fa -threads ${threads} &> $output_f/dereplication_unclassified.log
    echo -e "\n...done dereplicating reads.\n"
fi

if [ ! -e $output_f/unclassified_uniques.fa ]; then echo -e "\n${output_f}/unclassified_uniques.fa was not created. Dereplication failed. Exiting...\n"; exit 1; fi

## CLUSTERING unclassified reads with UPARSE
if [ -e $output_f/otus_unclassified.fa ]
then
    echo -e "\nClustered unclassified sequences already exist. Skip this step.\n"
else
    echo -e "\nClustering unclassified reads and de-novo chimera checking (UPARSE algorithm)...\n"
    $usearch -cluster_otus $output_f/unclassified_uniques.fa -minsize ${minsize} -otus $output_f/otus_unclassified.fa -relabel Unclass &> $output_f/clustering_unclassified.log
    cat $ref $output_f/otus_unclassified.fa > $output_f/final_references.fa
    echo -e "\n...done clustering reads.\n"
fi

if [ ! -e $output_f/otus_unclassified.fa ]; then echo -e "\n${output_f}/otus_unclassified.fa was not created. Clustering failed. Exiting...\n"; exit 1; fi

## QUANTIFICATION of unclassified otus
if [ -e $output_f/otutab_unclassified.txt ]
then
    echo -e "\nUnclassified OTU table already exists. Skip this step.\n"
else
    echo -e "\nQuantifying vs unclassified otus...\n"
    $usearch -otutab $output_f/unclassified.fa -otus $output_f/otus_unclassified.fa -strand both -id $ref_id -otutabout $output_f/otutab_unclassified.txt -biomout $output_f/otutab_unclassified.json -mothur_shared_out $output_f/otutab_unclassified.mothur -sample_delim . -threads ${threads} &> $output_f/make_otutab_unclassified.log
    echo -e "\n...done quantifying reads.\n"
fi

if [ ! -e $output_f/otutab_unclassified.txt ]; then echo -e "\n{$output_f}/otutab_unclassified.txt was not created. Quantification failed. Exiting...\n"; exit 1; fi

## REALIGNMENT to reference sequences plus unclassified otus
if [ -e $output_f/final_classification.txt ]
then
    echo -e "\nFinal classification file already exists. Skip this step.\n"
else
    echo -e "\nClassifying reads according to given reference sequences and unclassified otus (USEARCH_GLOBAL $ref_id ID)\n"
    $usearch -usearch_global $output_f/uniques.fa -db $output_f/final_references.fa -strand both -id $ref_id -top_hit_only -output_no_hits -blast6out $output_f/final_classification.txt -threads ${threads} &> $output_f/final_classification.log
    echo -e "\n...done classifying reads.\n"
fi

if [ ! -e $output_f/final_classification.txt ]; then echo -e "\n{$output_f}/final_classification.txt was not created. Classification failed. Exiting...\n"; exit 1; fi

## QUANTIFICATION of reference sequences plus unclassified otus
if [ -e $output_f/otutab_final_classified.txt ]
then
    echo -e "\nFinal classified OTU table already exists. Skip this step.\n"
else
    echo -e "\nQuantifying vs reference sequences...\n"
    $usearch -otutab $output_f/filtered.fa -otus $output_f/final_references.fa -strand both -id $ref_id -otutabout $output_f/otutab_final_classified.txt -biomout $output_f/otutab_final_classified.json -mothur_shared_out $output_f/otutab_final_classified.mothur -sample_delim . -threads ${threads} &> $output_f/make_otutab_final_classified.log
    echo -e "\n...done quantifying reads.\n"
fi

if [ ! -e $output_f/otutab_final_classified.txt ]; then echo -e "\n{$output_f}/otutab_final_classified.txt was not created. Quantification failed. Exiting...\n"; exit 1; fi

## Last-common-ancestor function
function lca(){ cat $@ | sed -e '$!{N;s/^\(.*\).*\n\1.*$/\1\n\1/;D;}' | awk -F ";" '{$NF=""; OFS=";"; print $0}'; return; }

# TAXONOMY of unclassified otus
if [ -e $output_f/taxsearch_unclassified.tax ]
then
    echo -e "\nTaxonomy search file for unclassified OTUS already exist. Skip this step.\n"
else
    if [[ -z ${db+x} ]]
    then
        echo -e "\nTaxonomical database not provided. Skipping taxonomy assignment.\n"
    else
        echo -e "\nSearching unclassified OTUs...\n"
        $usearch -usearch_global $output_f/otus_unclassified.fa -db ${db} -id ${tax_id} -maxaccepts 500 -maxrejects 500 -strand both -top_hits_only -output_no_hits -blast6out $output_f/taxsearch_unclassified.tax -threads ${threads} &> $output_f/taxsearch_unclassified.log
        echo -e "\n...done annotating OTUs.\n"

        if [ ! -e $output_f/taxsearch_unclassified.tax ]; then echo -e "\n${output_f}/taxsearch_unclassified.tax was not created. Taxonomy search for unclassified OTUS failed. Exiting...\n"; exit 1; fi
    fi
fi

# LCA of unclassified otus
if [ -e $output_f/taxonomy_unclassified_lca.txt ]
then
    echo -e "\nTaxonomy assignment file for unclassified OTUS already exist. Skip this step.\n"
else
    if [[ -z ${db+x} ]]
    then
        echo -e "\nTaxonomical database not provided. Skipping taxonomy assignment.\n"
    else
        echo -e "\nAnnotating unclassified OTUs with LCA...\n"
        for i in $(cut -f 1 -d $'\t' $output_f/taxsearch_unclassified.tax | sort | uniq); do id=$(grep -m 1 -P $i'\t' $output_f/taxsearch_unclassified.tax | cut -f 3 -d$'\t'); res=$(grep -P $i'\t' $output_f/taxsearch_unclassified.tax | cut -f 2 -d$'\t' | cut -f 1 -d ' ' --complement | lca); echo -e $i'\t'$res'\t'$id; done > $output_f/taxonomy_unclassified_lca.txt
        echo -e "\n...done annotating OTUs.\n"

        if [ ! -e $output_f/taxonomy_unclassified_lca.txt ]; then echo -e "\n${output_f}/taxonomy_unclassified_lca.txt was not created. Taxonomy assignment for unclassified OTUS failed. Exiting...\n"; exit 1; fi
    fi
fi

## END OF DEFINED COMMUNITY ROUTINE
else

## CLUSTERING with UPARSE
if [ -e $output_f/otus_uparse.fa ]
then
    echo -e "\nClustered (UPARSE) sequences already exist. Skip this step.\n"
else
    echo -e "\nClustering reads and de-novo chimera checking (UPARSE algorithm)...\n"
    $usearch -cluster_otus $output_f/uniques.fa -minsize ${minsize} -otus $output_f/otus_uparse.fa -relabel Otu &> $output_f/clustering.log
    echo -e "\n...done clustering reads.\n"
fi

if [ ! -e $output_f/otus_uparse.fa ]; then echo -e "\n${output_f}/otus.fa was not created. Clustering failed. Exiting...\n"; exit 1; fi

## DENOISING with UNOISE
if [ -e $output_f/otus_unoise.fa ]
then
    echo -e "\nDenoised sequences already exist. Skip this step.\n"
else
    echo -e "\nDenoising reads and de-novo chimera checking (UNOISE3 algorithm)...\n"
    $usearch -unoise3 $output_f/uniques.fa -zotus $output_f/otus_unoise.fa -relabel Otu &> $output_f/denoising.log
    sed -i 's/Zotu/Otu/g' $output_f/otus_unoise.fa
    echo -e "\n...done denoising reads.\n"
fi

if [ ! -e $output_f/otus_unoise.fa ]; then echo -e "\n${output_f}/otus_unoise.fa was not created. Denoising failed. Exiting...\n"; exit 1; fi

## TAXONOMICAL ANNOTATION
# Last-common-ancestor function
function lca(){ cat $@ | sed -e '$!{N;s/^\(.*\).*\n\1.*$/\1\n\1/;D;}' | awk -F ";" '{$NF=""; OFS=";"; print $0}'; return; }

# Uparse OTUs
if [ -e $output_f/taxsearch_uparse.tax ]
then
    echo -e "\nTaxonomy search file for OTUS (UPARSE algorithm) already exist. Skip this step.\n"
else
    if [[ -z ${db+x} ]]
    then
        echo -e "\nTaxonomical database not provided. Skipping taxonomy assignment.\n"
    else
        echo -e "\nSearching OTUs (UPARSE algorithm)...\n"
        $usearch -usearch_global $output_f/otus_uparse.fa -db ${db} -id ${tax_id} -maxaccepts 500 -maxrejects 500 -strand both -top_hits_only -output_no_hits -blast6out $output_f/taxsearch_uparse.tax -threads ${threads} &> $output_f/taxsearch_uparse.log
        echo -e "\n...done annotating OTUs.\n"

        if [ ! -e $output_f/taxsearch_uparse.tax ]; then echo -e "\n${output_f}/taxsearch_uparse.tax was not created. Taxonomy search for OTUS (UPARSE algorithm) failed. Exiting...\n"; exit 1; fi
    fi
fi

# LCA for UPARSE OTUs
if [ -e $output_f/taxonomy_uparse_lca.txt ]
then
    echo -e "\nTaxonomy assignment file for OTUS (UPARSE algorithm) already exist. Skip this step.\n"
else
    if [[ -z ${db+x} ]]
    then
        echo -e "\nTaxonomical database not provided. Skipping taxonomy assignment.\n"
    else
        echo -e "\nAnnotating OTUs (UPARSE algorithm) with LCA...\n"
        for i in $(cut -f 1 -d $'\t' $output_f/taxsearch_uparse.tax | sort | uniq); do id=$(grep -m 1 -P $i'\t' $output_f/taxsearch_uparse.tax | cut -f 3 -d$'\t'); res=$(grep -P $i'\t' $output_f/taxsearch_uparse.tax | cut -f 2 -d$'\t' | cut -f 1 -d ' ' --complement | lca); echo -e $i'\t'$res'\t'$id; done > $output_f/taxonomy_uparse_lca.txt
        echo -e "\n...done annotating OTUs.\n"

        if [ ! -e $output_f/taxonomy_uparse_lca.txt ]; then echo -e "\n${output_f}/taxonomy_uparse_lca.txt was not created. Taxonomy assignment for OTUS (UPARSE algorithm) failed. Exiting...\n"; exit 1; fi
    fi
fi

# Unoise OTUs
if [ -e $output_f/taxsearch_unoise.tax ]
then
    echo -e "\nTaxonomy search file for OTUS (UNOISE3 algorithm) already exist. Skip this step.\n"
else
    if [[ -z ${db+x} ]]
    then
        echo -e "\nTaxonomical database not provided. Skipping taxonomy assignment.\n"
    else
        echo -e "\nSearching OTUs (UNOISE3 algorithm)...\n"
        $usearch -usearch_global $output_f/otus_unoise.fa -db ${db} -id ${tax_id} -maxaccepts 500 -maxrejects 500 -strand both -top_hits_only -output_no_hits -blast6out $output_f/taxsearch_unoise.tax -threads ${threads} &> $output_f/taxsearch_unoise.log
        echo -e "\n...done annotating OTUs.\n"
        
        if [ ! -e $output_f/taxsearch_unoise.tax ]; then echo -e "\n${output_f}/taxsearch_unoise.tax was not created. Taxonomy search for OTUS (UNOISE3 algorithm) failed. Exiting...\n"; exit 1; fi
    fi
fi

# LCA for UNOISE OTUs
if [ -e $output_f/taxonomy_unoise_lca.txt ]
then
    echo -e "\nTaxonomy assignment file for OTUS (UNOISE3 algorithm) already exist. Skip this step.\n"
else
    if [[ -z ${db+x} ]]
    then
        echo -e "\nTaxonomical database not provided. Skipping taxonomy assignment.\n"
    else
        echo -e "\nAnnotating OTUs (UNOISE3 algorithm) with LCA...\n"
        for i in $(cut -f 1 -d $'\t' $output_f/taxsearch_unoise.tax | sort | uniq); do id=$(grep -m 1 -P $i'\t' $output_f/taxsearch_unoise.tax | cut -f 3 -d$'\t'); res=$(grep -P $i'\t' $output_f/taxsearch_unoise.tax | cut -f 2 -d$'\t' | cut -f 1 -d ' ' --complement | lca); echo -e $i'\t'$res'\t'$id; done > $output_f/taxonomy_unoise_lca.txt
        echo -e "\n...done annotating OTUs.\n"
        
        if [ ! -e $output_f/taxonomy_unoise_lca.txt ]; then echo -e "\n${output_f}/taxonomy_unoise_lca.txt was not created. Taxonomy assignment for OTUS (UNOISE3 algorithm) failed. Exiting...\n"; exit 1; fi
    fi
fi

## QUATIFICATION OF OTUs

# Uparse OTUs
if [ -e $output_f/otutab_uparse.txt ]
then
    echo -e "\nOTU tables (UPARSE algorithm) already exist. Skip this step.\n"
else
    echo -e "\nQuantifying vs OTUs (UPARSE algorithm) using all filtered reads...\n"
    $usearch -otutab $output_f/filtered.fa -otus $output_f/otus_uparse.fa -strand both -id 0.97 -otutabout $output_f/otutab_uparse.txt -biomout $output_f/otutab_uparse.json -mothur_shared_out $output_f/otutab_uparse.mothur -sample_delim . -threads ${threads} &> $output_f/make_otutab_uparse.log
    echo -e "\n...done quantifying vs OTUs using al reads.\n"
fi

if [ ! -e $output_f/otutab_uparse.txt ]; then echo -e "$\n{output_f}/otutab_uparse.txt was not created. Quantification of OTUs (UPARSE algorithm) failed. Exiting...\n"; exit 1; fi

# Unoise OTUs
if [ -e $output_f/otutab_unoise.txt ]
then
    echo -e "\nOTU tables (UNOISE3 algorithm) already exist. Skip this step.\n"
else
    echo -e "\nQuantifying vs OTUs (UNOISE3 algorithm) using all filtered reads...\n"
    $usearch -otutab $output_f/filtered.fa -zotus $output_f/otus_unoise.fa -strand both -id 0.97 -otutabout $output_f/otutab_unoise.txt -biomout $output_f/otutab_unoise.json -mothur_shared_out $output_f/otutab_unoise.mothur -sample_delim . -threads ${threads} &> $output_f/make_otutab_unoise.log
    echo -e "\n...done quantifying vs OTUs using al reads.\n"
fi

if [ ! -e $output_f/otutab_unoise.txt ]; then echo -e "$\n{output_f}/otutab_unoise.txt was not created. Quantification of OTUs (UNOISE3 algorithm) failed. Exiting...\n"; exit 1; fi

fi

function fasta_length_hist(){
cat $1 | awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}' | awk -F '\t' '{print $2}' | awk '{print length($1)}' | sort -n | uniq -c | awk '{print $2 " " $1}'
return
}

## REPORT for defined community run
function dcreport
{
cat <<EOF

-- Results --

Number of R1 reads:                         $(find $input_f/*R1*.fastq -exec wc -l {} \; | awk '{total += $1} END {print total}' | awk '{print $1/4}')
Number of R2 reads:                         $(find $input_f/*R2*.fastq -exec wc -l {} \; | awk '{total += $1} END {print total}' | awk '{print $1/4}')
Number of merged reads:                     $(wc -l $output_f/merged.fq | awk '{print $1/4}')
Number of quality-filtered reads:           $(grep "^>" -c $output_f/filtered.fa)
Number of primer-matched reads:             $(grep "^>" -c $output_f/filtered_primermatch.fa)
Number of dereplicated sequences:           $(grep "^>" -c $output_f/uniques.fa)
Number of reads mapping to references:      $(grep "mapped to OTUs" $output_f/make_otutab_initial_classified.log | awk '{print $1}')
Number of unclassified OTUs:                $(grep "^>" -c $output_f/otus_unclassified.fa)
Number of reads mapping to both:            $(grep "mapped to OTUs" $output_f/make_otutab_final_classified.log | awk '{print $1}')

Length distribution of dereplicated sequences:
$(fasta_length_hist $output_f/uniques.fa)

Length distribution of dereplicated unclassified sequences:
$(fasta_length_hist $output_f/unclassified_uniques.fa)

Length distribution of unclassified OTUs:
$(fasta_length_hist $output_f/otus_unclassified.fa)

EOF
}

## REPORT for normal run
function report
{
cat <<EOF

-- Results --

Number of R1 reads:                         $(find $input_f/*R1*.fastq -exec wc -l {} \; | awk '{total += $1} END {print total}' | awk '{print $1/4}')
Number of R2 reads:                         $(find $input_f/*R2*.fastq -exec wc -l {} \; | awk '{total += $1} END {print total}' | awk '{print $1/4}')
Number of merged reads:                     $(wc -l $output_f/merged.fq | awk '{print $1/4}')
Number of quality-filtered reads:           $(grep "^>" -c $output_f/filtered.fa)
Number of primer-matched reads:             $(grep "^>" -c $output_f/filtered_primermatch.fa)
Number of dereplicated sequences:           $(grep "^>" -c $output_f/uniques.fa)
Number of OTUs (UPARSE):                    $(grep "^>" -c $output_f/otus_uparse.fa)
Number of OTUs (UNOISE3):                   $(grep "^>" -c $output_f/otus_unoise.fa)
Number of reads mapping to OTUs (UPARSE):   $(grep "mapped to OTUs" $output_f/make_otutab_uparse.log | awk '{print $1}')
Number of reads mapping to OTUs (UNOISE3):  $(grep "mapped to OTUs" $output_f/make_otutab_unoise.log | awk '{print $1}')
    
Length distribution of OTUs (UPARSE):
$(fasta_length_hist $output_f/otus_uparse.fa)

Length distribution of OTUs (UNOISE3):
$(fasta_length_hist $output_f/otus_unoise.fa)

Length distribution of dereplicated sequences:
$(fasta_length_hist $output_f/uniques.fa)

EOF
}

if [ -e $output_f/report.txt ]
then
    echo -e "\nReport already exist. Skip this step.\n"
else
    echo -e "\nMaking report...\n"
    show_params > $output_f/report.txt
    if [[ -v ref ]]
    then
        dcreport >> $output_f/report.txt
    else
        report >> $output_f/report.txt
    fi
    echo -e "\n...Report done\n"

fi

if [ ! -e $output_f/report.txt ]; then echo -e "$\n{output_f}/report.txt was not created. The report failed. Exiting...\n"; exit 1; fi

echo -e "\nPipeline successfully finished\n"

