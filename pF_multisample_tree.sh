#!/bin/bash
# Script to collect the SPADES assembled rRNAS and references in a common tree
# Script by Harald Gruber-Vodicka, 2025 with help from ChatGPT
# Thanks to N.Leisch for pushing this again and again
# This is work in progress and might not work as intended
#
# The script expects a folder with .tar.gz compressed phyloFlash 3.3 and above output files.
# Easiest way to generate the necessary input is phyloFlash.pl [...] with the -almosteverything switch

# Dependencies:
# - phyloFlash output in .tar.gz format
# - vsearch         (for clustering)
# - mafft           (for alignment)
# - reformat.sh     (part of bbsuite, used for deduplication and renaming)
# - fasttree        (for tree calculation)
# - tar             (for extracting files)

# Expected command-line options: 
#   -p PROJECTNAME (default: test)
#   -c CPUs        (default: 4)

#------------------------------------------------
# Functions for usage and exit


#subs
#post usage
usage() {                                      # Function: Print a help message.
  echo "\nUsage: $0 [ -p PROJECTNAME ] [ -c CPUs ]" 1>&2 
}
exit_abnormal() {                              # Function: Exit with error.
  usage
  exit 1
}
#variables
cpus=4
project=test

#------------------------------------------------
# Parse command line options
# The options are the project name and the number of CPUs to be used.
while getopts ":p:c:" opts; do              # Loop: Get the next option;
                                               # use silent error checking;
                                               # options p and c take arguments.
  case $opts in                         # 
    p)                                         # If the option is c,
      echo "argument -p called with parameter $OPTARG"
      project=${OPTARG}                        # set $project to specified value.
      	;;
    c)                                         # If the option is c,
       echo "argument -c called with parameter $OPTARG"
       cpus=${OPTARG}                          # Set $cpus to specified value.
        ;;
    *)                                         # If expected argument omitted:
      echo "invalid command"
      exit_abnormal                            # Exit abnormally.
      ;;
  esac
done

echo "\n Detected commandline parameters"
echo "Using $cpus CPUs to get the job done"
echo "Will use $project as prefix for all files"
echo "\n"

#------------------------------------------------
# MAIN
# Check for required dependencies
for cmd in vsearch mafft fasttree reformat.sh tar; do
  command -v $cmd >/dev/null 2>&1 || { echo >&2 "$cmd is required but not installed. Aborting."; exit 1; }
done

#------------------------------------------------
# report the number of found files to process
# or abort if no tar.gz files are found

file_count=$(ls *.tar.gz 2>/dev/null | wc -l)
echo "Found $file_count .tar.gz files for processing"

shopt -s nullglob
tar_files=(*.tar.gz)
if [ ${#tar_files[@]} -eq 0 ]; then
  echo "No .tar.gz files found in the current directory. Exiting."
  exit 1
fi


#------------------------------------------------
# Extract assembled SSU rRNA gene sequences and related SILVA dbhits
# and concatenate PF hits and dbhits into files in the project folder

echo "Extracting assembled SSUs and dbhits into project folder\n"
mkdir -p "${project}_data"

for file in *tar.gz; do  
  	sample=$(basename "$file" .tar.gz)
 	echo "Processing $sample\n"
	echo "collecting assembled SSUs\n"
	# Corrected tar command (replace `-0` with `-O`)
	tar -xzf "$file" -C "${project}_data/" "${file%phy*}spades_rRNAs.final.fasta" 2>/dev/null
	echo "Extracting reference SSUs\n"
 	tar -xzf "$file" -C "${project}_data/" "${file%phy*}all.dbhits.NR97.fa" 2>/dev/null
done

##concatenate the files
cat ${project}_data/*.spades_rRNAs.final.fasta > "${project}_data/${project}.collection_spades.fasta"
cat ${project}_data/*.all.dbhits.NR97.fa > "${project}_data/${project}.collection_references.fasta"


#------------------------------------------------
# Remove duplicate sequences from the references using vsearch
echo "Creating a nonredundant reference SSUs database\n"
vsearch --threads $cpus --cluster_fast "${project}_data/$project.collection_references.fasta" --id 0.97 --relabel_keep --notrunclabels --centroids "${project}_data/$project.collection_references_NR97.fasta"

#------------------------------------------------
# Concatenate PF hits and nonredundant references
echo "Creating a NR-Ref and assembled SSUs database\n"
cat "${project}_data/$project.collection_spades.fasta" "${project}_data/$project.collection_references_NR97.fasta" > "${project}_data/$project.spades_and_refsNR97_collection.fasta"

#------------------------------------------------
# Check for duplicated sequence headers and rename if necessary
echo "Checking for duplicated sequence headers\n"
reformat.sh uniquenames=t in="${project}_data/$project.spades_and_refsNR97_collection.fasta" out="${project}_data/$project.spades_and_refsNR97_collection_uniqed.fasta"

#------------------------------------------------
# Calculate multiple sequence alignment using Mafft-Ginsi 
echo "Calculating alignment using Mafft-Ginsi\n"
mafft-ginsi --thread $cpus --adjustdirection "${project}_data/$project.spades_and_refsNR97_collection_uniqed.fasta" > "${project}_data/$project.spades_and_refsNR97_collection.mafft_ginsi.fasta"

#------------------------------------------------
# Build the phylogenetic tree using FastTree with GTR model
echo "Calculating tree using Fasttree\n"
fasttree -gtr -nt -quote "${project}_data/$project.spades_and_refsNR97_collection.mafft_ginsi.fasta" > "${project}_data/$project.spades_and_refsNR97_collection.mafft_ginsi.fasttree_GTR.tre"

#------------------------------------------------
# Final messages
echo "\nJob completed at: $(date)"
echo "\nPlease find your final outputs here:"
echo " - Dereplicated references: ${project}_data/$project.collection_references_NR97.fasta"
echo " - Combined & uniqued sequences: ${project}_data/$project.spades_and_refsNR97_collection_uniqed.fasta"
echo " - MAFFT alignment: ${project}_data/$project.spades_and_refsNR97_collection.mafft_ginsi.fasta"
echo " - Phylogenetic tree: ${project}_data/$project.spades_and_refsNR97_collection.mafft_ginsi.fasttree_GTR.tre"
echo "\nJob completed at: $(date)\n"

exit 0
