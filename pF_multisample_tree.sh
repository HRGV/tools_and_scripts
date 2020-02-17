#!/bin/bash
#script to collect the SPADES assembled rRNAS and references in a common tree
#script by Harald Gruber-Vodicka, 2020
#thanks to N.Leisch for pushing this
#this is work in progress and might not work as intended
#
#the script expects a folder with .tar.gz compressed phyloFlash 3.3 and above output files.
#easiest way to generate the necessary input is phyloFlash.pl [...] with the -almosteverything switch

#dependencies
#phyloFlash output in .tar.gz format
#vsearch - clustering
#mafft - alignment
#fasttree - tree calculation

#expected commandline options are -p PROJECTNAME [defaut: pf_collection] and -c CPUs for computing [default:1]


#post usage
usage() {                                      # Function: Print a help message.
  echo -e "\nUsage: $0 [ -p PROJECTNAME ] [ -c CPUs ]" 1>&2 
}
exit_abnormal() {                              # Function: Exit with error.
  usage
  exit 1
}

#get options from command line
#the options are the project name and the number of CPUs to be used.
while getopts ":p:c:" options; do              # Loop: Get the next option;
                                               # use silent error checking;
                                               # options n and t take arguments.
  case "${options}" in                         # 
    p)                                         # If the option is c,
      project=${OPTARG}                        # set $project to specified value.
      re_isaword='^[a-zA-Z_]+$'                 # Regex: match whole words only, a-z and AZ as well as - and _ allowed
      if ! [[ $project =~ $re_isaword ]] ; then   # if $TIMES not a whole number:
       echo -e "\n\e[91mError\e[0m: project must be a single word, only a-z, A-Z and _ allowed."
      	exit_abnormal
        exit 1
      fi
	;;
    c)                                         # If the option is c,
      cpus=${OPTARG}                          # Set $cpus to specified value.
      re_isanum='^[0-9]+$'                     # Regex: match whole numbers only
      if ! [[ $cpus =~ $re_isanum ]] ; then   # if $TIMES not a whole number:
        echo -e "\n\e[91mError\e[0m: CPUs must be a positive, whole number."
        exit_abnormal
        exit 1
      elif [ $cpus -eq "0" ]; then            # If it's zero:
        echo -e "\n\e[91mError\e[0m: CPUs must be greater than zero."
        exit_abnormal                          # Exit abnormally.
      fi
      ;;
    :)                                         # If expected argument omitted:
      echo -e "\n\e[91mError\e[0m: -${OPTARG} requires an argument."
      exit_abnormal                            # Exit abnormally.
      ;;
    *)                                         # If unknown (any other) option:
      exit_abnormal                            # Exit abnormally.
      ;;
  esac
done

echo -e "\n Detected commandline parameters"
echo -e "Using $cpus CPUs to get the job done"
echo -e "Will use $project as prefix for all files"
echo -e "\n"

#MAIN
#Extract all PF spades assembled SSU rRNA gene sequencens and the related SILVA dbhits (and possibly also emirge hits)
echo -e "Extracting assembled SSUs\n"
for i in *tar.gz; 
	do  tar -xf $i ${i%phy*}spades_rRNAs.final.fasta; 
done
echo -e "Extracting reference SSUs\n"
for i in *tar.gz; 
	do  tar -xf $i ${i%phy*}all.dbhits.NR97.fa; 
done

#concatenate PF hits and dbhits
echo -e "collecting all SSUs\n"
cat *.spades_rRNAs.final.fasta > $project.collection_spades.fasta
cat *.all.dbhits.NR97.fa > $project.collection_references.fasta

#run the dbits through vsearch to remove duplicates 
echo -e "Creating a nonredundant reference SSUs database\n"
vsearch --threads $cpus --cluster_fast $project.collection_references.fasta --id 0.97 --relabel_keep --notrunclabels --centroids $project.collection_references_NR97.fasta

#concatenate PF hits and references
echo -e "Creating a NR-Ref and assembled SSUs database\n"
cat $project.collection_spades.fasta $project.collection_references_NR97.fasta > $project.spades_and_refsNR97_collection.fasta

#check for duplicated seqs and rename
echo -e "Checking for duplicated sequence headers\n"
reformat.sh uniquenames=t in=$project.spades_and_refsNR97_collection.fasta out=$project.spades_and_refsNR97_collection_uniqed.fasta

#calculate alignment 
echo -e "Calculating alignment using Mafft-Ginsi\n"
mafft-ginsi --thread $cpus --adjustdirection $project.spades_and_refsNR97_collection_uniqed.fasta > $project.spades_and_refsNR97_collection.mafft_ginsi.fasta

#calculate tree
echo -e "Calculating tree using Fasttree\n"
fasttree -gtr -nt -quote $project.spades_and_refsNR97_collection.mafft_ginsi.fasta > $project.spades_and_refsNR97_collection.mafft_ginsi.fasttree_GTR.tre

#Done
echo -e "\n"
echo -e "All jobs completed - please find your files named $project.FILES\n"
echo -e "The tree is in a file named $project.spades_and_refsNR97_collection.mafft_ginsi.fasttree_GTR.tre\n"
exit 0  
