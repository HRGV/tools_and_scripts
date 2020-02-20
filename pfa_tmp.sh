#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#change to reasonable amount of CPUs - small assemblies (bacteria) 12, normal assemblies (100 Mb - 1 Gb) 24, super big jobs (>1Gb) 64. spades will run out of memory for assemblies >4Gb even on highmem nodes
#$ -V
#$ -q main.q
# change queue to main.q@@himem if you need more than 400 gb of RAM - usually necessary from ~300 Mb total assembly size.
#
echo "job started: "
echo "job ID: $JOB_ID"
date
hostname
echo "shell: $SHELL"
####@@himem - add this to main.q if you need more than 400 gb

####general variables
project=YOUR-PROJET # this is a folder to find the data in, should be in your user space
lib=YOUR-Library-NAME;

####amount of reads - set to -1 to run with all reads, 100k is a good number for testing
reads=100000;

####create scratch folder
mkdir /scratch/$USER/tmp.$JOB_ID -p; 

####sync files
rsync -a /YOURFOLDER/$USER/$project/$lib/ /scratch/$USER/tmp.$JOB_ID; #change dir according to project

####move to scratch for computations - scratch is a local harddiks on the compute node.
cd /scratch/$USER/tmp.$JOB_ID/

####this assumes you are working from reads in a format XXX_f.fq.gz and XXX_r.fq.gz
for n in *_f*.gz; do #this picks up the Forward read file name that is assumed to be in the folder you have just synced to scratch
   ####run phyloFlash on raw reads that should be in the base folder of the library
	phyloFlash.pl -lib ${lib}_pf -read1 $n -read2 ${n%_f.fq.gz}_r.fq.gz -everything -CPUs 12 -readlength 150 -readlimit=$reads;
   ####qualtiy and adapter trim reads into a folder ../library/trimmed
	mkdir  /scratch/$USER/tmp.$JOB_ID/trimmed;
	cd /scratch/$USER/tmp.$JOB_ID/trimmed;
	bbduk.sh ref=/opt/extern/bremen/symbiosis/tools_HGV/adapters.fa ktrim=l minlength=36 mink=11 hdist=1 in=../$n in2=../${n%_f.fq.gz}_r.fq.gz out=${lib}_ktriml.fq.gz reads=$reads;
	bbduk.sh ref=/opt/extern/bremen/symbiosis/tools_HGV/adapters.fa ktrim=r trimq=2 qtrim=rl minlength=36 mink=11 hdist=1 in=${lib}_ktriml.fq.gz interleaved=t out=${lib}_q2_ktrimmed.fq.gz;
   ####kmer filter reads into a folder ../library/filtered
	higher=3; #define cutoff here 
	mkdir  /scratch/$USER/tmp.$JOB_ID/filtered;
	cd /scratch/$USER/tmp.$JOB_ID/filtered;
	# - you can optionally find the cutoff for your particular library using bbnorms kmer histogram analysis support
	# bbnorm.sh -Xmx400g in=../trimmed/${lib}_q2_ktrimmed.fq.gz hist=${lib}_all_khist.txt peaks=${lib}_all_khist_peaks.txt threads=12 passes=1 interleaved=t;
	bbnorm.sh -Xmx400g in=../trimmed/${lib}_q2_ktrimmed.fq.gz lowbindepth=1 highbindepth=$higher outhigh=${lib}_q2_ktrimmed_k31higher${higher}.fq.gz passes=1 threads=24 interleaved=t; #cutoff is the highbindepth
   ####assemble with megahit into a folder ../library/assemblies
	as=${lib}_all_kfilt${higher}_MH_21_141; #rename for the parameters you are using
	mkdir /scratch/$USER/tmp.$JOB_ID/assemblies -p;
	cd /scratch/$USER/tmp.$JOB_ID/assemblies;
	megahit --k-min 21 --k-max 141 --k-step 20 -m 0.4 -t 24 --out-prefix $as --12 ../filtered/${lib}_q2_ktrimmed_k31higher${higher}.fq.gz -o $as;
	megahit_toolkit contig2fastg 141 ./${as}/intermediate_contigs/k141.contigs.fa > ./${as}/${as}.fastg;
   ####fastgfish megahit assembly into ../library/fastgfish
	mkdir /scratch/$USER/tmp.$JOB_ID/fastgfish -p;
	cd /scratch/$USER/tmp.$JOB_ID/fastgfish;
	ln -sf ../assemblies/${as}/intermediate_contigs/k141.contigs.fa ./${as}.k141.contigs.fasta;
	ln -sf ../assemblies/${as}/${as}.fastg .;
	phyloFlash_fastgFishing.pl --fasta ${as}.k141.contigs.fasta --fastg ${as}.fastg --out ${as}_fastgbin --compare-zip ../${lib}_pf.phyloFlash.tar.gz --assembler megahit --outfasta Yes --min-SSU-frac 0.6
   ####bin using metabat into ../library/metabat
	mkdir /scratch/$USER/tmp.$JOB_ID/metabat -p;
	cd /scratch/$USER/tmp.$JOB_ID/metabat;
	bbmap.sh in=../filtered/${lib}_q2_ktrimmed_k31higher${higher}.fq.gz ref=../assemblies/${as}/${as}.contigs.fa outm=${lib}_on_${as}.bam outputunmapped=f nodisk=t unpigz=t pigz=t threads=24;
	samtools sort --threads 24 outm=${lib}_on_${as}.bam -o ${lib}_on_${as}.sorted.bam ;
	rm ${lib}_on_${as}.bam;
	jgi_summarize_bam_contig_depths --outputDepth ${lib}_on_${as}.depth.txt ${lib}_on_${as}.sorted.bam;
	rm ${lib}_on_${as}.sorted.bam;
  	metabat2 -i ../assemblies/${as_sp}/scaffolds.fasta -a ${lib}_on_${as_sp}.depth.txt -o ./bin_${as} -t 24	
   ####assemble with spades into a folder ../library/assemblies
	cd /scratch/$USER/tmp.$JOB_ID/assemblies;
	as_sp=${lib}_all_kfilt${higher}_SPm_21_127; #rename 
	spades.py -k 21,33,55,77,99,127 -m 400 -t 24 -o $as_sp --12 ../filtered/${lib}_q2_ktrimmed_k31higher${higher}.fq.gz; #use -meta switch for communities with more than 1 bug!
   ####fastgfish spades assembly into ../library/fastgfish
	cd /scratch/$USER/tmp.$JOB_ID/fastgfish;
	ln -sf ../assemblies/${as_sp}/scaffolds.fasta ./${as_sp}.scaffolds.fasta;
	ln -sf ../assemblies/${as_sp}/assembly_graph.fastg ./${as_sp}.assembly_graph.fastg;
	ln -sf ../assemblies/${as_sp}/scaffolds.paths ./${as_sp}.scaffolds.paths;
	phyloFlash_fastgFishing.pl --fasta ${as_sp}.scaffolds.fasta --fastg ${as_sp}.assembly_graph.fastg --paths ${as_sp}.scaffolds.paths --out ${as_sp}_fastgbin --compare-zip ../${lib}_pf.phyloFlash.tar.gz --assembler spades --outfasta Yes --min-SSU-frac 0.6
   ####bin using metabat into ../library/metabat
	mkdir /scratch/$USER/tmp.$JOB_ID/metabat -p;
	cd /scratch/$USER/tmp.$JOB_ID/metabat;
	bbmap.sh in=../filtered/${lib}_q2_ktrimmed_k31higher${higher}.fq.gz ref=../assemblies/${as_sp}/scaffolds.fasta outm=${lib}_on_${as_sp}.bam outputunmapped=f nodisk=t unpigz=t pigz=t threads=24;
	samtools sort --threads 24 ${lib}_on_${as_sp}.bam -o ${lib}_on_${as_sp}.sorted.bam ;
	rm ${lib}_on_${as_sp}.bam;
	jgi_summarize_bam_contig_depths --outputDepth ${lib}_on_${as_sp}.depth.txt ${lib}_on_${as_sp}.sorted.bam;
	rm ${lib}_on_${as_sp}.sorted.bam;
  	metabat2 -i ../assemblies/${as_sp}/scaffolds.fasta -a ${lib}_on_${as_sp}.depth.txt -o ./bin_${as_sp} -t 24
#final words
done;

####clean up the data structure before syncing back (this is explicit for each folder to be able to deselect)
cd /scratch/$USER/tmp.$JOB_ID/

####processing pipeline files
find . -name 'trimmed' -type d -exec rm {} -R \;
find . -name 'filtered' -type d -exec rm {} -R \;

####spades temporary files
find . -name 'split_input' -type d -exec rm {} -R \;
find . -name 'K21' -type d -exec rm {} -R \;
find . -name 'K33' -type d -exec rm {} -R \;
find . -name 'K55' -type d -exec rm {} -R \;
find . -name 'K77' -type d -exec rm {} -R \;
find . -name 'K99' -type d -exec rm {} -R \;
find . -name 'K127' -type d -exec rm {} -R \;
find . -name 'misc' -type d -exec rm {} -R \;
find . -name 'tmp' -type d -exec rm {} -R \;

###rsync back
rsync -a /scratch/grmber/tmp.$JOB_ID/ /opt/extern/bremen/symbiosis/$USER/$project/$lib/; #change dir according to project
if [ "$?" -eq "0" ] #checks status of rsync - $? gives the error code of the last command, and if it worked it will clean up, otherwise it will not
	then
		rm /scratch/$USER/tmp.$JOB_ID -R;
		echo "Rsync done"
	else
		echo "Error while running rsync"
fi
echo "job finished: "
date

