Fierst Lab protocol for ONT assembly and annotation of nematode genomes.

## Contents

 - [**1. Wet lab protocols**](#part-1-Wet-lab-protocols)
	- [1.1 Nematode culture](#11-Nematode-culture)
	- [1.2 High Molecular Weight DNA extraction](#12-High-Molecular-Weight-DNA-extraction)
	-   [1.2.1 Phenol Chloroform](#121-Phenol-Chloroform)
	-   [1.2.2 Circulomics Nanobind](#122-Circulomics-Nanobind)
	- [1.3 Short Read Eliminator (SRE)](#13-Short-Read-Eliminator-(SRE))
	- [1.4 Library Preparation](#14-Library-Preparation)
	- [1.5 ONT Sequencing](#15-ONT-Sequencing)
 - [**2. Library Analysis**](#part-2-library-analysis)
	- [2.1 NanoStat](#21-nanostat)
 - [**3. Assembly**](#part-3-assembly)
  	- [3.1 NextDenovo](#31-NextDenovo)
  	- [3.2 Read correction with Canu](#32-Read-correction-with-canu)
  	- [3.2 Assembly software](#33-Assembly-software)
  	- [3.3 Assembly polishing](#34-Assembly-polishing)
 - [**4. Evaluation**](#part-4-evaluation)
 	-  [4.1 QUAST](#41-QUAST)
 	-  [4.2 BUSCO](#42-BUSCO)
 	-  [4.3 Decontamination](#43-Decontamination)
    		- [4.3.1 SIDR](#SIDR)
    		- [4.3.2 Blobology](#432-Blobology)
 - [**5. Annotation**](#part-5-annotation)
   	- [5.1 Characterizing Repeats and Transposable Elements](#51-Characterizing-Repeats-and-Transposable-Elements)
    		- [5.1.1 RepeatModeler](#511-RepeatModeler)
    		- [5.1.2 EDTA](#512-EDTA)
  	- [5.2 Protein-coding gene annotation with BRAKER2](#52-Protein-coding-gene-annotation-with-BRAKER2)
    		- [5.2.1 Align RNASeq with STAR](#521-align-rnaseq-with-star)
    		- [5.2.2 Run BRAKER](#522-run-braker)
    	- [5.3 Protein-coding gene annotation with MAKER2](#53-Protein-coding-gene-annotation-with-MAKER2)
        	- [5.3.1 Create profile with BUSCO or CEGMA ](#531-create-profile-with-busco-or-cegma)
        	- [5.3.2 Create GeneMark hmm file](#532-create-genemark-hmm-file)
        	- [5.3.3 Run MAKER2](#533-run-maker2)
        	- [5.3.4 Run MAKER2 iteratively](#534-run-maker2-iteratively)
		- [5.3.5 Run EVidenceModeler with the braker2 and MAKER2 output](#535-run-evidencemodeler-with-the-braker2-and-maker2-output)
   	- [5.4 Annotation statistics with AGAT](#Annotation-statistics-with-AGAT)
   		- [5.4.1 Install AGAT via conda](#541-install-agat-via-conda)
   		- [5.4.2 Count genes and other features](#542-count-genes-and-other-features)
  - [**6. Upload data to NCBI](#6-Upload-data-to-NCBI)

## PART 1: Wet lab protocols

Assembling a genome sequence starts with the organism. Here, we are targeting nematodes that live in culture in the lab.

### 1.1 Nematode culture 
(from [Sutton et al. 2021](https://www.biorxiv.org/content/10.1101/2020.05.05.079327v2))

1) Nematodes grow on two 100mm NGM (nematode growth medium) plates seeded with E. coli OP50. 
2) Worms are harvested by washing plates with M9 minimal media into 15mL conical tubes. 
3) Samples are rocked on a tabletop rocker for 1 hour before being centrifuged to pellet worms. 
4) The supernatant is removed and tubes refilled with sterile M9, mixed and pelleted by centrifugation again. 
5) This process needs to be repeated five times or until the supernatant is clear after centrifugation. 
6) The pellet needs to be moved to 2mL tubes and frozen at ???20??C until extraction. 

### 1.2 High Molecular Weight DNA extraction

There are two methods we have used with good success. Phenol chloroform is a standard protocol and the Nanobind kit is available from [Circulomics](https://www.circulomics.com/).

### 1.2.1 Phenol Chloroform

1) Worm pellets need to thaw to room-temperature then flash freeze them in liquid nitrogen; repeat three times. 
2) Place worm pellets in 1.2mL of lysis buffer solution (100mM EDTA, 50mM Tris, and 1%SDS) and 20uL of Proteinase K (100mg/mL).
3) Place tubes on a 56C heat block for 30 minutes with shaking. 
4) Extract genomic DNA was then isolated using a modified phenol-chloroform extraction. 
5) FILL IN HERE

### 1.2.2 Circulomics Nanobind

Follow the kit protocols.

After you have HMW DNA:
1) Measure DNA concentrations and purity with a Qubit and NanoDrop?? 1000 spectrophotometer, respectively. The Qubit does high-accuracy quantification and the NanoDrop tells us about the quality of the DNA extract.
2) Visualize your DNA on a 0.8% agarose gel to verify high-molecular weight gDNA. 

### 1.3 Short Read Eliminator (SRE)

Short fragments will be preferentially sequenced but do not provide good information for assembly. We use the [Short Read Eliminator](https://www.circulomics.com/) for DNA size selection according to manufacturer guidelines.

### 1.4 Library Preparation

1) Prepare DNA libraries using the SQK-LSK109 ligation sequencing kit and load on to R9.4.1 RevD flow cells. 
2) Modify the recommended protocol from ONT by replacing the first AmpureXP bead clean step with an additional treatment with the Short Read Eliminator Kit. 

### 1.5 ONT Sequencing

1) Load approximately 700ng of gDNA from each library on to a flow cell and sequenced for 48 hours on a GridION X5 platform.
2) Basecalling is performed by Guppy v.4.0.11 set to high-accuracy mode. 

## PART 2: Library Analysis

After sequencing we need to move our data from the gridION to another computing system for assembly and analysis.
We will start to use UAHPC, the central computing server. Below is an example job submission script, we will discuss how the system works.

	#!/bin/bash

	#SBATCH -J #job name

	#SBATCH -p #partition (a group of nodes). We will use main, owners, or highmem

	#SBATCH --qos #This tells what sort of resource allowances you have for a particular partition. We will use main, or jlfierst

	#SBATCH -N #number of nodes (1 node = 1 machine/computer) We typically will only use 1)

	#SBATCH -c #cpus per task. A cpu is a processor and most nodes have numerous cpus. We typically won???t need to use this option because most of our programs don???t need a certain number of cpus allocated for specific tasks.

	#SBATCH -n #number of cpus/cores/threads/processors (these terms are often used interchangeably) required for the job to run. We often use between 1-16. Note that certain combinations of partition/quos will have a maximum number of cpus allowed for use.

	#SBATCH --mem #Amount of memory(RAM) requested per node. You must specify a value and unit (20G(20gigabytes)). Default unit is M(megabyte) We can specifically ask for a certain amount of memory per cpu, but usually don???t need to.

	#SBATCH -o %A.%a.out #STDOUT output

	#SBATCH -e %A.%a.err #STDERR output

	#SBATCH ???mail-type #emails your if certain things happen with your job. the type ???ALL??? is typically used and will inform you when your job begins, ends, fails, has invalid dependencies, times out, or gets requeued for some reason.

	#SBATCH ???mail-user #the username and/or email address to send notifications to. 

Logging on from your terminal:

	$ ssh -l username uahpc.ua.edu

Change to the /jlf drive, make a directory for your data and change into the directory

	$ cd /jlf
	
	$ mkdir {your mybama name}

	$ cd {your my bama name}

UAHPC is a large system with a single log on point (the head node). You will need to write job submission scripts for longer running jobs. For shorter jobs you can request an interactive node

	$ srun --pty -p [partition_name] /bin/bash
	
This will move you from the head node to a compute node but you will interact at the command prompt. The partition_name tells the system which resources you need.

UNIX refresher here: https://github.com/BamaComputationalBiology/IntroToBioinfo/blob/master/1.LifeInTheShell.md

### 2.1 [Nanostat]

https://github.com/wdecoster/nanostat

NanoStat measures some statistics about our ONT library. It is fairly quick but we need to get on an interactive node so the head node doesn't get 'bogged down' (way too slow).

Use vi to create your job script

	$ vi nanostat.sh

Type 'i' for insertion and enter the following:

	#!/bin/bash
	
	#SBATCH -J nanoStat #job name
	#SBATCH -p long
	#SBATCH --qos long
	#SBATCH -n 8
	#SBATCH -o %A.%a.out #STDOUT output
	#SBATCH -e %A.%a.err #STDERR output
	#SBATCH ???mail-user {your mybama email} 

	/jlf/jdmillwood/anaconda3/bin/NanoStat --fastq [ONT_reads.fastq]  --name [name_for_output] --outdir [Directory_for_output] --readtype 1D --threads 16

Hit the 'esc' key, type ':wq' (no quotation marks) to exit and save. To submit your job type

	$ sbatch < nanostat.sh

## PART 3: Assembly

ONT libraries have large numbers of incorrectly called nucleotides, insertions and deletions. Our group has found the best protocol is to correct ONT sequence reads, assemble and then polish. If you have high long read coverage (>30x >10kb) you can use NextDenovo/NextPolish; otherwise we use Canu/Flye/Pilon. It's worth experimenting with different protocols too and evaluating your assembled sequence.

### 3.1 NextDenovo (https://github.com/Nextomics/NextDenovo)

You will need drmaa

	$ pip install drmaa
	$ wget https://github.com/natefoo/slurm-drmaa/releases/download/1.1.0/slurm-drmaa-1.1.0.tar.gz
	$ tar -vxzf slurm-drmaa-1.1.0.tar.gz
	$ cd slurm-drmaa-1.1.0
	$ ./configure && make && make install
	$ export DRMAA_LIBRARY_PATH=/home/{mybama name}/slurm-drmaa-1.1.0/slurm_drmaa/.libs/libdrmaa.so.1

Tell NextDenovo where the ONT libraries are

	$ ls {ONT-libraries} > input.fofn

The output should be a file with a list of libraries, one per line. To view it:

	$ more input.fofn
	
Create a run.cfg (configuration) file. Use vi to create a new file

	$ vi run.cfg
	
Once you are in the file hit 'i' for insertion mode, copy and paste the configurations below:

	[General]
	job_type = local # local, slurm, sge, pbs, lsf
	job_prefix = nextDenovo
	task = all # all, correct, assemble
	rewrite = yes # yes/no
	deltmp = yes 
	parallel_jobs = 8 # number of tasks used to run in parallel
	input_type = raw # raw, corrected
	read_type = ont # clr, ont, hifi
	input_fofn = input.fofn # your file created above
	workdir = {species name} # output

	[correct_option]
	read_cutoff = 1k # Discard anything below 1000bp
	genome_size = 80M # estimated genome size, we're not sure for Oscheius but we'll use this as a first estimate
	sort_options = -m 20g -t 15
	minimap2_options_raw = -t 8
	pa_correction = 8 # number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4 bytes of memory usage.
	correction_options = -p 15

	[assemble_option]
	minimap2_options_cns = -t 8 
	nextgraph_options = -a 1
	
If you need to edit you can tab around and change your file. To save hit 'esc' (the escape key) for command mode then ':wq' to save and exit.
	
Create a job submission script to assemble with NextDenovo

	$ vi assemble.sh
	
Remember to hit 'i' for insertion mode and type the following:

	#!/bin/bash
	
	#SBATCH -J nextDenovo #job name
	#SBATCH -p highmem
	#SBATCH --qos jlfierst
	#SBATCH -n 8
	#SBATCH --mem 150G
	#SBATCH -o %A.%a.out #STDOUT output
	#SBATCH -e %A.%a.err #STDERR output
	#SBATCH ???mail-user {your mybama email} 
	
	module load bio/nextdenovo
	
 	nextDenovo run.cfg

To exit hit the 'esc' key then type ':wq'. To submit to the UAHPC queue system type

	$ sbatch < assemble.sh

To check your job type

	$ squeue -u {mybama name}

The output will be in the {species name}/03.ctg_graph directory. There are two files, nd.asm.fasta is the assembled genome sequence and nd.asm.fasta.stat gives you some statistics regarding the contiguity and quality of the assembled sequence.

Polishing with NextPolish follows a similar workflow. First, tell NextPolish where the ONT libraries are

	$ ls {ONT libraries} > lgs.fofn

Then tell NextPolish where the Illumina short reads are

	$ ls {Illumina libraries} {Illumina libraries} > sgs.fofn

Edit the run.cfg file for your sequence

	[General]
	job_type = local
	job_prefix = nextPolish
	task = default
	rewrite = yes
	rerun = 3
	parallel_jobs = 2
	multithread_jobs = 3
	genome = {species name}/03.ctg_graph/nd.asm.fasta
	genome_size = auto
	workdir = {species name}
	polish_options = -p {multithread_jobs}

	[sgs_option]
	sgs_fofn = ./sgs.fofn
	sgs_options = -max_depth 100

	[lgs_option]
	lgs_fofn = ./lgs.fofn
	lgs_options = -min_read_len 5k -max_depth 100
	lgs_minimap2_options = -x map-ont

Your assembled, polished sequence will be at {species name}/genome.nextpolish.fasta

See https://nextdenovo.readthedocs.io/en/latest/OPTION.html for a detailed introduction about all the parameters

If you have completed this portion with NextDenovo/NextPolish you can skip down to (4) and begin your evaluations.

### 3.2 Read correction with Canu, Assemble with Flye

Read correction with Canu using the canu-correct module. Canu assembly can be very slow and our group has found Flye to be much faster with similar or improved accuracy

	$ vi canu_assemble.sh

Hit 'i' for insertion and type the following:

	#!/bin/bash
	
	#SBATCH -J canu_flye #job name
	#SBATCH -p highmem
	#SBATCH --qos jlfierst
	#SBATCH -n 8
	#SBATCH --mem=128G
	#SBATCH -o %A.%a.out #STDOUT output
	#SBATCH -e %A.%a.err #STDERR output
	#SBATCH ???mail-user {your mybama email}
	
	module load bio/canu/2.1
	module load java/1.8.0
	
	canu -correct -p [outfile_name_prefix] -d [out_directory] genomeSize=80M useGrid=false -nanopore [reads].fastq
        
	module load bio/bioinfo-gcc
	module load python/python3/3.6.5
	/jlf/jdmillwood/Flye/bin/flye --nano-corr /jlf/[mybamausername]/[out_directory]/[prefix].correctedReads.fasta.gz -o flye_try -t 8 --genome-size 80M
	
Exit by hitting the 'esc' button, typing ':wq'.

Submit your job

	$ sbatch < canu_assemble.sh

### 3.4 Assembly polishing

ONT assemblies have small errors that can be addressed through iterative polishing with Illumina libraries.
Copy illumina data from /jlf/CRE_UA/raw_reads/illumina/ if you haven't already

	$ vi pilon.sh
Hit 'i' for insertion and type the following (edit FORWARD, REVERSE, LINE NAME, and GENOME):


	#!/bin/bash
	#SBATCH -n 8
	#SBATCH -p highmem
	#SBATCH --qos jlfierst
	#SBATCH --mem=128G
	#SBATCH --job-name=pilon
	#SBATCH -e %J.err
	#SBATCH -o %J.out


	#Load Modules
	module load bio/bioinfo-gcc
	module load bio/samtools/1.9


	FORWARD=[PATH TO FASTQ_1]
	REVERSE=[PATH TO FASTQ_2]
	LINE_NAME=PB127 ## YOUR LINE
	mkdir ./pilon_out/

	## ROUND 1 ##
	GENOME=[ path to assembled genome]
	#index genome 
	bwa index ${GENOME}
	#align reads
	bwa mem -t 8 -M ${GENOME} ${FORWARD} ${REVERSE}  > ./pilon_out/bwa.sam
	#sam to bam
	samtools view -Sb ./pilon_out/bwa.sam > ./pilon_out/bwa.bam
	##Sort and index the BAM 
	samtools sort ./pilon_out/bwa.bam -o ./pilon_out/bwa.sort
	samtools index ./pilon_out/bwa.sort

	module load java/1.8.0
	module load bio/bioinfo-java
	##Pilon it 
	java -Xmx12G -jar /share/apps/bioinfoJava/pilon-1.22.jar --genome ${GENOME} --frags ./pilon_out/bwa.sort --output ./pilon_out/${LINE_NAME}_pilon1
	module unload java/1.8.0
	module unload bio/bioinfo-java
	
	## ROUND 2 ##
	GENOME=./pilon_out/${LINE_NAME}_pilon1.fasta 
	#index genome 
	bwa index ${GENOME}
	#align reads
	bwa mem -t 8 -M ${GENOME} ${FORWARD} ${REVERSE}  > ./pilon_out/bwa.sam
	#sam to bam
	samtools view -Sb ./pilon_out/bwa.sam > ./pilon_out/bwa.bam
	##Sort and index the BAM 
	samtools sort ./pilon_out/bwa.bam -o ./pilon_out/bwa.sort
	samtools index ./pilon_out/bwa.sort

	module load java/1.8.0
	module load bio/bioinfo-java
	##Pilon it 
	java -Xmx12G -jar /share/apps/bioinfoJava/pilon-1.22.jar --genome ${GENOME} --frags ./pilon_out/bwa.sort --output ./pilon_out/${LINE_NAME}_pilon2
	module unload java/1.8.0
	module unload bio/bioinfo-java

	## ROUND 3 ##
	GENOME=./pilon_out/${LINE_NAME}_pilon2.fasta 
	#index genome 
	bwa index ${GENOME}
	#align reads
	bwa mem -t 8 -M ${GENOME} ${FORWARD} ${REVERSE}  > ./pilon_out/bwa.sam
	#sam to bam
	samtools view -Sb ./pilon_out/bwa.sam > ./pilon_out/bwa.bam
	##Sort and index the BAM 
	samtools sort ./pilon_out/bwa.bam -o ./pilon_out/bwa.sort
	samtools index ./pilon_out/bwa.sort

	module load java/1.8.0
	module load bio/bioinfo-java
	##Pilon it 
	java -Xmx12G -jar /share/apps/bioinfoJava/pilon-1.22.jar --genome ${GENOME} --frags ./pilon_out/bwa.sort --output ./pilon_out/${LINE_NAME}_pilon3
	module unload java/1.8.0
	module unload bio/bioinfo-java


	## ROUND 4 ##
	GENOME=./pilon_out/${LINE_NAME}_pilon3.fasta 
	#index genome 
	bwa index ${GENOME}
	#align reads
	bwa mem -t 8 -M ${GENOME} ${FORWARD} ${REVERSE}  > ./pilon_out/bwa.sam
	#sam to bam
	samtools view -Sb ./pilon_out/bwa.sam > ./pilon_out/bwa.bam
	##Sort and index the BAM 
	samtools sort ./pilon_out/bwa.bam -o ./pilon_out/bwa.sort
	samtools index ./pilon_out/bwa.sort

	module load java/1.8.0
	module load bio/bioinfo-java
	##Pilon it 
	java -Xmx12G -jar /share/apps/bioinfoJava/pilon-1.22.jar --genome ${GENOME} --frags ./pilon_out/bwa.sort --output ./:pilon_out/${LINE_NAME}_pilon4
	
	
Exit by hitting the 'esc' button, typing ':wq'.

Submit your job

	$ sbatch < pilon.sh


## PART 4: Evaluation

### 4.1 QUAST

[QUAST](http://quast.sourceforge.net/quast.html): Quality Assessment Tool for Genome Assemblies. QUAST can be used to evaluate genome statistics like N50 and misassemblies relative to a reference. If you don't have a reference it will estimate statistics like N50, N90, L50 and L90.

$ vi quast.sh
Hit 'i' for insertion and type the following (edit the locations with your line name):
QUAST without a reference

	#!/bin/bash
	#SBATCH -n 4
	#SBATCH -p highmem
	#SBATCH --qos jlfierst
	#SBATCH --job-name=quast
	#SBATCH --mem=100G

	module load bio/quast/5.0

	#quast.py -t 4 --eukaryote --plots-format pdf ./pilon_out/PB127_pilon4.fasta -o ./PB127_quast/
	quast.py -t 4 --eukaryote --plots-format pdf ./pilon_out/DF5018_pilon4.fasta -o ./DF5018_quast/


	
Exit by hitting the 'esc' button, typing ':wq'.

Submit your job

	$ sbatch < quast.sh
	
	
QUAST with a reference

	$ python {quast location}/quast.py -t 12 --plots-format pdf -r {reference genome} {assembled sequence} -o {output name}
	
### 4.2 BUSCO

[BUSCO](http://busco.ezlab.org/) searches assembled genome sequences for a set of genes thought to be conserved in single copy in a group of organisms.

First step is to copy the augustus config directory to your directory 

Canu-flye:
 	
	cp -r /share/apps/augustus/augustus-3.3.2/config/ /jlf/USERNAME/augustus_config/
Nextdenovo:

	cp -r /share/apps/augustus/augustus-3.3.2/config/ /home/USERNAME/augustus_config/

	
$ vi busco.sh
Hit 'i' for insertion and type the following (edit the locations with your line name):

		

	#!/bin/bash
	#SBATCH -n 4
	#SBATCH -p highmem
	#SBATCH --qos jlfierst
	#SBATCH --job-name=busco
	#SBATCH --mem=100G


	module load bio/busco
	export AUGUSTUS_CONFIG_PATH="/jlf/USERNAME/augustus_config/" # replace username

	busco -c 4 -m genome -i [PATH_To_POLISHED_GENOME] -o busco_[LINE] --lineage_dataset nematoda_odb10 --augustus_species caenorhabditis


### 4.3 Decontamination

Blast
vi blast.sh

Hit 'i' for insertion and type the following:

	#!/bin/bash
	#SBATCH -J BLAST
	#SBATCH --qos jlfierst
	#SBATCH -p highmem
	#SBATCH -n 4
	#SBATCH -o %J.out
	#SBATCH -e %J.err


	module load compilers/gcc/5.4.0
	module load bio/blast/2.9.0
	BLASTDB=/jlf/newblastdb
	echo $BLASTDB 

	blastn  -task megablast -query [PATH_TO_POLISHED_GENOME] -db nt -outfmt '6 qseqid staxids' -culling_limit 5 -evalue 1e-25 -out [LINE].blast.out



### 4.3.1 [SIDR]()

### 4.3.2 [Blobology](http://drl.github.io/blobtools)

Generating Blob plots 
=======
Blobtoolkit
 
Create Blob Directory

	$ /path/to/blobtoolkit/blobtools2/blobtools create --fasta /path/to/genome.fasta --taxdump /path/to/blobtoolkit/taxdump {assembly}

Run blastn

	$ /path/to/blastn -db nt -query genome.fasta -outfmt '6 qseqid staxids bitscore std' \
	-max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 16 -out blast.out

Run blastx

	$ /path/to/diamond blastx --query genome.fasta --db reference_proteomes.fasta.gz \
	--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
	--sensitive --max-target-seqs 1 --evalue 1e-25 --threads 16 > diamond.out

Map Illumina reads 

	$ minimap2 -ax map-ont -t 16 genome.fasta /path/to/raw/reads/reads.fastq | samtools sort -@16 -O BAM -o assembly.reads.bam 

Add all data files to Blob Directory

	$ path/to/blobtools2/blobtools add --hits blast.out --hits diamond.out \
	--taxrule bestsumorder --taxdump ~/taxdump --cov assembly.reads.bam --busco /path/to/busco/run/full_table.tsv {assembly}

Open BlobToolKit Viewer to view BlobPlot

	$ path/to/blobtools2/blobtools host /path/to/working/blob/info/directory

Filtering Genomes with Whitelist

	## awk to split .fasta 

	#!/bin/bash

	awk '/^>scaffold/ {OUT=substr($0,2) ".fasta"}; OUT {print >OUT}' final.assembly.fasta

	cat scaffolds_to_keep.txt | while read line; do cat "scaffold_"$line".fasta" ; done > kept_scaffolds.fasta

	rm "scaffold"*".fasta"

Checking the output

	$ cat kept_scaffolds.fasta | grep '^>' | tr -d '\>scaffold\_' | sort -n | diff - scaffolds_to_keep.txt 

## PART 5: Annotation
    
### 5.1 Characterizing Repeats and Transposable Elements

### 5.1.1 [RepeatModeler](http://www.repeatmasker.org/)

RepeatModeler creates a custom library of repeats found in your assembled genome sequence. We are using RepeatModeler through [TE-Tools](https://github.com/Dfam-consortium/TETools)

	Launch the container

	$ ./dfam-tetools.sh

	Build the database

	$ BuildDatabase -name [species_name] [genome.fasta]
	
	Run RepeatModeler for de novo repeat identification and characterization

	$ RepeatModeler -pa 8 -database [species_name]

	Use the queryRepeatDatabase.pl script inside RepeatMasker/util to extract Rhabditida repeats

	$ queryRepeatDatabase.pl -species rhabditida | grep -v "Species:" > Rhabditida.repeatmasker

	Combine the files to create a library of de novo and known repeats

	$ cat RM*/consensi.fa.classified Rhabditida.repeatmasker > [species_name].repeats

Finally, RepeatMasker creates a masked version of our assembled genome sequence. There are two versions of masking. Hard-masking means we replace every nucleotide in a repeat region with 'N' and soft-masking means we replace the normally capitalized nucleotides with lower-case nucleotides in repeat regions. Here we soft-mask (-xsmall) and do not mask low complexity elements (-nolow).

	$ RepeatMasker -lib [species_name].repeats -pa 8 -xsmall -nolow [genome.fasta] 

### 5.1.2 [EDTA](https://github.com/oushujun/EDTA)



### 5.2 Protein-coding gene annotation with BRAKER2

#### 5.2.1 Align RNASeq with [STAR](https://github.com/alexdobin/STAR)

```
# Generate genome index
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir [species_dir] --genomeSAindexNbases 12 --genomeFastaFiles [species_genome]

# Map the reads
STAR --runThreadN 12 --genomeDir [species_dir] --outSAMtype BAM Unsorted --twopassMode Basic --readFilesCommand zcat \ # if reads are zip-compressed
--readFilesIn [File_1.fastq.gz] [File_2.fastq.gz] 
```

#### 5.2.2 Run [BRAKER](http://exon.gatech.edu/genemark/braker1.html)

```
braker.pl \
--workingdir=[output_dir] \
--species=[species_name] \
--cores=8 \
--genome=[assembly.masked] \
--bam=Aligned.out.bam \
--prot_seq=[protein.fasta] \
--prg=gth \
--GENEMARK_PATH=[GENEMARK_dir] \
--softmasking

gtf2gff.pl < braker.gtf --gff3 --out=braker.gff3
```

After braker2 is finished we will create protein and fasta files from the native .gtf output. We can align the protein and/or mRNA sequences to the NCBI BLAST databases to check the validity of our annotations.

```
#!/bin/bash


DIR="[BRAKER2 .gtf location]"

conda activate agatenv

agat_sp_extract_sequences.pl -f ${DIR}_filtered.fasta --mrna -g braker.gff3 -o braker.mRNA.fasta
agat_sp_extract_sequences.pl -f ${DIR}_filtered.fasta -p -g braker.gff3 -o braker.protein.fasta

tblastn -db nt -query braker.protein.fasta \
-outfmt '6 qseqid qlen staxids bitscore std sscinames sskingdoms stitle' \
-num_threads 4 -evalue 0.01 -max_target_seqs 2 -out blast.out
```

### 5.3 Protein-coding gene annotation with MAKER2

We can use the braker2-trained Augustus and GeneMark in our MAKER2 annotations or follow the protocol below for training.

### 5.3.1 

Use BUSCO to create Augustus profile **(new)** or CEGMA to create SNAP profile.

#### 5.3.1.1 Train [Augustus](http://bioinf.uni-greifswald.de/augustus/) with [BUSCO](http://busco.ezlab.org/)

Run BUSCO
```
python run_BUSCO.py -i [assembly.masked] -o [busco_euk] -l [eukaryota_odb9 dir] -m genome -c 16 -sp caenorhabditis --long
```
Rename the parameter files
```
# Files will be in run_busco_euk/augustus_output/retraining_parameters
# Rename the folder to [species_name] 
mv retraining_parameters [species_name]
# Rename the files to [species_name] from [Busco_name] (e.g just before _metapars)
cd [species_name]
for a in * ; do rename 's/[Busco_name]/[species_name]/' $a ; done
for a in * ; do sed -i 's/[Busco_name]/[species_name]/' $a ; done
# Move the folder to augustus config/species
```

#### 5.3.1.2 Train [SNAP](http://korflab.ucdavis.edu/software.html) with [CEGMA](http://korflab.ucdavis.edu/Datasets/cegma/)

Run CEGMA
```
cegma -g [assembly.masked] -T 8
```

Create the SNAP hmm file
```
cegma2zff [output.cegma.gff] [assembly.masked] # (from MAKER software)
fathom genome.ann genome.dna -categorize 1000
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl [Cv] . > [Cv.cegmasnap.hmm]
```

#### 5.3.2 Create [GeneMark](http://exon.gatech.edu/GeneMark/) hmm file

Run GeneMark, the output will be named ```gmhmm.mod```.
```
gmes_petap.pl --ES --max_intron 20000 --cores 8 --sequence [assembly.masked]
```

#### 5.3.3 Run [MAKER2](http://www.yandell-lab.org/software/maker.html)

Create the config files
```
maker -CTL
```

Edit ```maker_opts.ctl``` (an example is below)
```
#-----Genome (these are always required)
genome=/data/Caenorhabditis_genomes/C_remanei/356/RepeatMasker/356_filtered_SIDR.fasta.masked #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff=# MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=356_Trinity.fasta#set of ESTs or assembled mRNA-seq in fasta format
altest= /data/PX506/GCA_010183535.1_CRPX506_cds_from_genomic.fna, Caenorhabditis_briggsae.CB4.20.cdna.all.fa, Caenorhabditis_elegans.WBcel235.20.cdna.all.fa #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=uniprot_sprot.fasta#protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=#aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=#select a model organism for RepBase masking in RepeatMasker
rmlib= /data/Caenorhabditis_genomes/C_remanei/356/RepeatMasker/356.repeats #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=#pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=/data/Caenorhabditis_genomes/C_remanei/356/Maker/snap/rnd1/rnd1.zff.length50_aed0.5.hmm #SNAP HMM file
gmhmm= /data/jlfierst/356/braker/GeneMark-ET/gmhmm.mod #GeneMark HMM file
augustus_species=356 #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=1 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=20 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=5000 #skip genome contigs below this length (under 10kb are often useless)
```

Run MAKER2 and Extract the gff3 fasta files
```
#!/bin/bash

export PATH="/data/maker/bin:$PATH"
export AUGUSTUS_CONFIG_PATH="data/jlfierst/anaconda3/config/species:$PATH"

MAKERDIR="356"

maker -base ${MAKERDIR} maker_opts.ctl maker_bopts.ctl maker_exe.ctl -fix_nucleotides
gff3_merge  -d ${MAKERDIR}.maker.output/${MAKERDIR}_master_datastore_index.log
fasta_merge -d ${MAKERDIR}.maker.output/${MAKERDIR}_master_datastore_index.log
```

#### 5.3.4 Run MAKER2 iteratively



#### 5.3.5 Run EVidenceModeler with the braker2 and MAKER2 output

### 5.4 Annotation statistics with AGAT

AGAT(https://github.com/NBISweden/AGAT#installation) is a tool for annotation editing and evaluation. We will install via conda and use it to evaluate annotation statistics. AGAT creates conflicts with some other aspects of conda and we will install/activate it into its own environment to manage the conflicts.

#### 5.4.1 Install AGAT via conda(https://anaconda.org/bioconda/agat)

	conda create -n agatenv
	conda activate agatenv
	conda install -c bioconda agat

#### 5.4.2 Count genes and other features

	conda activate agatenv
	agat_sp_statistics.pl --gff {file}.gff3

### 5.5 Create pseudo-chromosomes with Ragtag

Caenorhabditis have highly conserved chromosome structure and we expect that large-scale rearrangements have not occurred between closely related strains. This is an assumption and it may not be true but we will operate under this for now. If we assume this we can use RagTag (https://github.com/malonge/RagTag) to scaffold our fragmented assembly based on a chromosome-scale assembly for a close relative, C. remanei PX506 (https://www.ncbi.nlm.nih.gov/genome/253?genome_assembly_id=771236).

RagTag is very fast, below is a script that can align our fragmented assembly to the chromosome-scale PX506 and create a new set of coordinates for our annotated genes.

	#!/bin/bash
	
	GENOME=534

	ragtag.py scaffold GCA_010183535.1_CRPX506_genomic.fna {GENOME}.fasta -t 8

	ragtag.py updategff braker.gff3 ./ragtag_output/ragtag.scaffold.agp > ragtag.{GENOME}.gff3


## PART 6: Upload data to NCBI


