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
	- [2.1 PoreChop](#21-PoreChop)
	- [2.2 poretools](#22-poretools)
 - [**3. Assembly**](#part-3-assembly)
 - 	- [3.1 Read correction](#31-Read-correction)
 - 	- [3.2 Assembly software](#32-Assembly-software)
 - 	- [3.3 Assembly polishing](#33-Assembly-polishing)
 - [**4. Evaluate Assembly**](#part-4-evaluate-assembly)
 -  [4.1 QUAST](#41-QUAST)
 -  [4.2 Decontamination](#42-Decontamination)
    - [4.2.1 SIDR](#SIDR)
    - [4.2.2 Blobology](#23-Blobology)
 - [**5. Annotation**](#part-5-annotation)
    - [5.1 Characterizing Repeats and Transposable Elements](#51-Characterizing-Repeats-and-Transposable-Elements)
    	- [5.1.1 RepeatModeler](#511-RepeatModeler)
    	- [5.1.2 EDTA](#512-EDTA)
    - [5.2 Protein-coding gene annotation with BRAKER2](#52-Protein-coding-gene-annotation-with-BRAKER2)
        - [5.2.1 Align RNASeq with STAR](#521-align-rnaseq-with-star)
        - [5.2.2 Run BRAKER](#522-run-braker)
        - [5.2.3 Add UTR info to the GFF3 file](#523-add-utr-info-to-the-gff3-file)
    - [5.3 Protein-coding gene annotation with MAKER2](#53-Protein-coding-gene-annotation-with-MAKER2)
        - [5.3.1 Create profile with BUSCO or CEGMA ](#531-create-profile-with-busco-or-cegma)
        - [5.3.2 Create GeneMark hmm file](#532-create-genemark-hmm-file)
        - [5.3.3 Run MAKER2](#533-run-maker2)
        - [5.3.4 Run Augustus](#534-run-augustus)
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
6) The pellet needs to be moved to 2mL tubes and frozen at –20°C until extraction. 

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
1) Measure DNA concentrations and purity with a Qubit and NanoDrop® 1000 spectrophotometer, respectively. The Qubit does high-accuracy quantification and the NanoDrop tells us about the quality of the DNA extract.
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

### 2.1 [Porechop](https://github.com/rrwick/Porechop)

We use porechop to discard reads with sequencing adapters.

Basic adapter trimming:
     
     	$ porechop -i input_reads.fastq.gz -discard-middle -o output_reads.fastq.gz

#### 2.2 [poretools](https://poretools.readthedocs.io/en/latest/index.html)

Poretools gives us a suite of utilities for assessing the quality, size and distribution of our library.

## PART 3: Assembly

### 3.1 Read correction

ONT libraries have large numbers of incorrectly called nucleotides, insertions and deletions. Our group has found the best protocol is to correct ONT sequence reads with Canu using the canu-correct module.

	$ canu -correct -p [sample] -d [out_directory_name] -nanopore reads.fastq genomeSize=XXXM

### 3.2 Assembly software

Canu assembly can be very slow and our group has found Flye to be much faster with similar or improved accuracy

	$ flye

### 3.3 Assembly polishing

ONT assemblies have small errors that can be addressed through iterative polishing with Illumina libraries.

	#!/bin/bash

	GENOME=/path/to/genome.fasta

	FORWARD=/path/to/reads1.fq

	REVERSE=/path/to/reads2.fq

	#index genome bwa index ${GENOME}

	bwa mem -M -t 48

	${GENOME} ${FORWARD} ${REVERSE} > bwa.sam

	samtools view -Sb bwa.sam > bwa.bam

	#Sort and index the BAM 
		
	samtools sort bwa.bam -o bwa.sort samtools index bwa.sort

	#Pilon it 
	
	java -Xmx300G -jar /data/jmsutton/anaconda3/share/pilon-1.23-0/pilon-1.23.jar --genome ${GENOME} --frags bwa.sort --output r1 --outdir r1

## PART 4: Evaluate Assembly

### 4.1 QUAST and BUSCO

[QUAST](http://quast.sourceforge.net/quast.html): Quality Assessment Tool for Genome Assemblies. QUAST can be used to evaluate genome statistics like N50 and misassemblies relative to a reference.

[BUSCO](http://busco.ezlab.org/) searches assembled genome sequences for a set of genes thought to be conserved in single copy in a group of organisms.

### 4.2 Decontamination

### 4.2.1 [SIDR]()

### 4.2.2 [Blobology](http://drl.github.io/blobtools)

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

	Use the queryRepeatDatabase.pl script inside RepeatMasker/util to extract Rhabditidia repeats

	$ queryRepeatDatabase.pl -species rhabditida | grep -v "Species:" > Rhabditida.repeatmasker

	Combine the files to create a library of de novo and known repeats

	$ cat RM*/consensi.fa.classified Rhabditida.repeatmasker > [species_name].repeats

Finally, RepeatMasker creates a masked verion of our assembled genome sequence. There are two versions of masking. Hard-masking means we replace every nucleotide in a repeat region with 'N' and soft-masking means we replace the normally capitalized nucleotides with lower-case nucleotides in repeat regions.

	$ RepeatMasker -lib [species_name].repeats -pa 8 [genome.fasta] 

### 5.1.2 [EDTA](https://github.com/oushujun/EDTA)



### 5.2 Protein-coding gene annotation with BRAKER2

#### 5.2.1 Align RNASeq with [STAR](https://github.com/alexdobin/STAR)

```
# Generate genome index
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir [genome_dir] \
--genomeFastaFiles [assembly.masked]

# Map the reads
STAR --runThreadN 8 --genomeDir [genome_dir] --outSAMtype BAM Unsorted \
--twopassMode Basic --readFilesCommand zcat \
--readFilesIn [RNASeq_1.sanfastq.gz] [RNASeq_2.sanfastq.gz] 
```

#### 5.2.2 Run [BRAKER](http://exon.gatech.edu/genemark/braker1.html)

```
braker.pl \
--workingdir=[output_dir] \
--species=[species_name] \
--cores=8 \
--genome=[assembly.masked] \
--bam=Aligned.out.bam
```

### 5.3 Protein-coding gene annotation with MAKER2

### 5.3.1 Create profile with BUSCO or CEGMA 

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

Edit ```maker_opts.ctl``` and the change the following
```
genome=[assembly.masked]
protein=[swissprot.fasta]
model_org=                      # Change value to blank
repeat_protein=                 # Change value to blank
snaphmm=[Cv.cegmasnap.hmm]      # If step 5.3.1.2
gmhmm=[genemark.mod]
augustus_species=[species_name] # If step 5.3.1.1
est2genome=1
protein2genome=1
pred_stats=1
keep_preds=1
single_exon=1
```

Run MAKER2 and Extract the gff3 fasta files
```
mpiexec -n 16 maker   # No. of threads

gff3_merge -d [datastore_index.log]
```

#### 5.3.4 Run [Augustus](http://bioinf.uni-greifswald.de/augustus/)

Create the training gff for Augustus (```maker2zff``` from MAKER, ```zff2gff3.pl``` from SNAP)
```
maker2zff -c 0 -e 0 [maker.gff3] 

zff2gff3.pl genome.ann | perl -plne 's/\t(\S+)$/\t\.\t$1/' >augustus.train.gff3
```

Run Augustus
```
autoAug.pl \
--species=[species_name] \
--genome=[assembly.masked] \
--trainingset=augustus.train.gff3 \
--useexisting -v --singleCPU \
--optrounds=3 &>augustus.out.txt
```
## PART 6: Upload data to NCBI


