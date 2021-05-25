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
    - [5.3 Protein-coding gene annotation with MAKER2](#53-Protein-coding-gene-annotation-with-MAKER2)
        - [5.3.1 Create profile with BUSCO or CEGMA ](#531-create-profile-with-busco-or-cegma)
        - [5.3.2 Create GeneMark hmm file](#532-create-genemark-hmm-file)
        - [5.3.3 Run MAKER2](#533-run-maker2)
        - [5.3.4 Run MAKER2 iteratively](#534-run-maker2-iteratively)
	- [5.3.5 Run EVidenceModeler with the braker2 and MAKER2 output](#535-run-evidencemodeler-with-the-braker2-and-maker2-output)
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
 
	$ flye --nano-corr [sample] --out-dir out_nano --threads 8 --genome-size [genome estimate] 

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

Finally, RepeatMasker creates a masked verion of our assembled genome sequence. There are two versions of masking. Hard-masking means we replace every nucleotide in a repeat region with 'N' and soft-masking means we replace the normally capitalized nucleotides with lower-case nucleotides in repeat regions. Here we soft-mask (-xsmall) and do not mask low complextiy elements (-nolow).

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
```
Run braker2 again to add the UTRs and generate a .gff3 file

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
--softmasking \
--addUTR=on \
--useExisting \
--gff3
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
altest= /data/PX506/GCA_010183535.1_CRPX506_cds_from_genomic.fna, /data/jlfierst/aciss_2016/MAKER/evidence_files/Caenorhabditis_briggsae.CB4.20.cdna.all.fa, /data/jlfierst/a
ciss_2016/MAKER/evidence_files/Caenorhabditis_elegans.WBcel235.20.cdna.all.fa #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein= /data/jlfierst/aciss_2016/MAKER/evidence_files/uniprot_sprot.fasta#protein sequence file in fasta format (i.e. from mutiple organisms)
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

## PART 6: Upload data to NCBI


