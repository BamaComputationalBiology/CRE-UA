## Preface

This is the Fierst Lab protocol for ONT assembly and annotation of nematode genomes.

## Contents

 - [**1. Wet lab protocols**](#part-1-wetlab)
	- [1.1 Nematode culture](#11-Nematode-culture)
	- [1.2 High Molecular Weight DNA extraction](#12-HMW-DNA-extraction)
	-   [1.2.1 Phenol Chloroform](#121-Phenol-Chloroform)
	-   [1.2.2 Circulomics Nanobind](#122-Circulomics-Nanobind)
	- [1.3 Short Read Eliminator (SRE)](#13-Short-Read-Eliminator-(SRE))
	- [1.4 Library Preparation](#14-Library-Preparation)
	- [1.5 ONT Sequencing](#15-ONT-Sequencing)
 - [**2. Library Analysis**](#part-2-library-analysis)
	- [2.1 PoreChop](#21-PoreChop)
	- [2.2 poretools](#22-poretools)
 - [**3. Assembly**](#part-3-assembly)
 - [**4. Evaluate Assembly**](#part-4-evaluate-assembly)
 -  [4.1 QUAST](#41-QUAST)
 -  [4.2 Decontamination](#42-Decontamination)
    - [4.2.1 SIDR](#SIDR)
	  - [4.2.2 Blobology](#23-blobology)
		  - [4.2.2.1 Search NCBI nt with blast](#231-search-ncbi-nt-with-blast)
		  - [4.2.2.3 Make the TAGC plot](#233-make-the-tagc-plot)
 - [**5. Annotation**](#part-5-annotation)
    - [5.1 Mask the repeats with RepeatModeler and RepeatMasker](#51-mask-the-repeats-with-repeatModeler-and-repeatMasker)
    - [5.2 Annotation Pipeline if we have RNASeq data](#52-annotation-pipeline-if-we-have-rnaseq-data)
        - [5.2.1 Align RNASeq with STAR](#521-align-rnaseq-with-star)
        - [5.2.2 Run BRAKER](#522-run-braker)
        - [5.2.3 Add UTR info to the GFF3 file](#523-add-utr-info-to-the-gff3-file)
    - [5.3 Annotation Pipeline if we do not have RNASeq data](#53-annotation-pipeline-if-we-do-not-have-rnaseq-data)
        - [5.3.1 Create profile with BUSCO or CEGMA ](#531-create-profile-with-busco-or-cegma)
        - [5.3.2 Create GeneMark hmm file](#532-create-genemark-hmm-file)
        - [5.3.3 Run MAKER2](#533-run-maker2)
        - [5.3.4 Run Augustus](#534-run-augustus)
  - [**6. Upload data to NCBI](#6-Upload-data-to-NCBI)

## PART 1: Wet lab protocols

### 1.1 

### 1.2 

### 1.3 

### 1.4 

### 1.5 

## PART 2: Library Analysis

### 2.1 [Porechop](https://github.com/rrwick/Porechop)

Basic adapter trimming:
     
     porechop -i input_reads.fastq.gz -discard-middle -o output_reads.fastq.gz

#### 2.2 [poretools](https://poretools.readthedocs.io/en/latest/index.html)

Library qc

## PART 3: Assembly

### 3.1 

#### 3.1.1 

#### 3.1.2

### 3.2 Assembly software


#### 3.2.1 

#### 3.2.2 

## PART 4: Evaluate Assembly

### 4.1 QUAST and BUSCO

[BUSCO](http://busco.ezlab.org/)



### 4.2 [Decontamination]()

### 4.2.2 [Blobology](http://drl.github.io/blobtools)

#### 4.2.2.1 Search NCBI nt with blast

    blastn \
    -task megablast \
    -query [assembly_se.fasta] \
    -db nt \
    -culling_limit 2 \
    -out [assembly_se_nt.blastn] \
    -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
    -num_threads 16 \
    -evalue 1e-25 

#### 4.2.2.2 Make the TAGC plot

Create the BlobDB JSON file

    blobtools create \
    -i [assembly_se.fasta] \
    --cas [lib1.cas] --cas [lib2.cas] \ # if clc, or -y velvet if velvet
    -t [assembly_se_uniref.daa.tagc] \
    -t [assembly_se_nt.blastn] \
    --nodes [path_to/nodes.dmb] \       # nodes.dmp (required once)
    --names [path_to/names.dmb]         # names.dmp (required once)

Generate a table

    blobtools view -i [BlobDB] --hits --rank all > [BlobDB.table]

Create the plot

    blobtools blobplot -i [BlobDB]

We can now inspect the image
![Blob](http://i.imgur.com/LilPNJI.png)

## PART 5: Annotation

Rename the assembly scaffolds to something sensible with [assembly.rename.sh](https://github.com/GDKO/CGP-scripts) (requires fastaqual_select.pl in path)

    assembly.rename.sh [assembly] [scaffold_name] # eg nCe.2.0
    
### 5.1 Mask the repeats with [RepeatModeler](http://www.repeatmasker.org/) and [RepeatMasker](http://www.repeatmasker.org/)
```
mkdir RepeatModeler
cd RepeatModeler

#Build the database
BuildDatabase -name [species_name] -engine ncbi ../[renamed_assembly]

# De novo model the repeats
RepeatModeler -engine ncbi -pa 16 -database [species_name]

# Extract repeats from Rhabditida
# queryRepeatDatabase.pl inside RepeatMasker/utils
queryRepeatDatabase.pl -species rhabditida | grep -v "Species:" > Rhabditida.repeatmasker

# Concatenate the 2 files
cat RM*/consensi.fa.classified Rhabditida.repeatmasker > all.repeats

# Mask the assembly
RepeatMasker -lib all.repeats -pa 16 ../[renamed_assembly]
```

### 5.2 Annotation Pipeline if we have RNASeq data

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

#### 5.2.3 Add UTR info to the GFF3 file

Section under testing


### 5.3 Annotation Pipeline if we do not have RNASeq data

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
