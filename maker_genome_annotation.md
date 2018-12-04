# Genome Annotation using MAKER

[MAKER](http://www.yandell-lab.org/software/maker.html) is a great tool for annotating a reference genome using empirical and *ab initio* gene predictions. [GMOD](http://gmod.org/wiki/Main_Page), the umbrella organization that includes MAKER, has some nice tutorials online for running MAKER. However, these were quite simplified examples and it took a bit of effort to wrap my head completely around everything. In addition, there are several pre-processing steps which need to be followed before working on the MAKER pipeline. Here I is my description of a *de novo* genome annotation for *Lycaeides melissa* in detail, so that there is a record and that it is easy to use this as a guide to annotate any genome. In addition, this can be referred to as a reseource for finding great resources for tutorials on installations and running sample datasets. 

## Software & Data

#### Software prerequisites:
1. [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) open-version 1.0.11 and [RepeatMasker](http://www.repeatmasker.org/RMDownload.html) version open-4.0.7 with all dependencies (I used NCBI BLAST version 2.2.28+ but did install all dependencies) and [RepBase](http://www.girinst.org/repbase/) (version used was 20150807).
1. [RepeatScout](http://bix.ucsd.edu/repeatscout/) Version 1.0.5 and [RepeatMasker](http://www.repeatmasker.org/RMDownload.html) version open-4.0.7 with all dependencies (I used NCBI BLAST version 2.2.28+ but did install all dependencies) and [RepBase](http://www.girinst.org/repbase/) (version used was 20150807).
2. MAKER MPI version 2.31.8 (though any other version 2 releases should be okay).
3. [Augustus](http://bioinf.uni-greifswald.de/augustus/) locally installed on UofU CHPC version 3.3.
4. [BUSCO](http://busco.ezlab.org/) version 2.0.1. (need to update version)
5. [SNAP](http://korflab.ucdavis.edu/software.html) version 2006-07-28.
6. [BEDtools](https://bedtools.readthedocs.io/en/latest/) version 2.17.0.

**RepeatScout was already present in my RepeatModeler folder since it is one of the prerequisite for this program. 

#### Raw data/resources:
1. `melissa_ref.fasta`: The *de novo* *Lycaeides melissa* reference genome assembled using Dovetail genomics methods. This is in FASTA format.
2. `melissa_Txome_assembly.fasta`: *Needs to be created* A *de novo* transcriptome assembly created using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) and mRNAseq data from 10 *L.melissa* tissues. I still need to build this and get this ready for MAKER. Trinity is pretty easy to run so I won't describe that here. This is in FASTA format.
2. `melissa_Txome_assembly.fasta`: *Needs to be created* A *de novo* transcriptome assembly created using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) and RNAseq data from 24 *L.melissa* tissues (paired-end libraries). Trinity is pretty easy to run so I won't describe that here. This is in FASTA format. I also created a genome reference based transcriptome assembly using TopHat.
3. Full protein amino acid sequences, in FASTA format, for three other Lepidopteran species from [NCBI](https://www.ncbi.nlm.nih.gov/genome/) or [Ensembl](http://www.ensembl.org/index.html): *Lep1*, *Lep2*, and *Lep3*.
4. A curated lepidopteran repeat library derived from XX butterfly species (from internal efforts). *Needs to be created*

#### Running commands:
I mostly run test the commands using some small sized sample data and then run the programs on the remote cluster using bash scripts. Whereever possible, I will provide the code for the bash script I used to submit jobs on the cluster.
I mostly test the commands using some small sized sample data and then run the programs on the remote cluster using bash scripts. Whereever possible, I will provide the code for the bash script I used to submit jobs on the cluster.

## Repeat Annotation

#### 1. *De Novo* Repeat Identification
The first, and very important, step to genome annotation is identifying repetitive content. Existing libraries from Repbase or from internal efforts are great, but it is also important to identify repeats *de novo* from your reference genome using `RepeatModeler`. This is pretty easy to do and normally only takes a couple days using 8-12 cores.
The first, and very important, step to genome annotation is identifying repetitive content. Existing libraries from Repbase or from internal efforts are great, but it is also important to identify repeats *de novo* from your reference genome using `RepeatScout`. This is pretty easy to do and normally only takes a couple days using 8-12 cores. I decided to use RepeatScout instead of RepeatModeler since this seemed to be a better method and the functioning of this program is more clear than RepeatModeler. 

Here are the steps I used to run RepeatScout:

**Step 1: Build a lmer table**
First, I used build_lmer_table from the RepeatScout program to create a tabulated file which contains the frequency of all l-mers in the sequence to be analyzed. Here is the bash code:

*Script is called: runRepeatScout_buildlmer.sh*

~~~
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 300:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J repscout

#sequence = genome file
#freq = output file name

cd workingdirectory

```bash
BuildDatabase -name melissa -engine ncbi melissa_mod.fasta
RepeatModeler -pa 8 -engine ncbi -database melissa
```
/uufs/chpc.utah.edu/common/home/u6007910/bin/RepeatModeler-open-1.0.11/RepeatScout-1/build_lmer_table -sequence melissa_ref.fasta -freq melissa_lmer.freq
~~~

The output file **melissa_lmer.freq** is used in the next step to run RepeatScout.

Here is a sample bash script I used to submit this as a job on the cluster:
**Step 2: Running RepeatScout**
Then I ran RepeatScout. RepeatScout takes the lmer table from previous step and the genome sequence and produces a fasta file that contains all the repetitive elements that it could find. Here is the bash code:

*Script is called: runRepeatScout.sh*

~~~
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 180:00:00
#SBATCH -t 300:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J repeatmodeller
#SBATCH -J repscout

#sequence = your genome file
#freq = lmer freq file from previous step
#output = name of output file

#change directory to the folder where your reference fasta file is and where you want to save the output files
cd /uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/repeatmodeller/
cd workingdirectory

#build a database
/uufs/chpc.utah.edu/common/home/u6007910/bin/RepeatModeler-open-1.0.11/BuildDatabase -name melissa -engine ncbi ./melissa_mod.fasta
#run repeatmodeler
/uufs/chpc.utah.edu/common/home/u6007910/bin/RepeatModeler-open-1.0.11/RepeatModeler -pa 8 -engine ncbi -database melissa
/uufs/chpc.utah.edu/common/home/u6007910/bin/RepeatModeler-open-1.0.11/RepeatScout-1/RepeatScout -sequence melissa_ref.fasta -output melissa.rs -freq melissa_lmer.freq
~~~

Further steps can be taken to annotate the resulting library, but the most important reason for this library is for downstream gene prediction. This *melissa* library can be combined with several other lepidopteran libraries and annotated.
The output file **melissa.rs** will be used in the next filtering step.

**Step 3: First filtering**
In this step, the "filter-stage-1.prl" script is run on the output of RepeatScout to remove low-complexity and tandem elements; RepeatMasker is run on the sequence of interest using this filtered RepeatScout library. Here is the bash code:

*Script is called: runRepeatScout_filter1.sh*

~~~
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 300:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J repscout

This script/command will create output folder and files with following extensions or names (substitute melissa with your given name for your database):
cd workingdirectory

/uufs/chpc.utah.edu/common/home/u6007910/bin/RepeatModeler-open-1.0.11/RepeatScout-1/filter-stage-1.prl  melissa.rs > melissa_filter1.fa
~~~

The output file **melissa_filter1.fa** is used run RepeatMasker to mainly count the number of times each repeat appears in the filtered library from RepeatScout.

**Step 3: RepeatMasker run 1**
I ran RepeatMasker using the output file from the previous step as the custom repeat library. ( I did this step in the repeatmasker library). Here is the bash code I used:

*Script is called: runRepeatMask_filter2.sh*

~~~
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 300:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J masking

#pa = number of threads
#lib = repeatscout filtered repeat file from previous step
#dir = output directory
#genome file = melissa_ref.fasta

cd workingdirectory

/uufs/chpc.utah.edu/common/home/u6007910/bin/RepeatMasker/RepeatMasker -pa 8 -lib melissa_filter1.fa -dir ./outputdir melissa_ref.fasta 
~~~

The output file **melissa_ref.fasta.out** is used run the second filtering step below.

**Step 4: Filter out more repeats**
In this step, I used the filter-stage-2.prl to filter out any repeats that appear less than 10 times (default setting). Here is the bash code:

*Script is called: runRepeatMask_filter2.sh*

~~~
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 300:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J repscout

cd workingdirectory

cat melissa_filter1.fa | /uufs/chpc.utah.edu/common/home/u6007910/bin/RepeatModeler-open-1.0.11/RepeatScout-1/filter-stage-2.prl  --cat=melissa_ref.fasta.out > melissa_filter2.fa
~~~

**Step 5: Repeatmasker run 2**
I ran repeatmasker with the new repeat library from the previous step. Here is the bash code:

*Script is called: runRepeatMask_filter2.sh*

~~~
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 300:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J masking

#module load perl/5.18.1

cd workingdirectory

/uufs/chpc.utah.edu/common/home/u6007910/bin/RepeatMasker/RepeatMasker -pa 8 -lib melissa_filter2.fa -dir ./outdir melissa_ref.fasta 
~~~

**BuildDatabase creates these files**
1. melissa_mod.fasta
2. melissa.nhr
3. melissa.nin
4. melissa.nnd
5. melissa.nni
6. melissa.noq
7. melissa.nsq
8. melissa.translation
9. unaligned.fa

**RepeatModeler creates this folder and files**
A folder named RM_ will be created in your output directory. This folder will contain 5 subfolders labelled round1-round5. These folders will contain several output files. The main file which is the repeat library and should be used for downstream stuff is called consensi.fa.classified. This file will have a list of *identified* and *unknown repeats* found in the genome.
A folder named RM_ will be created in your output directory. This folder will contain 5 subfolders labelled round1-round5. These folders will contain several output files. The main file which is the repeat library and should be used for downstream stuff is called *consensi.fa.classified*. This file will have a list of *identified* and *unknown repeats* found in the genome.

#### 2. Repeat masking or Full Repeat Annotation
I used the *de novo* repeat library from RepeatModeler to mask repeats to avoid issues with RepBase.

Here is the bash script I used for this:
Before I do repeatmasking, I wanted to filter proteins or coding regions from the repeatmodeler repeat library. This can be done in various ways. One of these ways is to use a proteome or transcriptome to identify putative coding regions which have been identified as repeats by repeatmodeler. Therefore, we can filter out these hits from the repeatmodeler library before we do repeatmasking so that we encounter lesser errors during genome annotation.

There are two approaches I used to filter coding regions or proteins from the repeatmodeler library. I mainly wanted to compare the proportion of reads identified as repeats using both these approaches.

**APPROACH 1**
This approach is described here: https://media.readthedocs.org/pdf/blaxter-lab-documentation/latest/blaxter-lab-documentation.pdf
Here, I mainly filter repeats by first blasting them against the default RepeatMasker TE database.

The folder in which all these files are saved on the cluster is: */uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/repeatmasker/filtering_blast*

Here are the steps I used alongwith the bash code for each of these steps:

A. Blast proteome against RepeatMasker TE database

##### 3. MAKER Analyses

**INSTALLATION**
I asked Anita to install MAKER globally so that I could have all dependencies running. I am using MAKER version 2.31.10. MAKER is pretty easy to get going and relies on properly completed control files. These can be generated by issuing the command maker -CTL. The only control file we will be the maker_opts.ctl file. In this first round, we will obviously providing the data files for the repeat annotation (rm_gff), the transcriptome assembly (est), and for the other Squamate protein sequences (protein). We will also set the model_org to 'simple' so that only simple repeats are annotated (along with RepeatRunner). Here is the full control file for reference.

I first tested if MAKER is running by testing on example data which is in the following directory /uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/maker/tem/data

This contains the same files as the data/ directory that comes with MAKER. I copied them here for convenience.

~~~
ls -l 
dpp_contig.fasta
dpp_contig.maker.output
dpp_est.fasta
dpp_protein.fasta
~~~

Next we need to tell MAKER all the details about how we want the annotation process to proceed. Because there can many variables and options involved in annotation, command line options will be too cumbersome. Instead MAKER uses a set of configuration files which guide each run. Before each MAKER run, we need to create a set of generic configuration files in the current wokring directory by typing the following:

~~~
module load MAKER
maker -CTL
~~~

This creates three files (type ls to see):

* *maker_exe.ctl* contains the path information for the underlying executables
* *maker_bopt.ctl* contains filtering statistics for BLAST and Exonerate (I did not modify this file)
* *maker_opt.ctl* contains all other information for MAKER, including the location of the input genome file.

Control files are run specific and a separate set of control files will need to be generated for each genome annotated with MAKER. MAKER will look for control files in the current working directory, so it is recommended that MAKER be run in a separate directory containing unique control files for each genome.

## 2. Round 1 of MAKER

Here is the modified maker_opts.ctl file for running the first round of MAKER. I did not use any AUGUSTUS information for this round. I just ran this will default parameters by only providing the location of the genome, the transcriptome and the protein sequences from LepBase.

~~~
 cat maker_opts.ctl 
#-----Genome (these are always required)
genome=/uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/repeatmasker/melissa_ref.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=/uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/matt_transcriptome/txome_assembly/trinity_denovo/trinity_out_dir/bowtieassembly/Trinity.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=/uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/lepbase/lepbase_protein_fasta/lepbase_prots.fa  #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=all #select a model organism for RepBase masking in RepeatMasker
rmlib=/uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/repeatscout/melissa_filter2.fa #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/uufs/chpc.utah.edu/sys/installdir/maker/2.31.10/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=#pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
~~~

Here is the bash script I used to run MAKER on the cluster (This round took roughly 13 days so be cautious of the time)

~~~
 cat runRound1_maker.sh 
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 320:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J maker

module load maker

cd /uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/maker/

maker -fix_nucleotides -base melissa_round1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl
~~~

This will create the folder *melissa_round1.maker.output* and will write the output in this folder. The name of the output directory is based on the input genomic sequence file.

Now let's see what is inside the output directory:

~~~
cd melissa_round1.maker.output
ls
~~~

Here are the files created by MAKER:

~~~
maker_bopts.log
maker_exe.log
maker_opts.log
melissa_mod.log
melissa_mod.log.all.gff
melissa_mod.log.noseq.gff
melissa_round1_datastore (folder)
melissa_round1.db
melissa_round1_master_datastore_index.log
mpi_blastdb (folder)
seen.dbm
~~~

* The maker_opts.log, maker_exe.log, and maker_bopts.log files are logs of the control files used for this run of MAKER.
* The mpi_blastdb directory contains FASTA indexes and BLAST database files created from the input EST, protein, and repeat databases.
* The melissa_round1_master_datastore_index.log contains information on both the run status of individual contigs and information on where individual contig data is stored.
* The melissa_round1_datastore (folder) directory contains a set of subfolders, each containing the final MAKER output for individual contigs from the genomic fasta file.

Once a MAKER run is finished the most important file to look at is the melissa_round1_master_datastore_index.log to see if there were any failures.

~~~
less -S melissa_round1_master_datastore_index.log
~~~

If everything proceeded correctly you should see the following:

~~~
contig-dpp-500-500      dpp_contig_datastore/05/1F/contig-dpp-500-500/  STARTED
contig-dpp-500-500      dpp_contig_datastore/05/1F/contig-dpp-500-500/  FINISHED
~~~

There are only entries describing a single contig because there was only one contig in the example file. These lines indicate that the contig contig-dpp-500-500 STARTED and then FINISHED without incident. Other possible entries include:

* FAILED - indicates a failed run on this contig, MAKER will retry these
* RETRY - indicates that MAKER is retrying a contig that failed
* SKIPPED_SMALL - indicates the contig was too short to annotate (minimum contig length is specified in maker_opt.ctl)
* DIED_SKIPPED_PERMANENT - indicates a failed contig that MAKER will not attempt to retry (number of times to retry a contig is specified in maker_opt.ctl)

The entries in the melissa_round1_master_datastore_index.log file also indicate that the output files for this contig are stored in the directory specified for each scaffold. For eg., for scaffold 1 the files are stored in melissa_round1_datastore/EA/C3/Scaffold_1/. Knowing where the output is stored may seem trivial, however, input genome fasta files can contain thousands even hundreds-of-thousands of contigs, and many file-systems have performance problems with large numbers of sub-directories and files within a single directory. Even when the underlying file-system handles things gracefully, access via network file-systems can still be an issue. To deal with this problem, MAKER creates a hierarchy of nested sub-directory layers, starting from a 'base', and places the results for a given contig within these datastore of possibly thousands of nested directories. The master_datastore_index.log file this is essential for identifying where the output for a given contig is stored.

Now let's take a look at what MAKER produced for Scaffold_1 in 'Scaffold_1'.

~~~
cd melissa_round1_datastore/EA/C3/Scaffold_1/
ls -1
~~~
The directory should contain a number of files and a directory.

~~~
Scaffold_1.gff
Scaffold_1.maker.proteins.fasta
Scaffold_1.maker.transcripts.fasta
Scaffold_1.maker.trnascan.transcripts.fasta
run.log
theVoid.Scaffold_1
~~~

* The Scaffold_1.gff contains all annotations and evidence alignments in GFF3 format. This is the important file for use with Apollo or GBrowse. 
* The Scaffold_1.maker.transcripts.fasta and Scaffold_1.maker.proteins.fasta files contain the transcript and protein sequences for MAKER produced gene annotations. 
* The run.log file is a log file. If you change settings and rerun MAKER on the same dataset, or if you are running a job on an entire genome and the system fails, this file lets MAKER know what analyses need to be deleted, rerun, or can be carried over from a previous run. One advantage of this is that rerunning MAKER is extremely fast, and your runs are virtually immune to all system failures.
* The directory theVoid.Scaffold_1 contains raw output files from all the programs MAKER runs (Blast, SNAP, RepeatMasker, etc.). You can usually ignore this directory and its contents.

The datastore directory contains one set of output files for each contig/chromosome from the input assembly, but at some point you're going to want merged files containing all of your output (i.e. a single GFF3 and FASTA file containing all genes). To do this we use two accessory scripts that come with MAKER: gff3_merge and fasta_merge. Both take the master_datastore_index.log file as input.

Before proceeding to this, I copied the melissa_round1_master_datastore_index.log file to melissa_mod.log file and changed the directory of scaffold 1631 & 1638 and changed their status to finished. This is because I ran MAKER seperately on these two scaffolds as MAKER failed to run on these scaffolds. These folders are called 1631 and 1638 respectively and are present in the melissa_round1_datastore folder, alongwith all other scaffolds files. 

## 3. Post-ROUND1 of MAKER

So let's move into the directory containing that file and run those scripts.

~~~
cd melissa_round1.maker.output

fasta_merge -d melissa_mod.log 
gff3_merge -d melissa_mod.log 
~~~

Now you will se a number of new files that represent the merged output for the entire assembly.

## 4. Training Gene Prediction Software

Besides mapping the empirical transcript and protein evidence to the reference genome and repeat annotation (not much of this in our example, given we've done so much up front), the most important product of this MAKER run is the gene models. These are what is used for training gene prediction software like augustus and snap.

*SNAP*

SNAP is pretty quick and easy to train. Issuing the following commands will perform the training. It is best to put some thought into what kind of gene models you use from MAKER. In this case, we use models with an AED of 0.25 or better and a length of 50 or more amino acids, which helps get rid of junky models. This is done in the folder: /uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/maker

Here are the basic steps:

~~~
mkdir snap
mkdir snap/round1
cd snap/round1

# export 'confident' gene models from MAKER and rename to something meaningful
#trying with -x 0.25
maker2zff -x 0.25 -l 50 ../../melissa_round1.maker.output/melissa_mod.log.all.gff

#gather some stats and validate using fathom 
fathom genome.ann genome.dna -gene-stats > gene-stats.log
fathom genome.ann genome.dna -validate > validate.log 2>&1
# collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom genome.ann genome.dna -categorize 1000 > categorize.log #creates alt, err, old, uni and wrn files
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log #creates export files
# create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
#assemble the HMM
hmm-assembler.pl melissa params > melissa_r1_length50_aed0.25.hmm
~~~

*Augustus*

I used the augustus training output from the BUSCO insecta run. I copied files from BUSCO folder to maker/augustus folder:

~~~
cp ../../busco/run_insecta/augustus_output/retraining_parameters/* ./
#change name of all the files
for file in * ; do mv $file ${file//BUSCO_insecta_1988554456/lycaeides_melissa}; done

awk '{ if ($2 == "est2genome") print $0}' melissa_round1.maker.output/melissa_mod.log.noseq.gff > melissa_r1_maker_est2genome.gff
awk '{ if ($2 == "protein2genome") print $0}' melissa_round1.maker.output/melissa_mod.log.noseq.gff > melissa_r1_maker_protein2genome.gff
awk '{ if ($2 ~ "repeat") print $0}' melissa_round1.maker.output/melissa_mod.log.noseq.gff > melissa_r1_maker_repeats.gff
#move files to the main maker folder
mv melissa_r1_maker*.gff ../../
~~~

Finally, I had to copy the melissa specific files to the augustus config folder so that augustus can use these custom species details through MAKER. I asked Anita to create a lycaeides_melissa folder in the config/species folder for augustus and then asked her to copy the retraining files to the folder: /uufs/chpc.utah.edu/sys/installdir/augustus/3.3/config/species/lycaeides_melissa

## 5. Round 2 of MAKER
Then I made the following changes:
* In the maker_opts_round2.ctl file: change augustus_species to lycaeides_melissa
* In maker_exe.ctl file: add the path to the agusutus executable: /uufs/chpc.utah.edu/sys/installdir/augustus/3.3/src/augustus
* In the bash script to submit to the cluster ad the config path: export AUGUSTUS_CONFIG_PATH="/uufs/chpc.utah.edu/sys/installdir/augustus/3.3/config"

Here are the control files for MAKER Round 2:

~~~
cat maker_opts_round2.ctl

#-----Genome (these are always required)
genome=/uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/repeatmasker/melissa_ref.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=/uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/maker/melissa_r1_maker_est2genome.gff #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=  #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=/uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/maker/melissa_r1_maker_protein2genome.gff  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=/uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/maker/melissa_r1_maker_repeats.gff #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=/uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/maker/snap/round1/melissa_r1_length50_aed0.25.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=lycaeides_melissa #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
~~~

Here is the maker_exe.ctl file:

~~~
cat maker_exe.ctl 
#-----Location of Executables Used by MAKER/EVALUATOR
makeblastdb=/uufs/chpc.utah.edu/sys/installdir/blast/ncbi-blast-2.7.1+/bin/makeblastdb #location of NCBI+ makeblastdb executable
blastn=/uufs/chpc.utah.edu/sys/installdir/blast/ncbi-blast-2.7.1+/bin/blastn #location of NCBI+ blastn executable
blastx=/uufs/chpc.utah.edu/sys/installdir/blast/ncbi-blast-2.7.1+/bin/blastx #location of NCBI+ blastx executable
tblastx=/uufs/chpc.utah.edu/sys/installdir/blast/ncbi-blast-2.7.1+/bin/tblastx #location of NCBI+ tblastx executable
formatdb= #location of NCBI formatdb executable
blastall= #location of NCBI blastall executable
xdformat= #location of WUBLAST xdformat executable
blasta= #location of WUBLAST blasta executable
RepeatMasker=/uufs/chpc.utah.edu/sys/installdir/repeatmasker/4.0.7/RepeatMasker/RepeatMasker #location of RepeatMasker executable
exonerate=/uufs/chpc.utah.edu/sys/pkg/exonerate/2.2.0/bin/exonerate #location of exonerate executable

#-----Ab-initio Gene Prediction Algorithms
snap=/uufs/chpc.utah.edu/sys/pkg/SNAP/snap/snap #location of snap executable
gmhmme3= #location of eukaryotic genemark executable
gmhmmp= #location of prokaryotic genemark executable
augustus=/uufs/chpc.utah.edu/sys/installdir/augustus/3.3/src/augustus #location of augustus executable
fgenesh= #location of fgenesh executable
tRNAscan-SE= #location of trnascan executable
snoscan= #location of snoscan executable

#-----Other Algorithms
probuild= #location of probuild executable (required for genemark)
~~~

Here is the bash script to run the round 2 of MAKER:

~~~
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 320:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J maker

module load maker
module load augustus/3.3
module load blast/2.7.1
module load hmmer3

export AUGUSTUS_CONFIG_PATH="/uufs/chpc.utah.edu/sys/installdir/augustus/3.3/config"

cd /uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/maker/

maker -fix_nucleotides -base melissa_round2 maker_opts_round2.ctl maker_bopts.ctl maker_exe.ctl
~~~

## 6. Post-round2 of MAKER

#generate an id mapping file using maker_map_ids

~~~
maker_map_ids --prefix melissa_ melissa_round2.all.gff > melissa_round2.all.map
~~~

This creates a two-column tab-delimited file with the original id in column 1 and the new
id in column 2. The --prefix is where you give your registered genome prefix; the value
following --justify determines the length of the number following the prefix (make
sure that you allow adequate places for the number of genes in the annotation set, e.g., if
you have 10,000 genes, --justify should be set to at least 5.

#use map file to change ids in gff3 and fasta file

~~~
cp melissa_round2.all.gff melissa_round2.all.ids.gff 
map_gff_ids melissa_round2.all.map melissa_round2.all.ids.gff 
map_fasta_ids melissa_round2.all.map melissa_round2.all.maker.proteins.fasta
map_fasta_ids melissa_round2.all.map melissa_round2.all.maker.transcripts.fasta 
~~~

#assigning putative gene function using maker and NCBI BLAST+

~~~
mkdir uniprot
#download the uniprot file from website (http://www.uniprot.org) to the desktop.Then copy it to the cluster.
scp ./Downloads/uniprot_sprot.fasta.gz ssh u6007910@kingspeak.chpc.utah.edu:/uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/maker/melissa_round2.maker.output/uniprot
~~~

1. Index the UniProt/Swiss-Prot multi-FASTA file using makeblastdb:

~~~
makeblastdb -in uniprot/uniprot_sprot.fasta -input_type fasta -dbtype prot
~~~

This creates 3 files:
* uniprot_sprot.fasta.phr 
* uniprot_sprot.fasta.pin
* uniprot_sprot.fasta.psq

2. BLAST the MAKER-generated protein FASTA file to UniProt/SwissProt with
BLASTP.

~~~
blastp -db uniprot/uniprot_sprot.fasta -query melissa_round2.all.maker.proteins.fasta -out maker2uni.blastp -evalue .000001 -outfmt 6 -num_alignments 1 -seg yes -soft_masking true -lcase_masking
~~~

Used the following bash script for this:

~~~
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 300:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J maker

module load maker

cd /uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/maker/melissa_round2.maker.output

blastp -db uniprot/uniprot_sprot.fasta -query melissa_round2.all.maker.proteins.fasta -out maker2uni.blastp -evalue .000001 -outfmt 6 -num_alignments 1 -seg yes -soft_masking true -lcase_masking
~~~

The key parts of this BLAST command line include the specification of the tabular format (-outfmt 6), and the -num_alignments 1.

3. Add protein homology data to the MAKER GFF3 and FASTA files 

~~~
maker_functional_gff uniprot/uniprot_sprot.fasta maker2uni.blastp melissa_round2.all.ids.gff > melissa_functional_blast.gff
maker_functional_fasta uniprot/uniprot_sprot.fasta maker2uni.blastp melissa_round2.all.maker.proteins.fasta > melissa_proteins_functional_blast.fasta
~~~

## INTERPROSCAN: Annotating using Interproscan

1. Download and extract interproscan (https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload) Also look at the how to install link on this page to see installations requirement.

~~~
mkdir interproscan
cd interproscan
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.32-71.0/interproscan-5.32-71.0-64-bit.tar.gz
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.32-71.0/interproscan-5.32-71.0-64-bit.tar.gz.md5
# Recommended checksum to confirm the download was successful:
md5sum -c interproscan-5.32-71.0-64-bit.tar.gz.md5
# Must return *interproscan-5.32-71.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.
~~~

Extract the tar ball:

~~~
tar -pxvzf interproscan-5.32-71.0-*-bit.tar.gz

# where:
#     p = preserve the file permissions
#     x = extract files from an archive
#     v = verbosely list the files processed
#     z = filter the archive through gzip
#     f = use archive file
~~~

2. Installing Panther Models

InterProScan 5 includes the Panther member database analysis.

Before Installing Panther Data, first ensure you have extracted the distribution of InterProScan 5

The path where this is extracted will be referred to below as [InterProScan5 home].

*Download the Panther model data:*
Download the latest Panther data files from the FTP site into the [InterProScan5 home]/data/ directory:

~~~
cd [InterProScan5 home]/data/

wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-12.0.tar.gz
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-12.0.tar.gz.md5

md5sum -c panther-data-12.0.tar.gz.md5
# This must return *panther-data-12.0.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.
~~~

Extract the Panther data files to the required location:

~~~
tar -pxvzf panther-data-12.0.tar.gz
~~~
3. Running interproscan

Testing if interproscan is running:

~~~
module load python3
cd interproscan-5.32-71.0
./interproscan.sh -i test_proteins.fasta -f tsv
./interproscan.sh -i test_proteins.fasta -f tsv -dp
~~~

The first test should create an output file with the default file name test_proteins.fasta.tsv, and the second would then create test_proteins.fasta_1.tsv (since the default filename already exists).





















