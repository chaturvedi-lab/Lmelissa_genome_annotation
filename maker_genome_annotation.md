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

##### 3. MAKER Analyses: ROUND 1

**INSTALLATION**
I asked Anita to install MAKER globally so that I could have all dependencies running. I am using MAKER version 2.31.10. I first tested if MAKER is running by testing on example data which is in the following directory /uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/maker/tem/data

As you will see this contains the same files as the data/ directory that comes with MAKER. I copied them here for convenience.

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

This creates three files (type ls to see)

MAKER is pretty easy to get going and relies on properly completed control files. These can be generated by issuing the command maker -CTL. The only control file we will be the maker_opts.ctl file. In this first round, we will obviously providing the data files for the repeat annotation (rm_gff), the transcriptome assembly (est), and for the other Squamate protein sequences (protein). We will also set the model_org to 'simple' so that only simple repeats are annotated (along with RepeatRunner). Here is the full control file for reference.

























