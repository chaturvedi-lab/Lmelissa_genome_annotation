# Genome Annotation using MAKER

[MAKER](http://www.yandell-lab.org/software/maker.html) is a great tool for annotating a reference genome using empirical and *ab initio* gene predictions. [GMOD](http://gmod.org/wiki/Main_Page), the umbrella organization that includes MAKER, has some nice tutorials online for running MAKER. However, these were quite simplified examples and it took a bit of effort to wrap my head completely around everything. In addition, there are several pre-processing steps which need to be followed before working on the MAKER pipeline. Here I is my description of a *de novo* genome annotation for *Lycaeides melissa* in detail, so that there is a record and that it is easy to use this as a guide to annotate any genome. In addition, this can be referred to as a reseource for finding great resources for tutorials on installations and running sample datasets. 

## Software & Data

#### Software prerequisites:
1. [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) open-version 1.0.11 and [RepeatMasker](http://www.repeatmasker.org/RMDownload.html) version open-4.0.7 with all dependencies (I used NCBI BLAST version 2.2.28+ but did install all dependencies) and [RepBase](http://www.girinst.org/repbase/) (version used was 20150807).
2. MAKER MPI version 2.31.8 (though any other version 2 releases should be okay).
3. [Augustus](http://bioinf.uni-greifswald.de/augustus/) locally installed on UofU CHPC version 3.3.
4. [BUSCO](http://busco.ezlab.org/) version 2.0.1. (need to update version)
5. [SNAP](http://korflab.ucdavis.edu/software.html) version 2006-07-28.
6. [BEDtools](https://bedtools.readthedocs.io/en/latest/) version 2.17.0.

#### Raw data/resources:
1. `melissa_ref.fasta`: The *de novo* *Lycaeides melissa* reference genome assembled using Dovetail genomics methods. This is in FASTA format.
2. `melissa_Txome_assembly.fasta`: *Needs to be created* A *de novo* transcriptome assembly created using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) and mRNAseq data from 10 *L.melissa* tissues. I still need to build this and get this ready for MAKER. Trinity is pretty easy to run so I won't describe that here. This is in FASTA format.
3. Full protein amino acid sequences, in FASTA format, for three other Lepidopteran species from [NCBI](https://www.ncbi.nlm.nih.gov/genome/) or [Ensembl](http://www.ensembl.org/index.html): *Lep1*, *Lep2*, and *Lep3*.
4. A curated lepidopteran repeat library derived from XX butterfly species (from internal efforts). *Needs to be created*

#### Running commands:
I mostly run test the commands using some small sized sample data and then run the programs on the remote cluster using bash scripts. Whereever possible, I will provide the code for the bash script I used to submit jobs on the cluster.

## Repeat Annotation

#### 1. *De Novo* Repeat Identification
The first, and very important, step to genome annotation is identifying repetitive content. Existing libraries from Repbase or from internal efforts are great, but it is also important to identify repeats *de novo* from your reference genome using `RepeatModeler`. This is pretty easy to do and normally only takes a couple days using 8-12 cores.

```bash
BuildDatabase -name melissa -engine ncbi melissa_mod.fasta
RepeatModeler -pa 8 -engine ncbi -database melissa
```

Here is a sample bash script I used to submit this as a job on the cluster:

~~~
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 180:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J repeatmodeller

#change directory to the folder where your reference fasta file is and where you want to save the output files
cd /uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/repeatmodeller/

#build a database
/uufs/chpc.utah.edu/common/home/u6007910/bin/RepeatModeler-open-1.0.11/BuildDatabase -name melissa -engine ncbi ./melissa_mod.fasta
#run repeatmodeler
/uufs/chpc.utah.edu/common/home/u6007910/bin/RepeatModeler-open-1.0.11/RepeatModeler -pa 8 -engine ncbi -database melissa
~~~

Further steps can be taken to annotate the resulting library, but the most important reason for this library is for downstream gene prediction. This *melissa* library can be combined with several other lepidopteran libraries and annotated.

This script/command will create output folder and files with following extensions or names (substitute melissa with your given name for your database):

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

