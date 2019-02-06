Following genome annotation using MAKER,we need to use the information to do SNP annotation for the dataset used in actual analyses. To do this I created custom python scripts and followed the following steps:

** Step 1: Use the create_snp_annotations.py script to make the snp_annotation file **
I did SNP annotation on my computer using the python script (*create_snp_annotations.py*). This script created the master SNP annotation file which contains functional and structural annotation information. I used this script to do the downstream statistical analyses (mainly randomizations).

Given a snp list file with scaffold and positions, and a genome annotation file (from MAKER), this script creates a SNP annotation file. Before giving the genome annotation file from maker (maker_ipr_go.txt, ), I removed the lines which had sequence data in this file and created the file genome_annotation.txt (lines=2310823). I used the mappos.txt (39193 SNPs) file as a list of SNPs with scaffold and position in the alignment from the hybrid zone data. 

Usage of python script: 
```
python create_snp_annotations.py --map mappos.txt --ann genome_annotation.txt --out snp_annotation_outtable.txt

```
I then moved the input, output files, script and dummy data to the folder: /uufs/chpc.utah.edu/common/home/gompert-group1/data/lycaeides/dovetail_melissa_genome/Annotation/maker/melissa_round2.maker.output/snp_annotation

** Step 2: Combine the Clines, popanc files with cutoff AIMS snps with the SNP annotation file. **

I first combined the clines and popanc files to create one file with all the hybrid zone analyses information.

In R 
```
## clines files #####
pct1<-read.table("clines/clines_pct1.txt", header=T)
pct2<-read.table("clines/clines_pct2.txt", header=T)
pct3<-read.table("clines/clines_pct3.txt", header=T)
### popanc files ####
aims1<-read.table("popanc/ancestrymeans_aims1.txt", header=T)
aims2<-read.table("popanc/ancestrymeans_aims2.txt", header=T)
aims3<-read.table("popanc/ancestrymeans_aims3.txt", header=T)
### combine the files ###
comb1<-cbind(pct1,aims1$ancestry.abs.)
comb2<-cbind(pct2,aims2$ancestry.abs.)
comb3<-cbind(pct3,aims3$ancestry.abs.)
#rename columns
colnames(comb1)[colnames(comb1) == "aims1$ancestry.abs."]<-"ancmeans"
colnames(comb2)[colnames(comb2) == "aims2$ancestry.abs."]<-"ancmeans"
colnames(comb3)[colnames(comb3) == "aims3$ancestry.abs."]<-"ancmeans"
#write out the files 
write.table(comb1, file = "clines_popanc1.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE) write.table(comb2, file = "clines_popanc2.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE) write.table(comb3, file = "clines_popanc3.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
```
On command line:

Use python script to analysesSnpAnnot_combine.py to combine the hybrid zone analyses file from above with annotation file.

Script usage:

```
python analysesSnpAnnot_combine.py  snp_annotation_outtable.txt analyses_file.txt outfile.txt
##Commands I ran:
### aims 1 ####
python analysesSnpAnnot_combine.py snp_annotation_outtable.txt clines_popanc1.txt comb1_clines_popanc_annot.txt
### aims 2 ####
python analysesSnpAnnot_combine.py snp_annotation_outtable.txt clines_popanc2.txt comb2_clines_popanc_annot.txt
### aims 3 ####
python analysesSnpAnnot_combine.py snp_annotation_outtable.txt clines_popanc3.txt comb3_clines_popanc_annot.txt

```
