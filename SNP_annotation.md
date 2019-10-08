Here SNPs are single nucleotide polymorphisms or genomic variants called using a specific genome alignment and variant calling pipeline.

Following genome annotation using MAKER, we can use the annotation information to annotate SNPs used in specific studies and know their functional properties. The following scripts can be used to annotate a SNP data set which has SNPs scaffold and position and the SNPs need to be aligned to the same genome whose annotation we are using.

** Step 1: Use the create_snp_annotations.py script to make the snp_annotation file **
I did SNP annotation on my computer using the python script (*create_snp_annotations.py*) present in this directory. This script created the master SNP annotation file which contains functional and structural annotation information. I used this script to do the downstream statistical analyses (mainly randomizations).

Given a SNP list file with scaffold and positions, and a genome annotation file (from MAKER), this script creates a SNP annotation file. Before giving the genome annotation file from maker (maker_ipr_go.txt, ), I removed the lines which had sequence data in this file and created the file genome_annotation.txt (lines=2310823). I used the mappos.txt (39193 SNPs) file as a list of SNPs with scaffold and position in the alignment from the hybrid zone data. 

Usage of python script: 

```
python create_snp_annotations.py --map mappos.txt --ann genome_annotation.txt --out snp_annotation_outtable.txt

```

** Step 2: Combine the data from specific statistical analyses with the SNP annotation file to find important SNPs **

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
I then used the analysesSNpAnnot_combine.py file (in this directory) to combine the SNP annotation and data files.

Usage:

```bash
python analysesSnpAnnot_combine.py  snp_annotation_outtable.txt analyses_file.txt outfile.txt
```
