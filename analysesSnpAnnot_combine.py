#!/usr/bin/python2.7
import csv
import sys

## ------------------------------------------------- ##
## This python scripts take the snp annotation information from the Annotation file  and combines it with the data file from different analyses. Both these files have scaffold and position for each SNP in common. Make sure to edit the column numbers according to input file.
## Author: Samridhi Chaturvedi
## Copyright: Copyright 2019
## Credits: [""]
## License: GPL3
## Version: 0.0.1
## Maintainer: Samridhi Chaturvedi
## Email: samridhi.chaturvedi@gmail.com
## -------------------------------------------- ##

'''
Usage:
	python2.7 analysesSnpAnnot_combine.py snp_annotation_outtable.txt analyses_snps.txt out_analyses_snpannot.txt 

Prerequisite:
Make sure the analyses file is tab separated. Here is a perl command to replace white space with tab:
	perl -p -i -e 's/ /\t/g' file.txt
'''

##read in the files
annotfile = open(sys.argv[1],'r')
clinefile = open(sys.argv[2], 'r')
outfile = open(sys.argv[3],'w')
outheader_str = "scaffold\tposition\talpha\tbeta\tancestry\tannstart\tannstop\tongene\toncds\tonmRNA\tonexon\tonte\tonprotein\tonmatch\tneargene\tnearcds\tnearmRNA\tnearexon\tnearte\tnearprotein\tnearmatch\tnumneargenes\tIPRnumber\tIPRtext\tGOnumber\tGOtext\n"
outfile.write(outheader_str)
			
scafpos = dict()
next(clinefile)
for line in annotfile:
	line = line.strip('\n')
	line = line.split('\t')
	#print line
	keyannot = line[0]+'-'+line[1]
	scafpos[keyannot] = line

for aline in clinefile:
	aline = aline.strip('\n')
	aline = aline.split('\t')
	keycline = aline[0]+'-'+aline[1]
	clinedict = [int(0)] * 26
	clinedict[0] = aline[0]
	clinedict[1] = aline[1]
	clinedict[2] = aline[2]
	clinedict[3] = aline[3]
	clinedict[4] = aline[4]
	if keycline in scafpos:
		match = scafpos[keycline]
		#print match[3]
		clinedict[5] = match[2]
		clinedict[6] = match[3]
		clinedict[7] = match[4]
		clinedict[8] = match[5]
		clinedict[9] = match[6]
		clinedict[10] = match[7]
		clinedict[11] = match[8]
		clinedict[12] = match[9]
		clinedict[13] = match[10]
		clinedict[14] = match[11]
		clinedict[15] = match[12]
		clinedict[16] = match[13]
		clinedict[17] = match[14]
		clinedict[18] = match[15]
		clinedict[19] = match[16]
		clinedict[20] = match[17]
		clinedict[21] = match[18]
		clinedict[22] = match[19]
		clinedict[23] = match[20]
		clinedict[24] = match[21]
		clinedict[25] = match[22]
		#print clinedict
		out = '\t'.join(clinedict)
		out_str = str(out)
		outfile.write(out_str + "\n")
