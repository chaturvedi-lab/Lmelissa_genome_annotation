#!/usr/bin/python2.7
import re
import argparse
import goatools
from goatools import obo_parser
goont = obo_parser.GODag('./go-basic.obo')

## -------------------------------------------- ##
## Author: Samridhi Chaturvedi
## Copyright: Copyright 2019
## Credits: [""]
## License: GPL3
## Version: 0.0.1
## Maintainer: Samridhi Chaturvedi
## Email: samridhi.chaturvedi@*.com
## Status: dev
## -------------------------------------------- ##

# -------------- #
# Pre-requisites #
# -------------- #
'''
Packages:
 - pip install git+git://github.com/tanghaibao/goatools.git

Data: 
 - wget http://purl.obolibrary.org/obo/go/go-basic.obo

Usage: 
 - python create_snp_annotations.py --map mappos_sub-copy.txt --ann annot_sub-copy.txt --out sub_snp_annotation_table-copy.out
 - python create_snp_annotations.py --help
'''

# --------------------- #
# Global variables conf #
# --------------------- #
''' 
Here we are getting SNPs in 1000bp of each start 
and stop position in the annotation table
'''
snpwindow = 1000			# snp position near or on this window range
snpdict = dict()			# snp pos from mappos_sub.txt file
snppos_tmp = dict()			# Keep map snp pos specific annotation
sep_snppos_tmp = "~,~"
sep_merge = '~'
annotarray = []				# annotation array from annot_sub.txt file
resultarray = []			# final output array
input_mapposition_file = "mappos_sub.txt"
input_annotations_file = "annot_sub.txt"
output_file_name = "sub_snp_annotation_table.out"

# --------------------- #
# Function Definitions  #
# --------------------- #
''' Utility to test integer '''
def is_integer(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

''' Load scaffold position to a dictionary '''
def load_scafpos():
	mapposfile = open(input_mapposition_file,'r')
	for line in mapposfile:
		if line[0] == '#':
			continue
		line = line.split(':')
		scaffold = int(line[0])
		snppos = int(line[1])
		if scaffold in snpdict:
			snpdict[scaffold].append(snppos)
		else:
			snpdict[scaffold] = [snppos]
	return

''' Load annotation table and sort the scaffolds and positions '''
def load_annot():
	annotfile = open(input_annotations_file, 'r')
	for line in annotfile:
		if line[0] == '#':
			continue
		line = line.strip('\n')
		line = line.split('\t')
		if line[2] == 'contig':
			continue
		annotarray.append(line)
	sorted(annotarray, key=lambda x: (x[0], x[3]))
	return

''' Group the position in annot file and write to a temp dict '''
def group_snppos(snppos, line, distance_type):
	if line is not None:
		info_str = line[2] + sep_snppos_tmp + line[3] + sep_snppos_tmp + line[4] + \
		sep_snppos_tmp + distance_type + sep_snppos_tmp + line[8].strip('\n')
	else:
		info_str = None
	if snppos in snppos_tmp:
		if info_str is not None:
			snppos_tmp[snppos].append(info_str)
	else:
		snppos_tmp[snppos] = [info_str]

''' Match position wrt ann start/stop range and cache '''
def match_position(snpscaf, annstart, annstop, line):
	global snpdict
	for snppos in snpdict[snpscaf]:
		if (snppos > annstart and snppos < annstop):
			group_snppos(snppos, line, 'on')
		elif (snppos < annstart and abs(annstart - snppos) < snpwindow):
			group_snppos(snppos, line, 'left')
		elif (snppos > annstop and abs(snppos - annstop) < snpwindow):
			group_snppos(snppos, line, 'right')
		else:
			group_snppos(snppos, None, '')

''' Fetch the IPR and GO terms from the GOA file '''
def fetch_ipr_go_info(ann_info):
	#print ann_info
	ipr_match = re.findall(r'IPR[0-9]+', ann_info) #19,20 (concat IPR using ;)
	ipr_match_str = ';'.join(map(str,ipr_match))
	go_match =  re.findall(r'GO:[0-9]+', ann_info) #21,22 (concat GO using ;)
	go_match_str = ';'.join(map(str,go_match))
	#print "IPR: ", ipr_match_str
	#print "GO : ", go_match_str
	go_term = []
	for g in go_match:
		go_term.append(goont[g].name)
	go_term_str = ';'.join(map(str,go_term))  # Go term for corresponding id
	return ipr_match_str, go_match_str, go_term_str


''' Flatten and merge the annotations with the same snp position '''
def compute_format_result(scaffold, snppos, annotations):
	out_array = [int(0)] * 23
	out_array[0] = scaffold
	out_array[1] = snppos
	for ann in annotations:
		if ann is not None:
			#type,start stop,distance_type,annotation
			ann_array = ann.split(sep_snppos_tmp)
			# Check for ongene only once
			if ann_array[3] == 'on': 
				if ann_array[0] == 'gene':
					# Use the annotation start/stop of 
					out_array[2] = ann_array[1] # annstart
					out_array[3] = ann_array[2] # annstop
					out_array[4] = 1			# ongene
					ann_info = str(ann_array[4])
					ipr_match_str, go_match_str, go_term_str = fetch_ipr_go_info(ann_info)
					out_array[19] = ipr_match_str  	# ipr
					out_array[21] = go_match_str	# go number
					out_array[22] = go_term_str		# go text
				elif ann_array[0] == 'CDS':
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[5] = 1 			# oncds
				elif ann_array[0] == 'mRNA':
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[6] = 1
				elif ann_array[0] == 'exon':
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[7] = 1
				elif ann_array[0] == 'five_prime_UTR' or ann_array[0] == 'three_prime_UTR':
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[8] = 1
				elif ann_array[0] == 'protein_match' or ann_array[0] == 'expressed_sequence_match':
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop		
					out_array[9] = 1
				elif ann_array[0] == 'match' or ann_array[0] == 'match_part' :
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[10] = 1
			# Go here only if immediately previous annotation is near and not on
			elif ann_array[3] == 'left' or ann_array[3] == 'right':
				# Go here only if type matches and is not on for that type
				if ann_array[0] == 'gene' and out_array[4] == 0:
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					# Keep the IPR/GO information for near gene
					ann_info = str(ann_array[4])
					ipr_match_str, go_match_str, go_term_str = fetch_ipr_go_info(ann_info)
					# Merge logic
					if out_array[19] is not 0:
						print "\nAttempt Merge, snp: ",snppos," - ", out_array, '\t-\t', ann_array
						out_array[2] = out_array[2] + sep_merge + ann_array[1] # annstart 
						out_array[3] = out_array[3] + sep_merge + ann_array[2] # annstop 
						out_array[19] = out_array[19] + sep_merge + ipr_match_str # ipr text
						out_array[21] = out_array[21] + sep_merge + go_match_str  # go number
						out_array[22] = out_array[22] + sep_merge + go_term_str 	# go text
					else:
						# Persist the current 
						out_array[19] = ipr_match_str  	# ipr
						out_array[21] = go_match_str	# go number
						out_array[22] = go_term_str		# go text
					out_array[11] = 1 				# neargene
					out_array[18] = out_array[18] + 1
				elif ann_array[0] == 'CDS' and out_array[5] == 0:	  # not onCDS
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[12] = 1 			# oncds
				elif ann_array[0] == 'mRNA' and out_array[6] == 0:
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[13] = 1
				elif ann_array[0] == 'exon' and out_array[7] == 0:
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[14] = 1
				elif (ann_array[0] == 'five_prime_UTR' or ann_array[0] == 'three_prime_UTR') and (out_array[8] == 0):
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[15] = 1
				elif (ann_array[0] == 'protein_match' or ann_array[0] == 'expressed_sequence_match') and (out_array[9] == 0):
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop		
					out_array[16] = 1
				elif (ann_array[0] == 'match' or ann_array[0] == 'match_part') and (out_array[10] == 0) :
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[17] = 1				

	#print "current out ---------->  SNP=", snppos, "\t", out_array
	# check the previous annot array item's position for concatenation
	# checking if type is gene, on & near distance type; only then proceed to concatenation
	resultarray.append(out_array)
	
''' Final computation for each scaffold and write to output resultarray '''
def compute_scafpos_variables(scaffold):
	global snppos_tmp
	#if snppos == 9804201 or snppos == 9804216:
	#print "Matches------> SCF=", scaffold, "\t Ann=",snppos_tmp[9804201],"\t\n"
	# look up snppos_tmp cache
	for key, value in snppos_tmp.iteritems():
		# compute all the results
		compute_format_result(scaffold,key,value)
	snppos_tmp.clear()

''' Search annotations for snp position '''
def search_snpposition_in_annotation():
	previous_scaffold = 0
	line_index = 0
	for line in annotarray:
		scaffold = int(line[0].split('_')[1])
		# if the scaffold changes then process output results
		if previous_scaffold != 0 and previous_scaffold != scaffold:
			compute_scafpos_variables(previous_scaffold)
			print "\nScaffold processed: ", previous_scaffold
			print "------------------------"
		if line[2] == "contig":
			print "\n Scaffold: ", scaffold
			continue
		# process every scaffold in map snp positions
		elif scaffold in snpdict:
			annstart = int(line[3])
			annstop = int(line[4])
			match_position(scaffold, annstart, annstop, line)
		# keep scaffold id state
		previous_scaffold = scaffold
		line_index = line_index + 1
	else:
		compute_scafpos_variables(previous_scaffold)
		print "Completed processing scaffold ID %d", previous_scaffold

#function to write result
def write_result(outfile):
	outheader_str = "scaffold\tposition\tannstart\tannstop\tongene\toncds\tonmRNA\t \
		onexon\tonte\tonprotein\tonmatch\tneargene\tnearcds\tnearmRNA\tnearexon\tnearte\t \
		nearprotein\tnearmatch\tnumneargenes\tIPRnumber\tIPRtext\tGOnumber\tGOtext\n"	
	outfile.write(outheader_str)
	for item in resultarray:
		out_str = '\t'.join([str(x) for x in item]) + '\n'
		outfile.write(out_str)	

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--map", help="Add scaffold position file full path", default="mappos_sub.txt")
	parser.add_argument("--ann", help="Add genome annotations file full path", default="annot_sub.txt")
	parser.add_argument("--out", help="Add snp annotation output file full path", default="sub_snp_annotation_table.out")
	args = parser.parse_args()
	if args.map:  input_mapposition_file = args.map
	if args.ann:  input_annotations_file = args.ann
	if args.out:  output_file_name = args.out
	
	print "Configs: \n\t MAP:", input_mapposition_file, "\n\t ANN:", input_annotations_file, "\n\t OUT:", output_file_name, "\n" 
	# Load the dataset
	load_scafpos()
	print "Loaded scaffold data: ", len(snpdict)
	load_annot()
	#for a in annotarray:
	#	print a[0], '\t' ,a[3], '\t' ,a[4], '\t',a[2]
	print "Loaded annotation data: ", len(annotarray)
	# compute the positions
	search_snpposition_in_annotation()

	
	# output file descriptor 
	outfile = open(output_file_name, 'w')
	# write results 
	write_result(outfile)
