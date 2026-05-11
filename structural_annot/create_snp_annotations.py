#!/usr/bin/python3
import re
import sys
import copy
import argparse

import pandas as pd
from tqdm import tqdm
import goatools
from goatools import obo_parser
goont = obo_parser.GODag('./go-basic.obo')
import datetime


## -------------------------------------------- ##
## Author: Samridhi Chaturvedi
## Copyright: Copyright 2023
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
Repository:
 - https://github.com/karwaan/lyc_rnaseq_transcript_annotations

Packages:
 - pip install git+git://github.com/tanghaibao/goatools.git

Data: 
 - wget http://purl.obolibrary.org/obo/go/go-basic.obo

Usage: 
 - python create_snp_annotations.py --pos mappos_sub.txt --ann annot1631_sub.txt --out sub_snp_annotations_table_1631.out
 - python create_transcript_annotations.py --help
'''

# --------------------- #
# Global variables conf #
# --------------------- #
''' 
Here we are getting SNPs in 1000bp of each start 
and stop position in the annotation table
'''
snp_window = 10000			# snp position near or on this window range
sep_merge = '~'

input_ipr_file = "./entry.list" 	# Load IPR reference list to a dataframe
ipr_df = pd.read_csv(input_ipr_file, sep='\t', skipinitialspace=True)
ipr_df = ipr_df.set_index('ENTRY_AC')

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

def printf(log):
	current_time = datetime.datetime.now()
	formatted_time = current_time.strftime('%Y-%m-%d %H:%M:%S.%f')
	print("{}: {}".format(formatted_time, log))

''' Load scaffold position to a dataframe '''
def load_scafpos(input_position_file):
	col_names = ["scaffold", "position"]
	pos_df =  pd.read_csv(input_position_file, sep=':', comment="#", names=col_names, header=None, skipinitialspace=True)
	pos_df["scaffold"] = pos_df["scaffold"].astype(int)
	pos_df["position"] = pos_df["position"].astype(int)
	pos_df.sort_values(['scaffold', 'position'], axis=0, ascending=True, inplace=True, na_position="first")
	return pos_df

''' Load annotation table and sort the scaffolds and positions '''
def load_annot(input_annotations_file):
	col_names = ["seq_id","source","type","start","end","score","strand","phase","attributes"]
	gene_df = pd.read_csv(input_annotations_file, sep='\t', comment="#", names=col_names, header=None, skipinitialspace=True)
	gene_df["scaffold"] = gene_df["seq_id"].str.split('_')
	gene_df["scaffold"] = gene_df["scaffold"].str[1]
	gene_df["scaffold"] = gene_df["scaffold"].astype(int)
	gene_df.sort_values(['scaffold', 'start', 'end'], axis=0, ascending=True, inplace=True, na_position="first")
	return gene_df




''' Fetch the IPR and GO terms from the GOA file '''
def fetch_ipr_go_info(ann_info):
	#Regex for IPR and GO terms
	ipr_match = re.findall(r'IPR[0-9]+', ann_info) #19,20 (concat IPR using ;)
	ipr_match_str = ';'.join(map(str,ipr_match))
	go_match =  re.findall(r'GO:[0-9]+', ann_info) #21,22 (concat GO using ;)
	go_match_str = ';'.join(map(str,go_match))
	#print "IPR: ", ipr_match_str
	#print "GO : ", go_match_str
	go_term = []
	for g in go_match:
		go_term.append(goont.get(g))
	go_term_str = ';'.join(map(str,go_term))  # Go term for corresponding id
	# Get IPR text
	ipr_term = []
	for ipr in ipr_match:
		try:
			iprv = ipr_df.loc[ipr]
			ipr_term.append(iprv["ENTRY_TYPE"]+"-"+ iprv["ENTRY_NAME"])
		except:
			continue
	ipr_term_str = ';'.join(map(str,ipr_term))  # IPR term for corresponding id
	return ipr_match_str, ipr_term_str, go_match_str, go_term_str





''' Check if the overlap falls in ON or NEAR category based on position for the same scaffold '''
def check_overlap_position(snp_df, ann_item_df):
	pos_snp = snp_df.position
	ann_start, ann_end = ann_item_df.start, ann_item_df.end
	# ON	: position value overlap on (a.) within ann window
	if (pos_snp >= ann_start and pos_snp <= ann_end):							# (a.)
		return "ON"
	# NEAR	: position (b.) upstream, c.) downstream within the acceptable annotation window=10K 
	elif ((pos_snp < ann_start and abs(pos_snp - ann_start) <= snp_window) or 	# (b.)
			(pos_snp > ann_end and abs(pos_snp - ann_end) <= snp_window)	    # (c.)
		):
		return "NEAR"
	# None
	return None


''' Compute cuff position window overlap in the annotation GFF data'''
def compute_position_overlap_annotation(snp_df, gff_df, debug_mode=False):
	out_row = {
				"scaffold"		: None,
				"position"		: None,
				"annstart"		: None,
				"annstop"		: None,
				"ongene"		: 0,
				"oncds"			: 0,
				"onmRNA"		: 0,
				"onexon"		: 0,
				"onte"			: 0,
				"onprotein"		: 0,
				"onmatch"		: 0,
				"neargene"		: 0,
				"nearcds"		: 0,
				"nearmRNA"		: 0,
				"nearexon"		: 0,
				"nearte"		: 0,
				"nearprotein"	: 0,
				"nearmatch"		: 0,
				"numneargenes"	: 0,
				"IPRnumber"		: None,
				"IPRtext"		: None,
				"GOnumber"		: None,
				"GOtext"		: None
			}
	out_row["scaffold"]  = snp_df.scaffold
	out_row["position"]  = snp_df.position
	
	#printf("Processing Snips -  ScaffId:{}, Pos:{}".format(cuff_df["scaffold"], cuff_df["position"]))

	for ann in gff_df.itertuples():	# Pandas Dataframe iterrows is very slow and needs to be changed
		idx = ann.Index
		_type = ann.type
		if _type not in ['gene', 'CDS', 'mRNA', 'exon', 
							'five_prime_UTR', 'three_prime_UTR', 'protein_match', 
							'expressed_sequence_match', 'match', 'match_part']:
			continue # Skip these types

		# Compare if cuff_window is (1.) ON or (2.) NEAR annotations for the same scaffold
		pos_scaff, ann_scaff = snp_df.scaffold, ann.scaffold
		overlap_state = None
		if pos_scaff == ann_scaff:
			overlap_state = check_overlap_position(snp_df, ann)

		if overlap_state is not None and overlap_state == "ON":
			if _type == "gene":
				out_row["ongene"] = 1											   					# ongene
				# Use the start/stop of ON gene annotation	
				out_row["annstart"] = str(ann.start) 												# annstart
				out_row["annstop"]  = str(ann.end)  												# annstop
				# Keep the IPR/GO information for near gene	
				ipr_match_str, ipr_term_str, go_match_str, go_term_str = fetch_ipr_go_info(ann.attributes)	
				out_row["IPRnumber"] = ipr_match_str  												# ipr
				out_row["IPRtext"]   = ipr_term_str													# ipr text - Read this from entry.list file 
				out_row["GOnumber"]  = go_match_str													# go numb
				out_row["GOtext"]    = go_term_str													# go text
			elif _type == "CDS":
				out_row["oncds"] = 1																# oncsd
			elif _type == "mRNA":
				out_row["onmRNA"] = 1																# onmRNA
			elif _type == "exon":
				out_row["onexon"] = 1																# onexon
			elif _type == "five_prime_UTR" or _type == "three_prime_UTR":
				out_row["onte"] = 1																	# onte
			elif _type == "protein_match" or _type == "expressed_sequence_match":
				out_row["onprotein"] = 1															# onprotein
			elif _type == "match" or _type == "match_part":
				out_row["onmatch"] = 1																# onmatch

			# In case gene type did not assign start stop. Use other types too
			if out_row["annstart"] is None: out_row["annstart"] = str(ann.start)				# annstar
			if out_row["annstop"] is None:  out_row["annstop"]  = str(ann.end) 					# annstop				
		elif overlap_state is not None and overlap_state == "NEAR":
			# Go here only if type matches and is not ON for that type
			if _type == "gene" and out_row["ongene"] == 0:
				out_row["neargene"] = 1											   					# neargene
				out_row["numneargenes"] = out_row["numneargenes"] + 1			   					# numneargenes
				# Keep the IPR/GO information for near gene
				ipr_match_str, ipr_term_str, go_match_str, go_term_str = fetch_ipr_go_info(ann.attributes)
				if out_row["IPRnumber"] is not None: 												# Merge logic
						out_row["annstart"]  = out_row["annstart"] + sep_merge + str(ann.start)  	# annstart
						out_row["annstop"]   = out_row["annstop"] + sep_merge + str(ann.end) 		# annstop
						out_row["IPRnumber"] = out_row["IPRnumber"] + sep_merge + ipr_match_str 	# ipr num
						out_row["IPRtext"]   = out_row["IPRtext"] + sep_merge + ipr_term_str		# ipr text 
						out_row["GOnumber"]  = out_row["GOnumber"] + sep_merge + go_match_str		# go numb
						out_row["GOtext"]    = out_row["GOtext"] + sep_merge + go_term_str			# go text						
				else:																				# Persist the current 
					if out_row["annstart"] is None: out_row["annstart"] = str(ann.start)			# annstar
					if out_row["annstop"] is None:  out_row["annstop"]  = str(ann.end) 				# annstop
					out_row["IPRnumber"] = ipr_match_str  											# ipr num
					out_row["IPRtext"]   = ipr_term_str												# ipr text
					out_row["GOnumber"]  = go_match_str												# go numb
					out_row["GOtext"]    = go_term_str												# go text
			elif _type == "CDS" and out_row["oncds"] == 0:
				out_row["nearcds"] = 1																# nearcds
			elif _type == "mRNA" and out_row["onmRNA"] == 0:
				out_row["nearmRNA"] = 1																# nearmRNA
			elif _type == "exon" and out_row["onexon"] == 0:
				out_row["nearexon"] = 1																# nearexon
			elif (_type == "five_prime_UTR" or _type == "three_prime_UTR") and (out_row["onte"] == 0):
				out_row["nearte"] = 1																# nearte
			elif (_type == "protein_match" or _type == "expressed_sequence_match") and (out_row["onprotein"] == 0):
				out_row["nearprotein"] = 1															# nearprotein
			elif (_type == "match" or _type == "match_part") and (out_row["onmatch"] == 0):
				out_row["nearmatch"] = 1															# nearmatch
			# In case near gene type did not assign start stop. Use other types for assignment
			if out_row["annstart"] is None: out_row["annstart"] = str(ann.start)				# annstar
			if out_row["annstop"] is None:  out_row["annstop"]  = str(ann.end) 					# annstop	

		if debug_mode and idx%10000 == 0:
			printf("\t{} Annotation Item: {}".format(idx, snp_df.position))
			printf("\tMatching position overlap:{}".format(out_row))
			
	return out_row

''' Search annotations for cuff window position '''
def search_snpposition_in_annotation(pos_df, gff_df):
	result_list = []
	# Operate on each cuff/transcript id 
	with tqdm(total = pos_df.shape[0]) as pbar: 
		for snpid in pos_df.itertuples():
			pbar.set_description("SnpPos({}), ScaffId({})".format(snpid.position, snpid.scaffold))
			pbar.update(1)
			out_row_dict = compute_position_overlap_annotation(snpid, gff_df, debug_mode=False)
			result_list.append(pd.DataFrame([out_row_dict]))
	# Accumulate and Report result on each cuffid
	result_df = pd.concat(result_list)
	return result_df




if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--pos", help="Add transcript scaffold positions file full path", default="mappos_sub.txt")
	parser.add_argument("--ann", help="Add genome annotations file full path", default="annot1631_sub.txt")
	parser.add_argument("--out", help="Add transcript annotation output file full path", default="sub_snp_annotations_table_1631.out")
	args = parser.parse_args()
	if args.pos:  input_position_file = args.pos
	if args.ann:  input_annotations_file = args.ann
	if args.out:  output_file_name = args.out

	printf("Configs: \n\t POS: {} \n\t ANN: {} \n\t OUT:  {} \n".format(input_position_file, input_annotations_file, output_file_name))
	# Load the input data
	pos_df = load_scafpos(input_position_file)
	#sample_tr_df = pos_df[pos_df["scaffold"] == 1631]
	printf("Loaded scaffold position data: \n {}".format(pos_df.head(4)))
	gff_df = load_annot(input_annotations_file)
	printf("Loaded GFF data: \n {}".format(gff_df.head(4)))

	# compute the positions
	result_df = search_snpposition_in_annotation(pos_df, gff_df)
	printf("Result data: \n {}".format(result_df.head()))

	# Wriite results to file
	result_df.to_csv(output_file_name, sep='\t', header=True, index=False)
	printf("Result writtent to file:  {}".format(output_file_name))