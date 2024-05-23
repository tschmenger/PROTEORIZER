#!/usr/bin/env/python2.7
#import matplotlib
#matplotlib.use('Agg') # or "TkAgg"
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from pprint import pprint as pp
import numpy as np
import gzip
import os
import sys
import subprocess
import re
nums = re.compile(r"[+-]?\d+(?:\.\d+)?")
whitespace_killer=re.compile(r"\s+")
#import seaborn as sns
#import scipy.cluster.hierarchy as sch
#import scipy.spatial.distance as ssd
from Bio import SwissProt as sp
from Bio import SeqIO
from Bio import AlignIO
import math
from scipy.stats import binom_test	### sth = binom_test(x,n,p) with x = number of successes, n = number of trials and p = probability and alternative = "two-sided"
import urllib
import logging ### use this in a try: exception: setup, write: logging.exception("message") under except:
import statistics
import time
from Bio import PDB
import svgwrite
import psycopg2
from multiprocessing import Process


start_time = time.time()

START_DIRECT = os.getcwd()
orig_stdout = sys.stdout
gapletters = [".","-"]
###########################################################	#### #  #     #  #   # #### #### ####	###########################################################
###########################################################     #    #  #    ##  ##  # #    #    #   	###########################################################
###########################################################     #    ####   #### # # # # ## ###   #	###########################################################
###########################################################     #    #  #  #   # #  ## #  # #      #	###########################################################
###########################################################     #### #  # #    # #   # #### #### ####	###########################################################
### 2024 04 25
# 1) Updated COSMIC and gnomAD data & corresponding code in section "STEP 4"
# 2) refined the verdict categories, adding "Ambiguous" as an outcome and adding medium and high confidence labels to the "No_Impact" category
### 2024 05 08
# 3) Introduced a shorthand command for the RHOA Y34C testcase

#################################################################################################################################################################################
###	### CALL THIS SCRIPT FROM THE "PDB_STRUCTURES" directory like this 'python ../CONNECTOR_alpha_Jan2023.py test_unip test_mutat test_gn project erledigt_one erledigt_two erledigt_three nololli
if "testcaseRHOA" in sys.argv:
	test_unip	=	"P61586"		### Uniprot ID	
	test_mutat_raw	=	"[Y34C,E40K,M134T]"	### Mutations, in python list format
	test_gn 	= 	"RHOA"			### Gene Name
	project		=	"Singular_Queries"	### Project Name
	erledigt_one	=	"FALSE"			### TRUE or FALSE, default FALSE
	erledigt_two	=	"FALSE"			### defaults to FALSE, can be used to supply a custom sequence alignment
	erledigt_three	=	"FALSE"			### TRUE or FALSE, default FALSE
	nololli		=	"FALSE"			### FALSE or TRUE, default is FALSE
	clustermethod	=	"RW"			### RW or HClust, default is RW
	databank	=	"Uniprot"		### Uniprot, Humsavar or Both, default should be Uniprot
else:
	test_unip	=	sys.argv[1]	### Uniprot ID	
	test_mutat_raw	=	sys.argv[2]	### Mutations, in python list format
	test_gn 	= 	sys.argv[3]	### Gene Name
	project		=	sys.argv[4]	### Project Name
	erledigt_one	=	sys.argv[5]	### TRUE or FALSE, default FALSE
	erledigt_two	=	sys.argv[6]	### defaults to FALSE, can be used to supply a custom sequence alignment
	erledigt_three	=	sys.argv[7]	### TRUE or FALSE, default FALSE
	nololli		=	sys.argv[8]	### FALSE or TRUE, default is FALSE
	clustermethod	=	sys.argv[9]	### RW or HClust, default is RW
	databank	=	sys.argv[10]	### Uniprot, Humsavar or Both, default should be Uniprot

try:
	r_counter = sys.argv[11]
except:
	r_counter = "none"

test_mutat 	= 	test_mutat_raw.strip('][').split(',')
for item in test_mutat:
	item = item.replace("'","")

test_length = 0

### python ToL_Orthofinder_AlphaFold_Any_Clusterer_batch_unlim_Mechismo.py P23946 [H66R] CMA1 Orthofound_ToL FALSE FALSE FALSE

#### making a protein specific directory and then supply the required files
os.chdir("/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/")
if erledigt_one == "FALSE":
	command = "mkdir "+project+"/"+test_unip
	os.system(command)
os.chdir("/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+test_unip)
########################################################### 	###	#####	#####	#####		#####	###########################################################
###########################################################	#	  #	#	#   #		#   #	###########################################################
###########################################################	 #	  #	####	#####		#   #	###########################################################
###########################################################	  #       #	#	#		#   #	###########################################################
###########################################################	###       #	#####	#		#####	###########################################################
start_time_0 = time.time()
### Here I would like to fetch Mechismo Results for our protein of interest.
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#mechismo = {}	# a dict which will hold the mechismo information
#fpaths=[]	# will store the location of mechismo files
def MECHISMO_GETTER(uniprot_identifier):
	mechismo_dict = {}
	filepaths = {}
	with open("/net/home.isilon/ag-russell/bq_tschmenger/PhD/NaturalVariants/1000Genomes/FILTERING/All_Paths_Cleaned.txt","r") as paths:
		for line in paths:
			if uniprot_identifier in line:
				if filepaths.has_key(uniprot_identifier)==False:
					filepaths[uniprot_identifier]=line
	try:	
		mechismofile = filepaths[uniprot_identifier].replace("\n","")
		with gzip.open(mechismofile,"r") as mech:
			for line in mech:
				if "abling" in line:
					protein = line.split("\t")[0]
					uniprot = line.split("\t")[1].replace(" ","")
					unimuta = line.split("\t")[6].replace(" ","")
					mut = unimuta.split("/")[1]
					muta = mut[-1:]
					partner = line.split("\t")[18]
					conf = line.split("\t")[26]
					if conf == "high":
						con = "h"
					elif conf == "medium":
						con = "m"
					else:
						con = "l"
					score = line.split("\t")[27]	### determination possible on whether its positive or negative
					posi = line.split("\t")[3].replace(" ","")
					pos = int(posi)
					string = partner+"/"+con
					if mechismo_dict.has_key(uniprot) == False:
						mechismo_dict[uniprot]={}
						mechismo_dict[uniprot][pos]={}
						mechismo_dict[uniprot][pos][string]=[float(score)]
					elif mechismo_dict[uniprot].has_key(pos)==False:	
						mechismo_dict[uniprot][pos]={}
						mechismo_dict[uniprot][pos][string]=[float(score)]
					elif mechismo_dict[uniprot][pos].has_key(string)==False:
						mechismo_dict[uniprot][pos][string]=[float(score)]
					else:
						mechismo_dict[uniprot][pos][string].append(float(score))
		#print pp(mechismo)
		for unip in mechismo_dict:
			for po in mechismo_dict[unip]:
				stringler_to_keep = ""
				for part in mechismo_dict[unip][po]:
					avg_score = sum(mechismo_dict[unip][po][part])/len(mechismo_dict[unip][po][part])
					stringler_to_keep = stringler_to_keep + part+"/"+str(avg_score)+"; "
				mechismo_dict[unip][po]=stringler_to_keep
		#print pp(mechismo)
	except:
		pass
	return mechismo_dict
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def fetch_BROT(inputa):
	uniprot = inputa
	try:
		content=urllib.urlopen("https://www.uniprot.org/uniprot/" + uniprot + ".txt")
		for line in content:
		### 		DR   PDB; 1A2B; X-ray; 2.40 A; A=1-181.
		###		DR   PDBsum; 1A2B; -.
			if "ID" and "_HUMAN" in line:		### ID   RHOA_HUMAN              Reviewed;         193 AA.
				lengelenge = int(line.split(";")[1].split("AA")[0].replace(" ","").replace("\n",""))
	except:
		logging.exception("message")
	return lengelenge
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mechismo = MECHISMO_GETTER(test_unip)
test_length  =  fetch_BROT(test_unip)

end_time_0 = time.time()
differ_0 = end_time_0 - start_time_0
print "Step 0 time: ",differ_0
########################################################### 	###	#####	#####	#####		  ##	###########################################################
###########################################################	#	  #	#	#   #		 # #	###########################################################
###########################################################	 #	  #	####	#####		#  #	###########################################################
###########################################################	  #       #	#	#		   #	###########################################################
###########################################################	###       #	#####	#		   #	###########################################################
start_time_1 = time.time()

### first, find the Orthofinder ID
### it is then required to get the uniprot IDs from only those model organisms as the tree of life is very large and has a lot of
### sequences that we do not need.
# OG0000967	P23946
### next, find the correct Orthofinder alignment file and save the aligned sequences of species I like in a dictionary
### uniprot ID stored in: 	test_unip

# alignments={}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def HOMOL_GETTER(uniprot_identifier, referencemutations):

	model_raw = os.listdir("/net/home.isilon/ds-russell/orthologs/tree_of_life/orthologs")
	#	ANOGA_orthologs.tsv.gz	
	#	CAEEL_orthologs.tsv.gz
	model_organisms = []
	for each in model_raw:
		organis = each.split("_")[0]
		if organis not in model_organisms:
			model_organisms.append(organis)
	with gzip.open("/net/home.isilon/ds-russell/orthologs/tree_of_life/orthologs/HUMAN_orthologs.tsv.gz","r") as eggid:
		checker = 0
		saeulen = []
		for line in eggid:
			if checker == 0:
				if "Orthogroup" in line:
					checker = 1
					index = 0
					for item in line.split("\t"):
						if item.replace(" ","") in model_organisms:
							saeulen.append(index)
						index += 1
#					print saeulen
			else:
				unip = line.split("\t")[1].replace("\n","").replace(" ","")	### since I am looking at the human file, human protein IDs are in column 1
				eggnog_id = line.split("\t")[0].replace("\n","").replace(" ","")
				if uniprot_identifier in unip:
					poi_eggnog_id = eggnog_id
					columns = line.split("\t")
					ids_of_interest = []
					for i in range(1, len(columns)):
#						print i
						if i in saeulen:
							subcolumn = columns[i]
							for item in subcolumn.split(","):
								if len(item) > 1:
									ids_of_interest.append(item.replace(" ",""))
					break
	
	#print ids_of_interest
	
	#align_dir = "/net/home.isilon/ds-russell/orthologs/tree_of_life/alignments"		#### old path
	align_dir = "/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/"+test_unip[0:4]
	
	alignments_to_keep = {}

	###			OG0000967.fa.gz	
	#alignments_of_interest = poi_eggnog_id+".fa.gz"	#### old way of finding the file
	alignments_of_interest = uniprot_identifier+".aln.gz"


	#print targetpath
	targetpath = align_dir+"/"+alignments_of_interest
	with gzip.open(targetpath,"r") as alignfile:
		#for record in SeqIO.parse(alignfile,"fasta"):	#### old way of handling the file
		for record in AlignIO.read(alignfile,"clustal"):
			seq_name = record.id.split("_")[1]
			if seq_name in ids_of_interest:	
			       	seq = record.seq
				taxa = seq_name.split(".")[0]
				if alignments_to_keep.has_key(seq_name)==False:
					alignments_to_keep[seq_name]=str(seq)
			if seq_name == test_unip:
				truemutstring = ""
				charcounter = 1
				testdictionary = {}
				for mutation in referencemutations.split(","):
					mutation = mutation.replace("[","").replace("]","")
					refaa = mutation[0]
					refpossi = nums.search(str(mutation)).group(0)
					if testdictionary.has_key(int(refpossi))==False:
						testdictionary[int(refpossi)]=[[refaa],[mutation[1:]]]
					else:
						testdictionary[int(refpossi)][1].append(mutation[1:])
#				print testdictionary
				for character in seq:
#					print character, charcounter
					if character != "-":	
						if charcounter in testdictionary:
#							print charcounter, character
							if character in testdictionary[charcounter][0]:
								pass
							else:
								for entry in testdictionary[charcounter][1]:
									if len(referencemutations.split(","))>1:
										truemutstring = character+entry+","
									else:
										truemutstring = character+entry
						charcounter+=1				
	return alignments_to_keep, truemutstring	
	### Step 1 should take less than 0.2 secs
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def FASTA_MAKER(uniprot_identifier, project_name, dictionary):
	outseqfile= "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project_name+"/"+uniprot_identifier+"/UsedSequences_unlim.fasta"
	fastafile_out = open(outseqfile,"w")
	sys.stdout = fastafile_out
	for k in alignments:
		print ">",k
		print alignments[k]
	sys.stdout = orig_stdout
	fastafile_out.close()
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def CUSTOM_ALIGN(targetfile, projectname, uniprot_identifier):
	alignments_to_keep = {}
	with open(targetfile,"r") as alignfile:
		for record in AlignIO.read(alignfile,"clustal"):
			if "." in record.id:	### UniprotID.Number
				seq_name = record.id.split(".")[0]
			elif "_" in record.id:	### HUMAN_UniprotID
				seq_name = record.id.split("_")[1]
			elif "|" in record.id:	### 00|P46734|MAP2K3
				seq_name = record.id.split("|")[1]
			else:
				seq_name = record.id
			seq = record.seq
			if alignments_to_keep.has_key(seq_name)==False:
				alignments_to_keep[seq_name]=str(seq)
	return alignments_to_keep	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if erledigt_two == "FALSE":
	alignments, truevariants = HOMOL_GETTER(test_unip, test_mutat_raw)
	FASTA_MAKER(test_unip, project, alignments)
else:
#	oldfile = erledigt_two
#	newcommand = "cp "+erledigt_two+" /net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+test_unip+"/UsedSequences_unlim.fasta"
#	os.system(newcommand)
	alignments = CUSTOM_ALIGN(erledigt_two, project, test_unip)
	FASTA_MAKER(test_unip, project, alignments)
	truevariants = ""

end_time_1 = time.time()
differ_1 = end_time_1 - start_time_1
print "Step 1 time: ",differ_1
########################################################### 	###	#####	#####	#####		#####	###########################################################
###########################################################	#	  #	#	#   #		    #	###########################################################
###########################################################	 #	  #	####	#####		 ####	###########################################################
###########################################################	  #       #	#	#		#	###########################################################
###########################################################	###       #	#####	#		#####	###########################################################
start_time_2 = time.time()
###	###	### Handle the Alphafold PDB-style structure
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def ALPHA_Handler(uniprot_identifier, project_name):

	pdb_counter = 0
	correct_wd = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/"
	#alpha_dir = "/net/home.isilon/ds-russell/alphafold/HUMAN/AF-"
	alpha_dir = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/DBs/Alphafold/AF-"
	alpha_folded = "AF-"+uniprot_identifier+"-F1-model_v3.pdb"


	#erledigt = "FALSE"	### I can change this to true if the loop has finished once, saving the intermediate results
	if erledigt_three == "FALSE":	
		targetdirectory = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project_name+"/"+uniprot_identifier
		os.chdir(targetdirectory)
		### AlphaFold file names: 	AF-P61586-F1-model_v1.pdb.gz
		copycommand = "cp "+alpha_dir+uniprot_identifier+"-F1-model_v3.pdb.gz "+targetdirectory
		os.system(copycommand)
		for filename in os.listdir(targetdirectory):
    			if filename.endswith(".pdb.gz"):
				pdb_counter += 1
				new_fn = filename.split(".gz")[0]
				command_in = "gzip -d /net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project_name+"/"+test_unip+"/"+filename	
				os.system(command_in)
			
			elif filename.endswith(".pdb"):
				pdb_counter += 1
			else:
				continue
		os.chdir(correct_wd)
	else:
		targetdirectory = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project_name+"/"+uniprot_identifier
		os.chdir(targetdirectory)
		for filename in os.listdir(targetdirectory):
			if filename.endswith(".pdb"):
				pdb_counter += 1
	return pdb_counter
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pdb_counter = ALPHA_Handler(test_unip, project)

end_time_2 = time.time()
differ_2 = end_time_2 - start_time_2
print "Step 2 time: ",differ_2
########################################################### 	###	#####	#####	#####		#####	###########################################################
###########################################################	#	  #	#	#   #		    #	###########################################################
###########################################################	 #	  #	####	#####		 ####	###########################################################
###########################################################	  #       #	#	#		    #	###########################################################
###########################################################	###       #	#####	#		#####	###########################################################
start_time_3 = time.time()
###	###	### Look for variant counts for this particular protein
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
end_time_3 = time.time()
differ_3 = end_time_3 - start_time_3
print "Step 3 time: ",differ_3
########################################################### 	###	#####	#####	#####		#	###########################################################
###########################################################	#	  #	#	#   #		#   #	###########################################################
###########################################################	 #	  #	####	#####		#####	###########################################################
###########################################################	  #       #	#	#		    #	###########################################################
###########################################################	###       #	#####	#		    #	###########################################################
start_time_4 = time.time()
### 	###	Now I will go ahead and check uniprot for any functional info of positions in those homol. proteins
### 	###	I will store this info in the dictionary homol_counts["UniProt"]
translatedictionary = {}
count_dict = {"Hereditary":{},	### Mutation: [Disease]
		"ClinVar":{},	### Mutation: [Disease]
		"COSMIC":{},	### Mutation: Count
		"UniProt":{},	### Mutation: [Observation]
		"Orthos":{}	### Mutation: HomolProtein/HomolProteinPosition/FunctionInfo
			}
homol_counts = {"Hereditary":{},	### Mutation: [Disease]
		"ClinVar":{},		### Mutation: [Disease]
		"COSMIC":{},		### Mutation: Count
		"UniProt":{}		### Mutation: [Observation]
}
extra_dictionary = {"Hereditary":{},	### Mutation: [Disease]
		"COSMIC":{},		### Mutation: Count
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def HOMOL_UNIP_GETTER(inputdictionary_one, inputdictionary_two, inputdictionary_three, uniprot_identifier, inputdictionary_four, canon_gene):
	conn = psycopg2.connect(database="proteorizer", 
            host="pevolution2.bioquant.uni-heidelberg.de",
            user="bq_tschmenger",
	    password="12345678910")
	cursor = conn.cursor()
	for protein in inputdictionary_one:
		if protein == uniprot_identifier:
			############
			try:
				commander_farsight = "SELECT namus FROM idnames WHERE(unip=\'"+protein+"');"
				cursor.execute(commander_farsight)
				results = cursor.fetchone()
				longname = results[0]
				if translatedictionary.has_key(protein)==False:
					translatedictionary[protein]=longname
			except:
				pass
			
			############
			try:
				commander_farsight = "SELECT namus FROM genenames WHERE(unip=\'"+protein+"');"
				cursor.execute(commander_farsight)
				results = cursor.fetchone()
				langname = results[0]
				translatedictionary[protein]=langname
			except:
				pass
			############
			try:
				commander_farsight = "SELECT variant, ac, het, hom, freq FROM naturalvariants_two WHERE(unip=\'"+protein+"');"
				cursor.execute(commander_farsight)
				for result in cursor.fetchall():
					variant, allelecount, hetcount, homcount, frequency = result
					kommentar = str(allelecount)+"-"+str(hetcount)+"-"+str(homcount)+"-"+str(frequency)
					if inputdictionary_four["Hereditary"].has_key(variant)==False:
						inputdictionary_four["Hereditary"][variant]=[]
						inputdictionary_four["Hereditary"][variant]=[kommentar.replace(" ","")]
					else:
						inputdictionary_four["Hereditary"][variant].append(kommentar.replace(" ",""))
			except:
				pass
			############
			try:
				commander_farsight = "SELECT count, variant FROM cosmic_two WHERE(unip=\'"+canon_gene+"');"
				cursor.execute(commander_farsight)
				for result in cursor.fetchall():
					count, variant = result
					if inputdictionary_four["COSMIC"].has_key(variant)==False:
						inputdictionary_four["COSMIC"][variant]=[]
						inputdictionary_four["COSMIC"][variant]=[int(count)]
					else:
						inputdictionary_four["COSMIC"][variant].append(int(count))

			except:
				logging.exception("message")
				pass
			############
			if databank == "Uniprot":
				commander_farsight = " SELECT typus, infotext, position FROM unipinformation WHERE(unip=\'"+protein+"\');"
			elif databank == "Humsavar":
				commander_farsight = " SELECT annotation, disease, mutation FROM humsavar WHERE(unip=\'"+protein+"\');"
			else:
				commander_farsight = [" SELECT typus, infotext, position FROM unipinformation WHERE(unip=\'"+protein+"\');",
						      " SELECT annotation, disease, mutation FROM humsavar WHERE(unip=\'"+protein+"\');"]
			if isinstance(commander_farsight,list):
				for query in commander_farsight:
					cursor.execute(query)
					for result in cursor.fetchall():
						checktype, comment, homol_mutat = result
						if databank in ["Humsavar","Both"]:
								if checktype != "LB/B " and checktype != "US ":
									kommentar = checktype+"/"+comment.replace("/","=")
									if inputdictionary_three["UniProt"].has_key(homol_mutat)==False:
										inputdictionary_three["UniProt"][homol_mutat]=[]
										inputdictionary_three["UniProt"][homol_mutat]=[kommentar]
									else:
										inputdictionary_three["UniProt"][homol_mutat].append(kommentar)
						else:
								kommentar = checktype+"/"+comment.replace("/","=")
								if inputdictionary_three["UniProt"].has_key(homol_mutat)==False:
									inputdictionary_three["UniProt"][homol_mutat]=[]
									inputdictionary_three["UniProt"][homol_mutat]=[kommentar]
								else:
									inputdictionary_three["UniProt"][homol_mutat].append(kommentar)
			else:
				cursor.execute(commander_farsight)
				for result in cursor.fetchall():
					checktype, comment, homol_mutat = result
					kommentar = checktype+"/"+comment.replace("/","=")
					if databank in ["Humsavar","Both"]:
						if checktype != "LB/B " and checktype != "US ":
							kommentar = checktype+"/"+comment.replace("/","=")
							if inputdictionary_three["UniProt"].has_key(homol_mutat)==False:
								inputdictionary_three["UniProt"][homol_mutat]=[]
								inputdictionary_three["UniProt"][homol_mutat]=[kommentar]
							else:
								inputdictionary_three["UniProt"][homol_mutat].append(kommentar)
					elif databank == "Uniprot":
						kommentar = checktype+"/"+comment.replace("/","=")
						if inputdictionary_three["UniProt"].has_key(homol_mutat)==False:
							inputdictionary_three["UniProt"][homol_mutat]=[]
							inputdictionary_three["UniProt"][homol_mutat]=[kommentar]
						else:
							inputdictionary_three["UniProt"][homol_mutat].append(kommentar)
		else:
			############
			try:
				commander_farsight = "SELECT namus FROM idnames WHERE(unip=\'"+protein+"');"
				cursor.execute(commander_farsight)
				results = cursor.fetchone()
				longname = results[0]
				if translatedictionary.has_key(protein)==False:
					translatedictionary[protein]=longname
			except:
				pass
			
			############
			try:
				commander_farsight = "SELECT namus FROM genenames WHERE(unip=\'"+protein+"');"
				cursor.execute(commander_farsight)
				results = cursor.fetchone()
				langname = results[0]
				translatedictionary[protein]=langname
			except:
				pass
			
			############
			if databank == "Uniprot":
				commander_farsight = " SELECT typus, infotext, position FROM unipinformation WHERE(unip=\'"+protein+"\');"
			elif databank == "Humsavar":
				commander_farsight = " SELECT annotation, disease, mutation FROM humsavar WHERE(unip=\'"+protein+"\');"
			else:
				commander_farsight = [" SELECT typus, infotext, position FROM unipinformation WHERE(unip=\'"+protein+"\');",
						      " SELECT annotation, disease, mutation FROM humsavar WHERE(unip=\'"+protein+"\');"]
			if isinstance(commander_farsight,list):
				for query in commander_farsight:
					cursor.execute(query)
					for result in cursor.fetchall():
						checktype, comment, homol_mutat = result
						if databank in ["Humsavar","Both"]:
							if checktype != "LB/B " and checktype != "US ":
								kommentar = checktype+"/"+comment.replace("/","=")
								if inputdictionary_two["UniProt"].has_key(protein)==False:
									inputdictionary_two["UniProt"][protein]={}
									inputdictionary_two["UniProt"][protein][homol_mutat]=kommentar
								else:
									inputdictionary_two["UniProt"][protein][homol_mutat]=kommentar
						else:
							kommentar = checktype+"/"+comment.replace("/","=")
							if inputdictionary_two["UniProt"].has_key(protein)==False:
								inputdictionary_two["UniProt"][protein]={}
								inputdictionary_two["UniProt"][protein][homol_mutat]=kommentar
							else:
								inputdictionary_two["UniProt"][protein][homol_mutat]=kommentar	
			else:
				cursor.execute(commander_farsight)
				for result in cursor.fetchall():
						checktype, comment, homol_mutat = result
						if databank in ["Humsavar","Both"]:
							if checktype != "LB/B " and checktype != "US ":
								kommentar = checktype+"/"+comment.replace("/","=")
								if inputdictionary_two["UniProt"].has_key(protein)==False:
									inputdictionary_two["UniProt"][protein]={}
									inputdictionary_two["UniProt"][protein][homol_mutat]=kommentar
								else:
									inputdictionary_two["UniProt"][protein][homol_mutat]=kommentar
						else:
							kommentar = checktype+"/"+comment.replace("/","=")
							if inputdictionary_two["UniProt"].has_key(protein)==False:
								inputdictionary_two["UniProt"][protein]={}
								inputdictionary_two["UniProt"][protein][homol_mutat]=kommentar
							else:
								inputdictionary_two["UniProt"][protein][homol_mutat]=kommentar			
	return inputdictionary_two, inputdictionary_three, inputdictionary_four
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
homol_counts, count_dict, extra_dictionary = HOMOL_UNIP_GETTER(alignments, homol_counts, count_dict, test_unip, extra_dictionary, test_gn)
#print pp(extra_dictionary)
#print pp(translatedictionary)
end_time_4 = time.time()
differ_4 = end_time_4 - start_time_4
print "Step 4 time: ",differ_4
########################################################### 	###	#####	#####	#####		#####	###########################################################
###########################################################	#	  #	#	#   #		#   	###########################################################
###########################################################	 #	  #	####	#####		#####	###########################################################
###########################################################	  #       #	#	#		    #	###########################################################
###########################################################	###       #	#####	#		#####	###########################################################
start_time_5 = time.time()
### We now map back functional positions in homologs to the protein of interest
### it is possible, but currently not done, to filter out positions based on sequence conservation thresholds
### currently we keep everything and maintain the ability to filter out residues with low conservation later, if needed
### Note, that keeping low-conserved residues of interest in the data pool also affects the community finding process
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def BACKMAPPER(inputdictionary, uniprot_identifier, einflussdictionary):
	storage = []
	idler = []
	length = []
	original = uniprot_identifier	
	for protein in  inputdictionary:
		sequence = inputdictionary[protein]
		storage.append(sequence)		### storing all sequences 
		idler.append(protein)			### storing all uniprot IDs
		if "." in protein:
			protein = protein.split(".")[0]
		if protein == uniprot_identifier:
			truesequence = sequence
		length.append(len(sequence))		### storing all sequence lengths 
	maximum = max(length)				### determining the max length (sequence + gaps)
	chk_dict = {}	
	for i in range(1,len(truesequence)+1):			### creating a dictionary for each possible position, AA or gap
		amino = truesequence[i-1]
		if chk_dict.has_key(i) == False:
			chk_dict[i]={}
			chk_dict[i][amino]=[0,[]]
		elif chk_dict[i].has_key(amino)==False:
			chk_dict[i][amino]=[0,[]]	
		else:
			pass
	cnt = 0
	for i in idler:	
		if "." in i:
			i = i.split(".")[0]				### determining which index the queried uniprot ID has to extract the correct sequence
		if str(i) == str(test_unip):
			resp_number = cnt
		else:
			cnt += 1
	ref_sequ = storage[resp_number]			### saving the original and therefore reference sequence, including gaps
	divider = float(len(storage))			### saving the number of sequences we can test to later calculate the % of conservation
	correctseqnumber = 0
	for protein in homol_counts["UniProt"]:
		cnt = 0
		for i in idler:					### determining which index the queried uniprot ID has to extract the correct sequence
			if str(i) == str(uniprot_identifier):
				protein_id = cnt
			else:
				cnt += 1
		for possi in homol_counts["UniProt"][protein]:	### then, I need to go through each of the functionally interesting positions in that protein
			functional 	= homol_counts["UniProt"][protein][possi]
			position 	= nums.search(str(possi)).group(0)	### task is to find the corresponding position in the alignment
			orig		= possi
			truechecker = 1
			wrongchecker = 1
			for character in inputdictionary[protein]:
				if character != "-":
	                		truechecker += 1
	                		wrongchecker += 1
	        		else:
	                		wrongchecker += 1
	       			if int(truechecker)==int(position):
	                		interresult = wrongchecker	### this gives me the position within the alignment for the mutation protein/X123Y
			try:
				ref_sequ_pos = ref_sequ[interresult-1]		### this gives me the reference position from the queried protein
										### Since I start the counter at 1, I might see Methionine at a count of 1, but as
										### the first entry in a list I need to use 1-1 = 0
			except:
				logging.exception("message")
				pass   
				
			stringler = protein+"/"+str(possi)+"/"+functional+";"
			try:
		        	if str(ref_sequ_pos) == str(orig):
		               		chk_dict[int(interresult)][ref_sequ_pos][0]+=1.0
					chk_dict[int(interresult)][ref_sequ_pos][1].append(stringler)
				else:
					chk_dict[int(interresult)][ref_sequ_pos][1].append(stringler)
			except:
				logging.exception("message")
				pass    
	zaehler = 1
	for k in range(1, len(chk_dict)):
		for a in chk_dict[k]:
			if "-" != str(a):
				amino = ref_sequ[k-1]
				stringler = ""
				for item in chk_dict[k][a][1]:
					stringler = stringler + item + " " 
				if float(len(chk_dict[k][a][1])) != 0.0:
					if (float(chk_dict[k][a][0])/float(len(chk_dict[k][a][1]))) >= 0.00:				### NOT GOING FOR A SPECIFIC PERCENTAGE OF SEQUENCE IDENTITY ANYMORE
						if einflussdictionary["Orthos"].has_key(zaehler)==False:
							einflussdictionary["Orthos"][zaehler]=[stringler]
						else:
							einflussdictionary["Orthos"][zaehler].append(stringler)
				zaehler+=1
	return einflussdictionary, zaehler
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
count_dict, counter = BACKMAPPER(alignments, test_unip, count_dict)

end_time_5 = time.time()
differ_5 = end_time_5 - start_time_5
print "Step 5 time: ",differ_5
########################################################### 	###	#####	#####	#####		#####	###########################################################
###########################################################	#	  #	#	#   #		#   	###########################################################
###########################################################	 #	  #	####	#####		#####	###########################################################
###########################################################	  #       #	#	#		#   #	###########################################################
###########################################################	###       #	#####	#		#####	###########################################################
start_time_6 = time.time()
###	###	Prepare the colleced data to be submitted to cluster analysis and then calculate the distances
###	###	I need a dictionary like this:			#fiftypercent_ug={"Q6P1J9":[70,88,100,110]}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def GATHERING(uniprot_identifier,projectname,inputdictionary,inputmutation, counter):
	gatherer = []	### I will simply use this to briefly see how the ratio between original residues and homol residues looks like
	workingdictionary={}
	workingdictionary[uniprot_identifier]=[]
	outfile = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+uniprot_identifier+"/PositionsOfInterest_unlim.txt"
	result_file = open(outfile,"w")
	sys.stdout = result_file
	new_resis = 0
	old_resis = 0
	for param in inputdictionary:
		for v in inputdictionary[param]:			# position
			if param != "Orthos":
				if v not in gatherer:
					gatherer.append(v)
			else:
				new_resis = len(inputdictionary[param])
			old_resis = len(gatherer)
			if v not in workingdictionary[uniprot_identifier]:
				try:
					va = nums.search(str(v)).group(0)
					if va not in workingdictionary[uniprot_identifier]:
						workingdictionary[uniprot_identifier].append(int(va))
				except:
					pass
#	['P06780/E45I/MUTAGEN/_Temperature_sensitive_growth_defect."_; P06780/E45V/MUTAGEN/_In_RHO1-2;_temperature_sensitive,_fails_to_/activate/PKC1."/; ']
##	P06780/E45I/MUTAGEN/_Temperature_sensitive_growth_defect."_|P06780/E45V/MUTAGEN/_In_RHO1-2;_temperature_sensitive,_fails_to_activate/PKC1."/|
#	['P06780/T22N/MUTAGEN/_Abolishes_GTP-binding."_; ']
##	P06780/T22N/MUTAGEN/_Abolishes_GTP-binding."_|
			infostring = str(inputdictionary[param][v]).replace("[","").replace("]","").replace("'","").replace("\n","").replace("_/","_").replace("; ","|").replace("\"/","")
#	P06780/E45I/MUTAGEN/Temperaturesensitivegrowthdefect.;P06780/E45V/MUTAGEN/InRHO1-2;temperaturesensitive,failstoactivate/PKC1."/; 
##	P06780/E45I/MUTAGEN/Temperaturesensitivegrowthdefect.|P06780/E45V/MUTAGEN/InRHO1-2;temperaturesensitive,failstoactivate/PKC1."/|
#	P06780/T22N/MUTAGEN/AbolishesGTP-binding.;
##	P06780/T22N/MUTAGEN/AbolishesGTP-binding.| 
			print param,"\t", v,"\t",infostring
	sys.stdout = orig_stdout
	result_file.close()
	for entries in inputmutation:
		test_pos = nums.search(entries).group(0)
		if test_pos not in workingdictionary[uniprot_identifier]:
			workingdictionary[uniprot_identifier].append(int(test_pos))
	print "Number of functionally interesting residues in queried protein",uniprot_identifier,"is",old_resis
	print "Number of functionally interesting residues added from homologous proteins is",new_resis
	############# CHECKPOINT 
	### do sys.exit() if we cannot find any PDB structures
	if counter <= 0:
		print "Not enough PDB structures"
		sys.exit()
	return workingdictionary
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
positions_of_interest = GATHERING(test_unip, project, count_dict, test_mutat, pdb_counter)

end_time_6 = time.time()
differ_6 = end_time_6 - start_time_6
print "Step 6 time: ",differ_6
########################################################### 	###	#####	#####	#####		#####	###########################################################
###########################################################	#	  #	#	#   #		    #  	###########################################################
###########################################################	 #	  #	####	#####		  ###	###########################################################
###########################################################	  #       #	#	#		  #  	###########################################################
###########################################################	###       #	#####	#		 #   	###########################################################
start_time_7 = time.time()
###	###	Prepare the colleced data to be submitted to cluster analysis and then calculate the distances
###	###	I need a dictionary like this:			#fiftypercent_ug={"Q6P1J9":[70,88,100,110]}	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def PDB_READY(uniprot_identifier,lenge):
	notification_message = "####Note: Distances calculated using only AlphaFold prediction"
	alpha_folded = "AF-"+uniprot_identifier+"-F1-model_v3.pdb"
	outdictionary = {}
	if outdictionary.has_key(uniprot_identifier)==False:
		notification = 0	
		outdictionary[uniprot_identifier]={}
		outdictionary[uniprot_identifier][alpha_folded]={}
		outdictionary[uniprot_identifier][alpha_folded]["A"]={}
		outdictionary[uniprot_identifier][alpha_folded]["A"]["PDB"]=[1,lenge]
		outdictionary[uniprot_identifier][alpha_folded]["A"]["Uniprot"]=[1,lenge]
	elif outdictionary[uniprot_identifier].has_key(alpha_folded)==False:
		notification = 1
		outdictionary[uniprot_identifier][alpha_folded]={}
		outdictionary[uniprot_identifier][alpha_folded]["A"]={}
		outdictionary[uniprot_identifier][alpha_folded]["A"]["PDB"]=[1,lenge]
		outdictionary[uniprot_identifier][alpha_folded]["A"]["Uniprot"]=[1,lenge]
	else:
		pass
	return outdictionary, notification, notification_message
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sifts, notification, notification_message = PDB_READY(test_unip, test_length)

end_time_7 = time.time()
differ_7 = end_time_7 - start_time_7
print "Step 7 time: ",differ_7
########################################################### 	###	#####	#####	#####		#####	###########################################################
###########################################################	#	  #	#	#   #		#   #  	###########################################################
###########################################################	 #	  #	####	#####		#####	###########################################################
###########################################################	  #       #	#	#		#   # 	###########################################################
###########################################################	###       #	#####	#		#####  	###########################################################
start_time_8 = time.time()
###	###	Now lets calculate the distances!
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def DISTANCE_MAKER(PDB_dict, uniprot_identifier, projectname, positiondictionary):
	parser = PDB.PDBParser()
	for protein in PDB_dict:
		distances = {}
		for i in range(1,counter+1):
			if i in positiondictionary[protein]:
				if distances.has_key(i)==False:
					distances[i]={}
		for k in distances:
			for i in range(1,counter+1):
				if i in positiondictionary[protein]:
					if distances[k].has_key(i)==False:
						distances[k][i]=[]	### basically creating a n*n matrix where n are the MUTATED positions of 'protein'
	
		for structure in PDB_dict[protein]:
			PDB_ID = ""
			if "AF-" in structure:
				PDB_ID = structure
			else:	### this only makes sense for non-AlphaFold PDB codes
				for i in range(0,4):
						testnumberletter = structure[i]
						try:
							numberletter = testnumberletter.upper()		#### this is required, since PDB codes in SIFTS are stored in lowercase but PDB files are in uppercase
							PDB_ID = PDB_ID+numberletter			#### I think doing structure.upper() would be enough, actually. No need for looping.
						except:
							#logging.exception("message")
							PDB_ID = PDB_ID+testnumberletter
	
			for chainler in PDB_dict[protein][structure]:
				#print "STRUCTURE","\t",structure										### DEBUGGING
				#print PDB_dict[protein][structure][chainler]									### DEBUGGING
				if "AF-" in structure:
					pdb1 = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+uniprot_identifier+"/"+str(PDB_ID)
				else:
					pdb1 = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+uniprot_identifier+"/"+str(PDB_ID)+".pdb"
				pdb_f = PDB_dict[protein][structure][chainler]["PDB"][0]
				pdb_l = PDB_dict[protein][structure][chainler]["PDB"][1]
				uni_f = PDB_dict[protein][structure][chainler]["Uniprot"][0]
				uni_l = PDB_dict[protein][structure][chainler]["Uniprot"][1]
				### first I need to check if the PDB structure covers mutated residues
				possible = []
				for a in range(int(uni_f), int(uni_l)+1):
					if a in positiondictionary[protein]:
						possible.append(a)

				if len(possible) >= 2:	### only in this case it would make sense to continue
					structurus = parser.get_structure(PDB_ID,pdb1) 
					model = structurus[0] 
					chain = model[chainler] 
					for index in range (int(pdb_f), int(pdb_l)+1):
						for b in possible:
							if b != index:
								if index in possible:
									try:	### some residues in the middle of the structure might be missing, jeopardizing the analysis
										#print index, "\t", b						### DEBUGGING			
										residue1 = chain[index] 
										residue2 = chain[int(b)]
	
										if "GLY" in str(residue1):
											atom1 = residue1['CA']		### the C-beta atom is the start of the side chain and therefore more interesting
											if "GLY" in residue2:		### to look at
												atom2 = residue2['CA']	### Glycin does not have a C-beta atom, though
												distance = atom1-atom2 
												if "AF-" in structure:
													bval1 = atom1.get_bfactor()
													bval2 = atom2.get_bfactor()
													if float(bval1) >= 50.0:		### "Confident" Alpha Zero predictions start at 70/100
														if float(bval2) >= 50.0: 	### I only check the CA of the main chain for confidence
															distances[index][b].append(distance)
												else:
													distances[index][b].append(distance)
	
	
											else:
												atom2 = residue2['CB'] 
												atom3 = residue2['CA']
												dist_one = atom1-atom2	### To avoid considering distances going through the proteins backbone
												dist_two = atom1-atom3	### I will compare CB-CB distance vs. CB-CA distance, and use the shorter of the two
												check_dist = [dist_one,dist_two]
												if "AF-" in structure:
													bval1 = atom1.get_bfactor()
													bval3 = atom3.get_bfactor()
													if float(bval1) >= 50.0:	### "Confident" Alpha Zero predictions start at 70/100
														if float(bval3) >= 50.0: ### I only check the CA of the main chain for confidence
															distances[index][b].append(min(check_dist))
												else:
													distances[index][b].append(min(check_dist))
	
	
	
										else:
											atom1 = residue1['CB']
											atom3 = residue1['CA']
											if "GLY" in residue2:
												atom2 = residue2['CA']	### Glycin does not have a C-beta atom
												dist_one = atom1-atom2	### To avoid considering distances going through the proteins backbone
												dist_two = atom3-atom2	### I will compare CB-CB distance vs. CB-CA distance, and use the shorter of the two
												check_dist = [dist_one,dist_two]
												if "AF-" in structure:
													bval3 = atom3.get_bfactor()
													bval2 = atom2.get_bfactor()
													if float(bval3) >= 50.0:	### "Confident" Alpha Zero predictions start at 70/100
														if float(bval2) >= 50.0: ### I only check the CA of the main chain for confidence
															distances[index][b].append(min(check_dist))
												else:
													distances[index][b].append(min(check_dist))
	
	
											else:
												atom2 = residue2['CA']
												atom4 = residue2['CB']
												dist_one = atom1-atom2
												dist_two = atom1-atom4
												check_dist = [dist_one,dist_two]
												if "AF-" in structure:
													bval3 = atom3.get_bfactor()
													bval2 = atom2.get_bfactor()
													if float(bval3) >= 50.0:	### "Confident" Alpha Zero predictions start at 70/100
														if float(bval2) >= 50.0: ### I only check the CA of the main chain for confidence
															distances[index][b].append(min(check_dist))
												else:
													distances[index][b].append(min(check_dist))
										
	
									except:
										pass
										
		#print pp(distances)
		for k in distances:
			for v in distances[k]:
				try:
					distance = sum(distances[k][v])/len(distances[k][v])
					if float(distance) <= 8.0:
						distances[k][v]=distance
					else:
						distances[k][v]=''		
				except:
					#logging.exception("message")
					pass
		df1 = pd.DataFrame.from_dict(distances)
		dataframeoutput ='/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/'+projectname+'/'+protein+'/distances_df_unlim'
		df1.to_csv(dataframeoutput, sep=',', mode='w')

	return dataframeoutput
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_output_path = DISTANCE_MAKER(sifts, test_unip, project, positions_of_interest)

end_time_8 = time.time()
differ_8 = end_time_8 - start_time_8
print "Step 8 time: ",differ_8
########################################################### 	###	#####	#####	#####		#####	###########################################################
###########################################################	#	  #	#	#   #		#   #  	###########################################################
###########################################################	 #	  #	####	#####		#####	###########################################################
###########################################################	  #       #	#	#		    # 	###########################################################
###########################################################	###       #	#####	#		#####  	###########################################################
start_time_9 = time.time()
###	###	Trying to call the 3D clustering R script, giving arguments, and prepare to get a good readout without having to go through tinkering with R-Studio
#filus <- args[1]	should be the input file, here dataframe of distances
#namus <- args[2]	should be the name of the protein
#outfilus <- args[3]	should be the target file for the output
#svgfilus <- args[4]	should be the target file for the svg output
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def RMAKER(uniprot_identifier, projectname, dataframeinput, genename):
	outfile = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+uniprot_identifier+"/"+"ClusterMemberships_unlim.txt"
	svgfilus_path = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+uniprot_identifier+"/"+"ClusterPlot_unlim.svg"
	
	if "RW" in clustermethod:
		rcmd = "Rscript --vanilla "+"/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/Network_Calculations.R "+dataframeinput+" "+genename+" "+outfile+" "+svgfilus_path
	elif "HClust" in clustermethod:
		rcmd = "Rscript --vanilla "+"/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/Network_Calculations_HClust.R "+dataframeinput+" "+genename+" "+outfile+" "+svgfilus_path
	else:	### defaults to RW
		rcmd = "Rscript --vanilla "+"/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/Network_Calculations.R "+dataframeinput+" "+genename+" "+outfile+" "+svgfilus_path
	print rcmd
	os.system(rcmd)
	return outfile
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
outfilus_path = RMAKER(test_unip, project, df_output_path, test_gn)

end_time_9 = time.time()
differ_9 = end_time_9 - start_time_9
print "Step 9 time: ",differ_9
########################################################### 	###	#####	#####	#####       ##	#####	###########################################################
###########################################################	#	  #	#	#   #	   # #	#   #  	###########################################################
###########################################################	 #	  #	####	#####	  #  #	#   #	###########################################################
###########################################################	  #       #	#	#	     #	#   # 	###########################################################
###########################################################	###       #	#####	#	     #	#####  	###########################################################
start_time_10 = time.time()
### ### combine the cluster membership output with the output in InterestingPositions
#ALPP		445	1					### cluster membership output
#Homols 	133 	['P05186/G129R/(in_HOPS)"__/; ']	### Positions of Interest Output
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def KLEBER(inputfile,uniprot_identifier, projectname,mutation, ppidict):
	memberships = {}
	try:
		with open(outfilus_path,"r") as clusterfile:
			for line in clusterfile:
				position = line.split("\t")[1].replace(" ","").replace("\n","")
				cluster = line.split("\t")[2].replace(" ","").replace("\n","")
				if memberships.has_key(int(position))==False:
					memberships[int(position)]=int(cluster)
	except:
		pass
	pymol_colors = ["red","salmon","ruby","forest","marine","density","yellow","sand","magenta","pink","deeppurple","cyan","teal","orange","lightorange"] #list of 15 colors for up to 15 communities/clusters
																		      # default color should be "grey70"
	colorfile =  "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+uniprot_identifier+"/"+"ClusterColors_for_Pymol.txt"
	colorlist = []
	for i in range(1,int(test_length)+1):
		try:
			colornumber = memberships[i]
			color = pymol_colors[int(colornumber)]
			colorlist.append(color)
		except:
			colorlist.append("grey70")
	colfile = open(colorfile,"w")
	sys.stdout = colfile
	print colorlist
	sys.stdout = orig_stdout
	colfile.close()	

	final_file = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+uniprot_identifier+"/FinalResults_unlim.txt"
	res_file = open(final_file,"w")
	sys.stdout = res_file
	print "Clusternumber","\t","Data_Source","\t","Functional_Information","\t","Mechismo_Predictions"
	for entries in mutation:
		try:
			inputposition = memberships[int(nums.search(entries).group(0))]
		except:
			inputposition = "NoClusterMembership"
		print inputposition,"\t", "Input","\t",nums.search(entries).group(0),"\t",entries
	sys.stdout = orig_stdout
	res_file.close()
	### mechismo:					mechismo[unip][po]=stringler_to_keep		&		stringler_to_keep = partner / confidence / average score
	### poifile:					Homols 	160 	['Q5R9Y9/140/MOD_RES/Phosphoserine; ']
	with open("/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+uniprot_identifier+"/PositionsOfInterest_unlim.txt","r") as poifile:
	#print "opened"	
		poifile.seek(0)															### DEBUGGING
		for line in poifile:
			linus = line.replace("\n","")
			#print linus														### DEBUGGING
			checkpos = line.split("\t")[1].replace(" ","").replace("\n","")

			try:
				predictions = ppidict[uniprot_identifier][int(nums.search(checkpos).group(0))]
			except:
				predictions = "No_Mechismo2b_Predictions"

			res_file = open(final_file,"a")
			sys.stdout = res_file
			try:
				clustership = memberships[int(nums.search(checkpos).group(0))]
				print clustership,"\t",linus,"\t",predictions
			except:
				#logging.exception("message")											### DEBUGGING
				print "NoClusterMembership","\t",linus,"\t",predictions
			sys.stdout = orig_stdout
			res_file.close()
	
		if notification == 0:
			res_file = open(final_file,"a")
			sys.stdout = res_file
			print notification_message
			sys.stdout = orig_stdout
			res_file.close()
	return final_file
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
try:
	final_file = KLEBER(outfilus_path, test_unip, project, test_mutat, mechismo)
except:
	pass
end_time_10 = time.time()
differ_10 = end_time_10 - start_time_10
print "Step 10 time: ",differ_10
########################################################### 	###	#####	#####	#####       ##	  ##	###########################################################
###########################################################	#	  #	#	#   #	   # #	 # #    ###########################################################
###########################################################	 #	  #	####	#####	  #  #	#  #	###########################################################
###########################################################	  #       #	#	#	     #	   # 	###########################################################
###########################################################	###       #	#####	#	     #	   #  	###########################################################
start_time_11 = time.time()
### Renaming the directory correctly
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def RENAMER(projectname, uniprot_identifier, genus, erledigt_two):
	if erledigt_two == "FALSE":
		new_dir_name = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+uniprot_identifier+"-"+genus
	else:
		new_dir_name = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+uniprot_identifier+"-"+genus+"-Custom"
	old_dir_name = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+uniprot_identifier
	os.rename(old_dir_name,new_dir_name)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
RENAMER(project, test_unip, test_gn, erledigt_two)

end_time_11 = time.time()
differ_11 = end_time_11 - start_time_11
print "Step 11 time: ",differ_11
########################################################### 	###	#####	#####	#####    ##   #####	###########################################################
###########################################################	#	  #	#	#   #   # #       #     ###########################################################
###########################################################	 #	  #	####	#####  #  #    ####	###########################################################
###########################################################	  #       #	#	#	  #   #    	###########################################################
###########################################################	###       #	#####	#	  #   #####  	###########################################################
start_time_12 = time.time()
### Making Lollipop plots
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def Lollipopper(project, uniprot_identifier, gene_name, COSMIC_data):
	if len(COSMIC_data) <= 0:
		COSMICdata = {}
		with gzip.open("/net/home.isilon/ds-russell/COSMIC/latest/CosmicCompleteTargetedScreensMutantExport.tsv.gz","r") as cosmicfile:
			for line in cosmicfile:
				if "Accession Number" not in line:
					if "Missense" in line:
						genus 		= line.split("\t")[0].split("_")[0]
						variant 	= line.split("\t")[20].split(".")[1]
						if genus == gene_name:
							if COSMICdata.has_key(genus)==False:
								COSMICdata[genus]={}
								COSMICdata[genus][variant]=1
							elif COSMICdata[genus].has_key(variant)==False:
								COSMICdata[genus][variant]=1
							else:
								COSMICdata[genus][variant]+=1
	else:
		pass

	inputs = []
	uniprots  = []
	cosmic  = []
	filus = workingdir+"/FinalResults_unlim.txt"
	with open(filus,"r") as datafile:
		for line in datafile:
			if "Clusternumber" not in line:
				if "#" not in line:
					typus = line.split("\t")[1].replace(" ","")
					position = line.split("\t")[2].replace(" ","")
					if typus == "Input":
						inputs.append(position)
					elif typus == "UniProt":
						uniprots.append(position)
					else:
						pass
	stringler = ""
	for variant in inputs:
		substringler = variant+"#3399ff"+"@5"
		stringler = stringler+" "+substringler.replace("'","").replace("\"","")
	for variant in uniprots:
		substringler = variant+"#00a86b"+"@1"
		stringler = stringler+" "+substringler.replace("'","").replace("\"","")
	try:
		for variant in COSMICdata[gene]:
			checkvalue = COSMICdata[gene][variant]
			if int(checkvalue) >= 15:
				substringler = variant+"#fcf030"+"@1"
				stringler = stringler+" "+substringler.replace("'","").replace("\"","")	
	except:
		pass
	command = "./lollipops -legend -labels -o="+gene_name+".svg -U "+uniprot_identifier+" "+stringler	### see this issue on github for using the -U option: https://github.com/joiningdata/lollipops/issues/59
	print command
	os.chdir("/home/bq_tschmenger/lollipops-v1.6.0-linux64")
	os.system(command)
	command_two = "mv "+gene_name+".svg "+workingdir+"/"+gene_name+".svg" 
	print command_two
	os.system(command_two)
	os.chdir(workingdir)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if erledigt_two == "FALSE":
	workingdir = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+test_unip+"-"+test_gn
else:
	workingdir = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+test_unip+"-"+test_gn+"-Custom"

if "FALSE" in nololli:	
	pass
elif "TRUE" in nololli:
	COSMIC_dataset = {}	### change this IF it is desired, a full dictionary may be submitted once for multiqueries when incorporated into another script
	Lollipopper(project, test_unip, test_gn, COSMIC_dataset)
else:
	pass

os.chdir(workingdir)

end_time_12 = time.time()
differ_12 = end_time_12 - start_time_12
print "Step 12 time: ",differ_12
########################################################### 	###	#####	#####	#####    ##   #####	###########################################################
###########################################################	#	  #	#	#   #   # #       #     ###########################################################
###########################################################	 #	  #	####	#####  #  #    ####	###########################################################
###########################################################	  #       #	#	#	  #       #    	###########################################################
###########################################################	###       #	#####	#	  #   #####  	###########################################################
start_time_13 = time.time()
### calculating sequence conservation of relevant residues
aminogroups = {	"positive":["R","H","K"],
		"negative":["D","E"],
		"polar,uncharged":["S","T","N","Q"],
		"others":["G","P","C"],
		"hydrophobic":["A","V","I","L","M","F","Y","W"]}

orig_stdout = sys.stdout
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def conservation_checker(identifier):
	for protein in  alignments:
		sequence = alignments[protein]
		if protein == identifier:
			truesequenzler = sequence
	sequencelength = len(truesequenzler)
	positionalcounter = 1

	conservational_dictionary = {}
	theforbiddenalignpos = []
	for i in range(0,sequencelength):
		identitycontainer = []
		for ident in alignments:
			if identifier in ident:
				orires = alignments[ident][i]
				identitycontainer.append(alignments[ident][i])
			else:
				identitycontainer.append(alignments[ident][i])
		identitypercentage = float(identitycontainer.count(orires))/float(len(identitycontainer))	### so far this also includes "-" as the original truesequence residue, be cautious
		oritype = "none"

		for typus in aminogroups:
			if orires.upper() in aminogroups[typus]:
				oritype = typus
		aminotypes = []
		for k in identitycontainer:
			if k in gapletters:
				aminotypes.append("none")
			else:
				for typus in aminogroups:
					if k.upper() in aminogroups[typus]:
						aminotypes.append(typus)
		typuspercentage = float(aminotypes.count(oritype))/float(len(aminotypes))
		if orires not in gapletters:
#			print positionalcounter,"\t", identitycontainer, "\t", orires,"\t",identitypercentage,"\t",typuspercentage ,"\t", aminotypes,"\t",oritype
			if conservational_dictionary.has_key(int(positionalcounter))==False:
				conservational_dictionary[int(positionalcounter)] = [float(identitypercentage), orires, float(typuspercentage), oritype]
			positionalcounter+=1
		elif orires in gapletters:
			if identitypercentage >= 0.90:
				theforbiddenalignpos.append(i+1)
		else:
			pass
	return conservational_dictionary, theforbiddenalignpos	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def conservation_maker(uniprot_identifier, projectname, gene_name):
	if erledigt_two == "FALSE":
		filus = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+uniprot_identifier+"-"+gene_name+"/"+"UsedSequences_unlim.fasta"
	else:
		filus = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+uniprot_identifier+"-"+gene_name+"-Custom/"+"UsedSequences_unlim.fasta"
	alignments = {}
	with open(filus,"r") as alignfile:
		for record in SeqIO.parse(alignfile,"fasta"):
			seq_name = record.id
       			seq = record.seq
			if alignments.has_key(seq_name)==False:
				alignments[seq_name]=str(seq)
	### creating the conservation data
	dict_of_interest, forbiddenpositions = conservation_checker(uniprot_identifier)
	#print pp(dict_of_interest)
	try:
		if erledigt_two == "FALSE":
			filus_two = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+uniprot_identifier+"-"+gene_name+"/"+"FinalResults_unlim.txt"
		else:
			filus_two = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+uniprot_identifier+"-"+gene_name+"-Custom/"+"FinalResults_unlim.txt"
		if erledigt_two == "FALSE":
			outfile_name = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+uniprot_identifier+"-"+gene_name+"/FinalResults_unlim_conservation.txt"
		else:
			outfile_name = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+uniprot_identifier+"-"+gene_name+"-Custom/FinalResults_unlim_conservation.txt"
		with open(filus_two,"r") as datafile:
			outfile = open(outfile_name,"w")
			sys.stdout = outfile
			for line in datafile:
				if "Clusternumber" not in line:
					if "#" not in line:
						cluster = line.split("\t")[0].replace(" ","")
						typus = line.split("\t")[1].replace(" ","")
						position = line.split("\t")[2].replace(" ","")
						funcinfo = line.split("\t")[3].replace(" ","").replace("'","").replace(" ","").replace("[","").replace("]","").replace("\n","").replace("="," ").replace("_"," ")
						if "Input" not in typus: 						
							mechismo = line.split("\t")[4].replace(" ","").replace("\n","")
						else:
							mechismo = "-"
						numeric_position = nums.search(position).group(0)
						try:
							identper, originalresidue, typeper, originaltype = dict_of_interest[int(numeric_position)]
							print cluster,"\t",typus,"\t",position,"\t",funcinfo,"\t",mechismo,"\t",int(float(identper)*100),"\t",int(float(typeper)*100)			
						except:
							identper = "-"
							typeper = "-"
							print cluster,"\t",typus,"\t",position,"\t",funcinfo,"\t",mechismo,"\t",identper,"\t",typeper	
				
					
					else:
						print line.replace("\n","")
				else:
 					print "Clusternumber","\t","Data_Source","\t","Position","\t","Functional_Information","\t","Mechismo_Predictions","\t","SeqIdent%","\t","TypeIdent%"
	except:
		#logging.exception("message")
		pass
			
	sys.stdout = orig_stdout
	return dict_of_interest, forbiddenpositions
	outfile.close()
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Konserve, TheForbiddenPositions = conservation_maker(test_unip, project, test_gn)

end_time_13 = time.time()
differ_13 = end_time_13 - start_time_13
print "Step 13 time: ",differ_13
########################################################### 	###	#####	#####	#####    ##   #		###########################################################
###########################################################	#	  #	#	#   #   # #   #   #     ###########################################################
###########################################################	 #	  #	####	#####  #  #   #####	###########################################################
###########################################################	  #       #	#	#	  #       #    	###########################################################
###########################################################	###       #	#####	#	  #   	  #  	###########################################################
start_time_14 = time.time()
### Creating an annotated alignmentfile as .svg
if databank == "Uniprot":
	colordict = {"VARIANT":"lightgreen",
	"MUTAGEN":"salmon",
	"BINDING":"yellow",
	"MOD_RES":"purple",
	"ACT_SITE":"gold"}
elif databank == "Humsavar":
	colordict = {"LB/B":"darkgreen",
	"LP/P":"red",
	"US":"blue"}
else:
	colordict = {"VARIANT":"lightgreen",
	"MUTAGEN":"salmon",
	"BINDING":"yellow",
	"MOD_RES":"purple",
	"ACT_SITE":"gold",
	"LB/B":"darkgreen",
	"LP/P":"red",
	"US":"blue"}

Clustalcolors = {"A":"hydrophobic",
		"I":"hydrophobic",
		"L":"hydrophobic",
		"M":"hydrophobic",
		"F":"hydrophobic",
		"W":"hydrophobic",
		"V":"hydrophobic",
		"C":"hydrophobic",
		"K":"positive",
		"R":"positive",
		"E":"negative",
		"D":"negative",
		"N":"polar",
		"Q":"polar",
		"S":"polar",
		"T":"polar",
		"G":"glycine",
		"P":"proline",
		"H":"aromatic",
		"Y":"aromatic"}
clustaltypes = {"hydrophobic":"blue",			
		"positive":"red",
		"negative":"magenta",
		"polar":"green",
		"glycine":"black",
		"proline":"orange",
		"aromatic":"cyan"}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def interprodownloader(identif):
	import urllib2
	import ast
	BASE_URL = "https://www.ebi.ac.uk/interpro/api/protein/UniProt/"+identif+"/?residues&page_size=200"
	req = urllib2.Request(BASE_URL)
	response = urllib2.urlopen(req)
	the_page = response.read()
	interpro = ast.literal_eval(the_page)	### just for the record, Interpro is the most over-engineered website I ever had to work with, it is basically useless for wet lab scientists
	interpro_processed = {}
	for k in interpro:
		for loca in interpro[k]["locations"]:
			descr = ""
			for b in loca["description"]:	###              'locations': [{'description': 'GEF (guanine nucleotide exchange factor) interaction site',
				descr = descr+b
			if "(" in descr:
				front_descr = descr.split("(")[0]
				back_descr = descr.split(")")[1]
				descr = front_descr+back_descr
			for categ in loca["fragments"]: ### 		  categ are dictionaries, again, because the nesting never ends here
				for element in categ:
					starting = categ["start"]
					ending = categ["end"]
					if int(starting) == int(ending):
						residue = int(starting)
					if interpro_processed.has_key(descr)==False:
						interpro_processed[descr]=[residue]
					else:
						if residue not in interpro_processed[descr]:
							interpro_processed[descr].append(residue)
	return	interpro_processed
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def SHOWORDER(seqs, doi, starti, endi, goi):
	# dictionary of interest, residue of interest, windowsize, gene of interest
	showtime = {}
	for k in doi:	### uniprot ID = k
		featurecount = []
		sequenzler = seqs[k]
		residue = 0
		for i, letter in enumerate(sequenzler,start = 1):
			if letter != "-":
				residue += 1
				if i >= starti:
					if i <= endi:
						for v in doi[k]: ### categories, i.e. VARIANT = v
							for vv in doi[k][v]:	### residue number = vv
								if int(vv) == residue:
									if int(vv) not in featurecount:
										featurecount.append(vv)
		if showtime.has_key(k)==False:
			showtime[k]=len(featurecount)
	ranking = sorted(showtime, key=lambda x: (-showtime[x], x))
	return ranking
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if erledigt_two == "FALSE":
	outfile_name = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+test_unip+"-"+test_gn+"/FinalResults_unlim_conservation.txt"
else:
	outfile_name = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+test_unip+"-"+test_gn+"-Custom/FinalResults_unlim_conservation.txt"
def clustergetter():
	clusterstorage = {}
	clustercolors = ["firebrick","skyblue","orchid","tan","plum","slateblue","peru","crimson"]
	acceptables = []
	with open(outfile_name,"r") as inputfilus:
		for line in inputfilus:
			if "#" not in line:
				if "Data_Source" not in line:
					if "Input" in line:
						
						position_raw	= line.split("\t")[2].replace(" ","")
						position 	= str(nums.search(position_raw).group(0))
						cluster		= line.split("\t")[0].replace(" ","")
						if cluster not in acceptables:
							acceptables.append(cluster)
						if clusterstorage.has_key(position)==False:
							clusterstorage[position]=cluster
					else:
						if "NoClusterMembership" not in line:
							position_raw	= line.split("\t")[2].replace(" ","")
							position 	= str(nums.search(position_raw).group(0))
							cluster		= line.split("\t")[0].replace(" ","")
							if cluster in acceptables:
								if clusterstorage.has_key(position)==False:
									clusterstorage[position]=cluster
	return clusterstorage, clustercolors
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def create_svg(sequences_dict, projectname, colordict, startposition, windowsize, poi, goi, forbidden):
### preparing the input dictionaries we need
### positions =  {"P61586":{"VARIANT":
	annotationsdictionary = {}
	if erledigt_two == "FALSE":
		openfilus = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+poi+"-"+goi+"/PositionsOfInterest_unlim.txt"
	else:
		openfilus = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+projectname+"/"+poi+"-"+goi+"-Custom/PositionsOfInterest_unlim.txt"

	clusters, clusterhighlights = clustergetter()
	clusterfilus = open('clusterfile.txt','w')
    	clusterfilus.write(str(clusters))
    	clusterfilus.close()
	with open(openfilus,"r") as poifile:
	#print "opened"	
		poifile.seek(0)															### DEBUGGING
		for line in poifile:	#UniProt 	Y34A 	['MUTAGEN/_Abolishes_interaction_with_DGKQ."_']
			typus 		= line.split("\t")[0].replace(" ","")
			position_raw	= line.split("\t")[1].replace(" ","")
			position	= int(nums.search(position_raw).group(0))
			komment		= line.split("\t")[2]
			if typus == "UniProt":
				try:
					for item in komment.split("|"):
						for ke in colordict:
							if ke in komment:
								if ke == "BINDING":
									firstbindingpos = int(nums.search(position_raw.split("..")[0]).group(0))
									lastbindingpos = int(nums.search(position_raw.split("..")[1]).group(0))
									for listposition in range(firstbindingpos,lastbindingpos+1):
										if annotationsdictionary.has_key(poi)==False:
											annotationsdictionary[poi]={}
											annotationsdictionary[poi][ke]=[listposition]
										elif annotationsdictionary[poi].has_key(ke)==False:
											annotationsdictionary[poi][ke]=[listposition]
										else:
											annotationsdictionary[poi][ke].append(listposition)
								else:
									if annotationsdictionary.has_key(poi)==False:
										annotationsdictionary[poi]={}
										annotationsdictionary[poi][ke]=[position]
									elif annotationsdictionary[poi].has_key(ke)==False:
										annotationsdictionary[poi][ke]=[position]
									else:
										annotationsdictionary[poi][ke].append(position)
				except:
					for ke in colordict:
						if ke in komment:
							if annotationsdictionary.has_key(poi)==False:
								annotationsdictionary[poi]={}
								annotationsdictionary[poi][ke]=[position]
							elif annotationsdictionary[poi].has_key(ke)==False:
								annotationsdictionary[poi][ke]=[position]
							else:
								annotationsdictionary[poi][ke].append(position)
			else:
### Orthos 	188 	['P63000/K186A/MUTAGEN/_Decreased_palmitoylation_by_the_V.cholerae_/toxin/RtxA."/; Q92930/180/MOD_RES/Phosphoserine; ']
				try:
					for item in komment.split("|"): 
						try:
							identifier 	= item.split("/")[0].replace("'","")
							position_raw 	= item.split("/")[1].replace("'","")
							position	= int(nums.search(position_raw).group(0))
							komment		= item.split("/")[2].replace("'","")
							for ke in colordict:
								if ke in komment:
									if annotationsdictionary.has_key(identifier)==False:
										annotationsdictionary[identifier]={}
										annotationsdictionary[identifier][ke]=[position]
									elif annotationsdictionary[identifier].has_key(ke)==False:
										annotationsdictionary[identifier][ke]=[position]
									else:
										annotationsdictionary[identifier][ke].append(position)
						except:
							pass
				except:
					try:
						identifier 	= item.split("/")[0].replace("'","")
						position_raw 	= item.split("/")[1].replace("'","")
						position	= int(nums.search(position_raw).group(0))
						komment		= item.split("/")[2].replace("'","")
						for ke in colordict:
							if ke in komment:
								if annotationsdictionary.has_key(identifier)==False:
									annotationsdictionary[identifier]={}
									annotationsdictionary[identifier][ke]=[position]
								elif annotationsdictionary[identifier].has_key(ke)==False:
									annotationsdictionary[identifier][ke]=[position]
								else:
									annotationsdictionary[identifier][ke].append(position)
					except:
						pass
	positioncolors = list(colordict.values())
	featurecolors = ["firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink",
			"firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink",
			"firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink"]
	colors = {}
	coloringcategories = []
	for k in annotationsdictionary:
		counter = 0
		for v in annotationsdictionary[k]:
			if colors.has_key(v)==False:
				if v not in coloringcategories:
					coloringcategories.append(v)
				colors[v]=positioncolors[counter]
				counter+=1
    	for item in annotationsdictionary:
       		for categ in colors:
            		if annotationsdictionary[item].has_key(categ)==False:
                		annotationsdictionary[item][categ]=[]
	#print pp(annotationsdictionary)	#<<<<<<<<<<<<<<<<<<<<<<<<<<<
	positionfile = open('positiondictionary.txt','w')
    	positionfile.write(str(annotationsdictionary))
    	positionfile.close()
	translationfile = open('translationfile.txt','w')
    	translationfile.write(str(translatedictionary))
    	translationfile.close()
	feature_dict = {}
	try:
		feature_dict = interprodownloader(poi)
	except:
		logging.exception("message")
		feature_dict = {}
	#print pp(feature_dict)			#<<<<<<<<<<<<<<<<<<<<<<<<<<<
	featurefile = open('featurefile.txt','w')
    	featurefile.write(str(feature_dict))
    	featurefile.close()
	#### prep phase done, now we generate the svg!
	dwg = svgwrite.Drawing('sequence.svg', profile='full')
   	x = 50
    	y = 80
    	sequence_of_interest = sequences_dict[poi]
    	non_minus_count = 0
    	distance_end = len(sequence_of_interest)+100		### to make sure it gets weeded out below, if none of the if statements directly below trigger
   	distance_start = 0					### to make sure it gets weeded out below, if none of the if statements directly below trigger
	multiinputchecker = "FALSE"				### We need this to enable this function to deal with multiple input positions
	if "," in startposition:
		multiinputchecker = "TRUE"
		distance_start  = 1
		distance_end	= len(sequences_dict[poi])
		for i, letter in enumerate(sequence_of_interest,start = 1):
			if letter not in gapletters:
				non_minus_count += 1

	else:
		for i, letter in enumerate(sequence_of_interest,start = 1):
			if letter not in gapletters:
				non_minus_count += 1
				if non_minus_count == int(startposition):
					startpos = i	### this is the alignment position that corresponds to the residue of interest. alignmpent position includes "-"
				if non_minus_count == int(startposition)+windowsize:
					distance_end = i
				if non_minus_count == int(startposition)-windowsize:
					distance_start = i
    	maxcharactercnt = non_minus_count		### should capture the true length of the sequence of interest
	roworder = SHOWORDER(sequences_dict, annotationsdictionary, distance_start, distance_end, poi)
	try:
		roworder = roworder[0:11]
	except:
		pass
    	if distance_start <= 0:
		distance_start = 1
    	if distance_end > len(sequence_of_interest):
		distance_end = len(sequence_of_interest)
	maximumdistance = distance_end - distance_start
	viewboxcounter = 1
    	highlightingID = 0
    	highlightsaver = {}
	for uniprot in roworder:
		seq 	= sequences_dict[uniprot]
		namus 	= uniprot
		if multiinputchecker == "FALSE":
			startingpoint = int(startposition) - windowsize	### this is required for the correct labeling according to the sequence of interest
		else:
			startingpoint = 0

		try:
			newname = translatedictionary[uniprot]
		except:
			newname = uniprot
		if poi in namus:	##### this if/else conditional can probably be put in yet another function to reduce the amount of code being used here
			old_x = x
			old_y = y
			x = 50
			y = 60
			dwg.add(dwg.rect((x-88, y), (80, 14), fill="yellow"))
			if len(newname) < 8:
				dwg.add(dwg.text(newname, insert = (x-55,y+7), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
			else:
				dwg.add(dwg.text(newname, insert = (x-55,y+7), text_anchor='middle', dominant_baseline='central', font_size='7px', font_family='Arial', font_weight='bold', fill='black'))
		else:
			if len(newname) < 8:
				dwg.add(dwg.text(newname, insert = (x-55,y+7), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
			else:
				dwg.add(dwg.text(newname, insert = (x-55,y+7), text_anchor='middle', dominant_baseline='central', font_size='7px', font_family='Arial', font_weight='bold', fill='black'))
		totalcount = 0
		if startingpoint <= 0:
			startnumberlabel = 1
		elif startingpoint >= maxcharactercnt:
			startnumberlabel = maxcharactercnt
		else:
			startnumberlabel = startingpoint 
		charactercount = 0
		tempfeat = {}
		featcount = 0
		firstdone = "false"#
		lastdone = "false"#
		forbidden_start = "false"#
		forbidden_end = "false"#
		gapcounter = 0#
		for i, letter in enumerate(seq, start=1):
		    totalcount += 1		#### gives the alignment position, including gaps
		    letter = seq[i-1]
		    #print letter, "\t", "okay1"
		    if letter not in gapletters:
		    	charactercount += 1
			#print letter, "\t", "okay2", "\t", charactercount,"\t",totalcount,"\t",distance_end
		    	if totalcount <= distance_end:	### distance_end refers to the last alignment position that will be considered, which is +20 non-gap residues from the input position
				#print letter, "\t", "okay3"
				endcounter = charactercount
				testlenge = int(distance_end)-int(totalcount)				
				if testlenge <= maximumdistance:	### checks that we still operate around the position of interest +/- residues only
				    #print letter, "\t", "okay4"
				    if totalcount >= distance_start: ###distance_start refers to the first alignment position that will be considered, which is +20 non-gap residues from the input position
					if totalcount >= distance_start:
						if firstdone == "false":
							forbidden_start = "true"
							startcounter = charactercount
							dwg.add(dwg.text(charactercount, insert=(35, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))					
							firstdone = "true"
					if totalcount not in forbidden:
						if poi in namus:
						    viewboxcounter += 1
						    if float(Konserve[startnumberlabel][0])>= 0.7:
						    	dwg.add(dwg.rect((x,y),(10,len(roworder)*20), fill= clustaltypes[Clustalcolors[letter.upper()]], opacity=0.2))

					    	    for feat in feature_dict:
							if startnumberlabel in feature_dict[feat]:
								if tempfeat.has_key(feat)==False:
									tempfeat[feat]=[featurecolors[featcount],featcount]
									featcount+=1
								elevator = tempfeat[feat][1]
								elevator_floor = 0
								if elevator >= 10:
									if elevator_floor <= 10:
										elevator = elevator_floor
										elevator_floor += 1
									else:
										elevator_floor = 0
										elevator = elevator_floor
								y_level = 0 + (elevator*3)
								y_level_text = -55 + (elevator*4.5)
								dwg.add(dwg.rect((x, y_level), (10, 2), fill=tempfeat[feat][0]))
								if "done" not in tempfeat[feat]:
									dwg.add(dwg.text(str(feat), insert=(x+15, y_level_text), text_anchor='middle', dominant_baseline='central', font_size='5px', font_family='Arial', font_weight='bold', fill=tempfeat[feat][0]))
									tempfeat[feat].append("done")

						    if multiinputchecker == "TRUE":
							inppossis = []
							for usersubmit in startposition.split(","):
								usernubmer = nums.search(usersubmit).group(0)
								if str(startnumberlabel) == str(usernubmer):
									inppossis.append(str(startnumberlabel))
									dwg.add(dwg.rect((x, 40), (10, 14*float(Konserve[startnumberlabel][0])), fill="black"))
        			        				dwg.add(dwg.text(str(startnumberlabel), insert=(x+5, 35), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill='red'))	
					    		else:
								if int(startnumberlabel)>= startingpoint:
									dwg.add(dwg.rect((x, 40), (10, 14*float(Konserve[startnumberlabel][0])), fill="black"))
									if int(startnumberlabel)%10 == False:
										if str(startnumberlabel) not in inppossis:
											dwg.add(dwg.text(str(startnumberlabel), insert=(x+5, 35), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill='black'))
							if clusters.has_key(str(startnumberlabel)):
								if clusters[str(startnumberlabel)] != "NoClusterMembership":
									dwg.add(dwg.text(str(startnumberlabel), insert=(x+5, 35), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill=clusterhighlights[int(clusters[str(startnumberlabel)])-1]))
					    	    else:
					    		if startnumberlabel == int(startposition):
								position_interest_x = x
								position_interest_y = y
        			        			dwg.add(dwg.rect((x, 40), (10, 14*float(Konserve[startnumberlabel][0])), fill="black"))
        			        			dwg.add(dwg.text(str(startnumberlabel), insert=(x+5, 35), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill='red'))	
					    		else:
								if int(startnumberlabel)>= startingpoint:
									dwg.add(dwg.rect((x, 40), (10, 14*float(Konserve[startnumberlabel][0])), fill="black"))
									if int(startnumberlabel)%10 == False:
										dwg.add(dwg.text(str(startnumberlabel), insert=(x+5, 35), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill='black'))
							if clusters.has_key(str(startnumberlabel)):
								if clusters[str(startnumberlabel)] != "NoClusterMembership":
									dwg.add(dwg.text(str(startnumberlabel), insert=(x+5, 35), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill=clusterhighlights[int(clusters[str(startnumberlabel)])-1]))	 
					    	    startnumberlabel+=1
					
						try:
							drawn = 0
							radius = 7
                                                        hightlightstring = ""
							for colorcateg in coloringcategories:
			        		  		if charactercount in annotationsdictionary[namus][colorcateg]:
                                            				if drawn != 1:
                                                				highlightingID += 1
                                                				hightlightstring = newname+"/"+str(charactercount) + "|" + colorcateg
                                            				else:
                                                				hightlightstring = hightlightstring + "/" + colorcateg 
	        				       			#dwg.add(dwg.rect((x, y), (10, 14), fill=colordict[colorcateg]))
									dwg.add(dwg.circle((x+5, y+7.5), (radius), fill=colordict[colorcateg]))
	        				        		dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
									drawn = 1
								radius -= 1
                                    			if drawn == 1:
                                        			if highlightsaver.has_key(str(highlightingID))==False:
                                            				highlightsaver[str(highlightingID)]=[x+5,y+7.5,hightlightstring]
	        					if drawn == 0:
	       					        	dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))			
				    		except:
							dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))	
        			    		x += 10
				    
		    else:	### will draw just a "-" for a gap in the alignment
			if totalcount >= distance_start:
				if totalcount <= distance_end:
					if totalcount not in forbidden:
						#dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='white'))
		        			x += 10
		viewboxcounter = x
		lastx = x
		lasty = y
		finalresidue = startcounter+gapcounter+(2*windowsize)
		dwg.add(dwg.text(endcounter, insert=(lastx+20, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))
		if poi in namus:
			dwg.add(dwg.rect((-38,y),(x+38,14), fill="none",stroke="black",stroke_width=1))	
			x = 50
			y = old_y
		else:
	      	 	x = 50
 	        	y += 20

		
    	viewboxwidth = (viewboxcounter+150)
    	viewboxheight = len(sequences_dict)*20+140
    	dwg.viewbox(-40, -80,viewboxwidth,viewboxheight)
	try:
    		dwg.add(dwg.rect((position_interest_x, position_interest_y), (10, len(roworder)*20),fill="none",stroke="black",stroke_width=1))
	except:
		pass

	x = 50
	y = 0
	for category in colordict:
		 dwg.add(dwg.rect((x-30, y-70), (60, 10), fill=colordict[category]))
	    	 dwg.add(dwg.text(category, insert=(x, y-65), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
	 	 x += 60
#	dwg.save()
#        newannoalignname = "AnnotatedAlignment_"+str(windowsize)+".svg"
#	os.rename("sequence.svg",newannoalignname)
#
#    	styletext = """<style>
#   		<![CDATA[
#    		text.moo {
#         		font-family: "arial";
#         		fill: black;
#         		font-size: 100%;
#    			}
#    		rect.hiss {
#         		fill:white;
#    			}
#   			]]>
#   		svg text.moo {display: none;}
#   		svg rect.hiss {display: none;}
#   		svg g:hover text {display: block;}
#   		svg g:hover rect {display: block;}
#		</style>"""
#
#    	imagefile = open(newannoalignname,"r")
#    	data= imagefile.read()
#    	data = data.replace("</svg>", styletext+"</svg>")
#	imagefile.close()
#    	writeFile = open(newannoalignname, "w")
#    	writeFile.write(data)
#	writeFile.close()
#
#    	circletext = ""
#    	for hlid in highlightsaver:
#        	cx = highlightsaver[hlid][0]
#        	cy = highlightsaver[hlid][1]
#        	txt = highlightsaver[hlid][2]
#        	uppertext = txt.split("|")[0]
#        	lowertext = txt.split("|")[1]
#        	circletext = circletext+"""<g xmlns="http://www.w3.org/2000/svg">
#          <circle xmlns="http://www.w3.org/2000/svg" cx='"""+str(cx)+"""' cy='"""+str(cy)+"""' r="7" style="fill:transparent;stroke:transparent;stroke-width:0.5;fill-opacity:0.25;stroke-opacity:0.25"/>      
#          <rect class="hiss" x='"""+str(cx)+"""' y='"""+str(cy-40)+"""' height='40' width='"""+str(len(lowertext)*10+10)+"""'></rect>
#          <text class="moo" x='"""+str(cx)+"""' y='"""+str(cy-28)+"""'><tspan class="text">"""+uppertext+"""</tspan> <tspan class="text" x='"""+str(cx)+"""' dy='20' >"""+lowertext+"""</tspan></text>
#          </g>"""
#
#    	imagefile = open(newannoalignname,"r")
#    	imagefile.seek(0)
#
#    	data = imagefile.read()
#	imagefile.close()
#    	data = data.replace("</svg>", circletext+"</svg>")
#
#    	writeFile = open(newannoalignname, "w")
#    	writeFile.write(data)
#	writeFile.close()
        
		

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
try:
#   def create_svg(sequences_dict, projectname, colordict, startposition, windowsize, poi, goi):
#### drawing the complete alignment may consume a lot of time if the proteins have a decent size and if there are a bunch of them in the alignment
	if "," in test_mutat_raw:
		create_svg(alignments, project, colordict, test_mutat_raw, 5, test_unip, test_gn, TheForbiddenPositions)
		annotationcommand = "/home/bq_tschmenger/anaconda2/bin/python /net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/STABLE_Annotate_Alignment_V9_Proteorizer.py "+test_unip+" "+test_mutat_raw+" 30 "+ "UsedSequences_unlim.fasta positiondictionary.txt featurefile.txt 15 clusterfile.txt translationfile.txt"
		os.system(annotationcommand)
	else:
		input_position = nums.search(str(test_mutat)).group(0)
		create_svg(alignments, project, colordict, input_position, 5, test_unip, test_gn, TheForbiddenPositions)
		annotationcommand = "/home/bq_tschmenger/anaconda2/bin/python /net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/STABLE_Annotate_Alignment_V9_Proteorizer.py "+test_unip+" "+input_position+" 30 "+ "UsedSequences_unlim.fasta positiondictionary.txt featurefile.txt 15 clusterfile.txt translationfile.txt"
		os.system(annotationcommand)


except:
	logging.exception("message")
	pass


end_time_14 = time.time()
differ_14 = end_time_14 - start_time_14
print "Step 14 time: ",differ_14
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########################################################### 	###	#####	#####	#####    ##   #####	###########################################################
###########################################################	#	  #	#	#   #   # #   #         ###########################################################
###########################################################	 #	  #	####	#####  #  #   #####	###########################################################
###########################################################	  #       #	#	#	  #       #    	###########################################################
###########################################################	###       #	#####	#	  #   #####  	###########################################################
start_time_15 = time.time()
# Mutation Scoring
all_types = np.array([['','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'],
                ['A',	   0,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2],
                ['C',	   0,  0, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
                ['D',	  -2, -3,  0,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3],
                ['E',	  -1, -4,  2,  0, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2],
                ['F',     -2, -2, -3, -3,  0, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3],
                ['G',	   0, -3, -1, -2, -3,  0, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3],
                ['H',	  -2, -3, -1,  0, -1, -2,  0, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2],
                ['I',	  -1, -1, -3, -3,  0, -4, -3,  0, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1],
                ['K',	  -1, -3, -1,  1, -3, -2, -1, -3,  0, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2],	
                ['L',	  -1, -1, -4, -3,  0, -4, -3,  2, -2,  0,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1],
                ['M',	  -1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  0, -2, -2,  0, -1, -1, -1,  1, -1, -1],
                ['N',	  -2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  0, -2,  0,  0,  1,  0, -3, -4, -2],
                ['P',	  -1, -2, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  0, -1, -2, -1, -1, -2, -4, -3],
                ['Q',	  -1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  0,  1,  0, -1, -2, -2, -1],
                ['R',	  -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  0, -1, -1, -3, -3, -2],
                ['S',	   1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  0,  1, -2, -3, -2],
                ['T',	   0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  0,  0, -2, -2],
                ['V',	   0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  0, -3, -1],
                ['W',	  -3, -2, -4, -3,  1, -3, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3,  0,  2],
                ['Y',	  -2, -2, -3, -2,  3, -2,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  0]])
allt = pd.DataFrame(data=all_types[1:,1:],index=all_types[1:,0],columns=all_types[0,1:])

def AAScorer(muta):
	try:
		first = muta[:1]
		second = muta[-1:]
		aachange = allt.loc[first, second]
	except:
		aachange = "?"
	return aachange



### Here we introduce bayesian & machine learning scoring
def bayesian_scoring(seqval, typval, protclustaval, clustmembval,evidtotalval,evidclustaval, mutation_value):
        seq_score       = 0
        typus_score     = 0
        protclusta_score= 0
        clustmemb_score = 0
        evid_total_score= 0
        evid_clusta_score = 0
	mutation_score = 0
	bayesian = {}
### https://www.nature.com/articles/s41525-022-00322-z#Sec2
### bayes score is sum of log2(TPR/FPR)
	if "RW" in clustermethod:
		metricfilepath = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/CONNECTOR_Scoring/Version5/3_RandomWalk_Metrics_20230818.txt"
	elif "HClust" in clustermethod:
		metricfilepath = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/CONNECTOR_Scoring/Version_HClust3/3_Metrics_20230821_HClust3fixed.txt"
	else:
		metricfilepath = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/CONNECTOR_Scoring/Version5/3_RandomWalk_Metrics_20230818.txt"

	with open(metricfilepath, "r") as metricfile:
    		for line in metricfile:
		##Category 	Threshold 	TP 	FP 	TN 	FN 	TPR 	FPR
		##SeqIdentity 	1.0 	19391 	1814 	104 	553 	0.967518211755 	0.945776850886
        		if "Threshold" not in line:
 		        	categ   =   line.split("\t")[0].replace(" ","").replace("\n","")
		          	grenze  =   float(line.split("\t")[1].replace(" ","").replace("\n",""))
  		          	tpr     =   float(line.split("\t")[6].replace(" ","").replace("\n",""))
  		          	fpr     =   float(line.split("\t")[7].replace(" ","").replace("\n",""))
  		          	if bayesian.has_key(categ)==False:
  		              		bayesian[categ]={}
  		              		bayesian[categ][grenze]=[tpr, fpr]
  		         	elif bayesian[categ].has_key(grenze)==False:
  		          		bayesian[categ][grenze]=[tpr, fpr]
   		          	else:
       			  		pass

	collected = {"TotalEvidence":[],"Clustermembers":[],"ProteinClusters":[],"SeqIdentity":[], "TypeIdentity":[],"ClusterEvidence":[], "Mutscore":[]}

	for cate in bayesian:
    		for thresh in bayesian[cate]:
        		collected[cate].append(thresh)

	for cate in collected:
    		sorts = sorted(collected[cate])
    		collected[cate] = sorts
	try:
		try:
                    for i in collected["SeqIdentity"]:
                        if float(seqval)>= i:
                            trueposi  = bayesian["SeqIdentity"][i][0]
                            falseposi = bayesian["SeqIdentity"][i][1]
                    seq_score = math.log((trueposi/falseposi),2)
                except:
                    pass
                try:
                    for i in collected["TypeIdentity"]:
                        if float(typval)>= i:
                            trueposi  = bayesian["TypeIdentity"][i][0]
                            falseposi = bayesian["TypeIdentity"][i][1]
                    typus_score = math.log((trueposi/falseposi),2)
                except:
                    pass
                try:
                    for i in collected["ProteinClusters"]:
                        if float(protclustaval)>= i:
                            trueposi  = bayesian["ProteinClusters"][i][0]
                            falseposi = bayesian["ProteinClusters"][i][1]
                    protclusta_score = math.log((trueposi/falseposi),2)
                except:
                    pass
		try:                                            
                    for i in collected["Clustermembers"]:
                        if "NoClusterMembership" in line:
                            clustmemb_score = 0.0
                        else:
                            if float(clustmembval)>= i:
                                trueposi  = bayesian["Clustermembers"][i][0]
                                falseposi = bayesian["Clustermembers"][i][1]
                            clustmemb_score = math.log((trueposi/falseposi),2)
                except:
                    pass
                try:
                    for i in collected["TotalEvidence"]:
                        if float(evidtotalval)>= i:
                            trueposi  = bayesian["TotalEvidence"][i][0]
                            falseposi = bayesian["TotalEvidence"][i][1]
                    evid_total_score = math.log((trueposi/falseposi),2)
                except:
                    pass
                try:
                    for i in collected["ClusterEvidence"]:
                        if "NoClusterMembership" in line:
                            evid_clusta_score = 0.0
                        else:                            
                            if float(evidclustaval)>= i:
                                trueposi  = bayesian["ClusterEvidence"][i][0]
                                falseposi = bayesian["ClusterEvidence"][i][1]
                            evid_clusta_score = math.log((trueposi/falseposi),2)
                except:
                    pass
                try:
                    for i in collected["Mutscore"]:
                        if "NoClusterMembership" in line:
                            mutation_score = 0.0
                        else:                            
                            if float(mutation_value)>= i:
                                trueposi  = bayesian["Mutscore"][i][0]
                                falseposi = bayesian["Mutscore"][i][1]
                            mutation_score = math.log((trueposi/falseposi),2)
                except:
                    pass



		bayesianscore = float(seq_score+typus_score+protclusta_score+clustmemb_score+evid_total_score+evid_clusta_score+mutation_score)
                return bayesianscore
	except:
		bayesianscore = 0.0
                return bayesianscore
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
vectorholder = {"Cluster":{},
	"No_Cluster":{} }
def vectormaker(file_of_interest):
	filus = test_unip+"-"+test_gn
	dataholder_input_clust = []
	dataholder_input_noclust = []
	with open(file_of_interest, "r") as foi:
			for line in foi:
		#	Clusternumber 	Data_Source 	Functional_Information 	Mechismo_Predictions
		#	NoClusterMembership 	Input 	167 	F167Y
				if "#" not in line:
					if "Data_Source" not in line:
						cluster		= line.split("\t")[0].replace(" ","")
						typus		= line.split("\t")[1].replace(" ","")
						position	= line.split("\t")[2].replace(" ","")
						entry           = line.split("\t")[3]
						if typus == "UniProt":
                                                    evidencecount = 1
						elif typus == "Input":
						    evidencecount = 0
                                                else:                                        
                                                    evidencecount = entry.count("|")-1
						seq_id          = line.split("\t")[5].replace(" ","")
						if "-" in seq_id:
                                                    seq_id = "0"
						type_id         = line.split("\t")[6].replace(" ","").replace("\n","")
						if "-" in type_id:
                                                        type_id = "0"
						info		= position+"-"+typus+"-"+seq_id+"-"+type_id+"-"+str(evidencecount)+"-"+entry
						if typus == "Input":
							if cluster != "NoClusterMembership":
								inputcluster = cluster

								if vectorholder["Cluster"].has_key(filus) == False:
									vectorholder["Cluster"][filus]={}
									vectorholder["Cluster"][filus][cluster]=[info]
								elif vectorholder["Cluster"][filus].has_key(cluster)==False:
									vectorholder["Cluster"][filus][cluster]=[info]
								else:
                                                                        if info not in vectorholder["Cluster"][filus][cluster]: 
                                                                            vectorholder["Cluster"][filus][cluster].append(info)
							else:
								if position not in dataholder_input_noclust:
									dataholder_input_noclust.append(position)
								if vectorholder["No_Cluster"].has_key(filus) == False:
									vectorholder["No_Cluster"][filus]={}
									vectorholder["No_Cluster"][filus][cluster]=[info]
								elif vectorholder["No_Cluster"][filus].has_key(cluster)==False:
									vectorholder["No_Cluster"][filus][cluster]=[info]
								else:
                                                                        if info not in vectorholder["No_Cluster"][filus][cluster]: 
                                                                            vectorholder["No_Cluster"][filus][cluster].append(info)
						else:
							if cluster != "NoClusterMembership":
								if vectorholder["Cluster"].has_key(filus) == False:
									vectorholder["Cluster"][filus]={}
									vectorholder["Cluster"][filus][cluster]=[info]
								elif vectorholder["Cluster"][filus].has_key(cluster)==False:
									vectorholder["Cluster"][filus][cluster]=[info]
								else:	
                                                                    if info not in vectorholder["Cluster"][filus][cluster]: 
                                                                            vectorholder["Cluster"][filus][cluster].append(info)
							else:
								if vectorholder["No_Cluster"].has_key(filus) == False:
									vectorholder["No_Cluster"][filus]={}
									vectorholder["No_Cluster"][filus][cluster]=[info]
								elif vectorholder["No_Cluster"][filus].has_key(cluster)==False:
									vectorholder["No_Cluster"][filus][cluster]=[info]
								else:	
                                                                    if info not in vectorholder["No_Cluster"][filus][cluster]: 
                                                                            vectorholder["No_Cluster"][filus][cluster].append(info)
	returnal = []
	#print vectorholder
	for k in vectorholder:
		for v in vectorholder[k]:
        		totalcount = 0
        		for clust in vectorholder[k][v]:           
            			for entry in vectorholder[k][v][clust]:
                			evid = entry.split("-")[4]
                			try:
                    				totalcount += int(evid)
                			except:
                    				pass
        		for clusta in vectorholder[k][v]:
                		counter = 0
                		for entr in vectorholder[k][v][clusta]:
					#print entr
                        		evidence = entr.split("-")[4]
                        		counter += int(evidence)
                		for ent in vectorholder[k][v][clusta]: 
					if "Input" in ent:
                        			outputvector = "Output-"+v+"-"+clusta+"-"+ent+"-"+str(len(vectorholder[k][v]))+"-"+str(len(vectorholder[k][v][clusta]))+"-"+str(totalcount)+"-"+str(counter)
						returnal.append(outputvector)
	return returnal
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###
asdrubael = {}
###
if erledigt_two == "FALSE":
	machinelearning_savefile = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+test_unip+"-"+test_gn+"/machine_learning_results.txt"
else:
	machinelearning_savefile = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+test_unip+"-"+test_gn+"-Custom/machine_learning_results.txt"

try:
	vectors =  vectormaker(outfile_name)
	for vect in vectors:
		position 	= vect.split("-")[4].replace(" ","").replace("\n","")
		seq             = vect.split("-")[6].replace(" ","").replace("\n","")
                typus           = vect.split("-")[7].replace(" ","").replace("\n","")
		evidence	= vect.split("-")[8].replace(" ","").replace("\n","")
		variant			= vect.split("-")[9].replace(" ","")
		mutscore		= AAScorer(variant)
                protclusta      = vect.split("-")[10].replace(" ","").replace("\n","")
                clustmemb       = vect.split("-")[11].replace(" ","").replace("\n","")
                evid_total      = vect.split("-")[12].replace(" ","").replace("\n","")
                evid_clusta     = vect.split("-")[13].replace(" ","").replace("\n","")
		bayes		= str(round(float(bayesian_scoring(seq,typus,protclusta,clustmemb,evid_total,evid_clusta, mutscore)),2))
		if databank == "Uniprot":
			randomforestvector = [seq,typus,evidence, protclusta,clustmemb,evid_total,evid_clusta,mutscore]
			rfcmd = "/usr/bin/python3 /net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/MachineLearning_Test/DevVersion_ML_Proteorizer_AppPredictor.py "+str(randomforestvector).replace(" ","")+" "+machinelearning_savefile+" "+str(position)+" "+str(clustermethod)
			try:
				os.system(rfcmd)
			except:
				pass
		if databank != "Uniprot":
			bayes = 0.0
		if asdrubael.has_key(position)==False:
			asdrubael[position]=[bayes]
except:
	logging.exception("message")
	bayes = "error"
	try:
		if asdrubael.has_key(position)==False:
			asdrubael[position]=[bayes]
	except:
		pass
	pass
try:
	with open(machinelearning_savefile,"r") as ml_result:
		for line in ml_result:
			position 	= line.split("\t")[4].replace(" ","").replace("\n","")
			classified 	= line.split("\t")[1].replace(" ","").replace("\n","")
			score		= line.split("\t")[3].replace(" ","").replace("\n","")
			if "1" in classified:
				rfresult = str(round(float(score),2))
			else:
				rfresult = str(round(float(score),2))
			asdrubael[position].append(rfresult)
except:
			pass
###
if "RW" in clustermethod:
	filus_addendum = "FinalResults_unlim_conservation_scored.txt"
elif "HClust" in clustermethod:
	filus_addendum = "FinalResults_unlim_conservation_scored_hclust.txt"
else:
	filus_addendum = "FinalResults_unlim_conservation_scored.txt"

if erledigt_two == "FALSE":
	outfile_name_FINAL = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+test_unip+"-"+test_gn+"/"+filus_addendum
else:
	outfile_name_FINAL = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+test_unip+"-"+test_gn+"-Custom/"+filus_addendum
###
try:
	with open(outfile_name,"r") as datafile:
		outfile = open(outfile_name_FINAL,"w")
		sys.stdout = outfile
		for line in datafile:
			if "Clusternumber" not in line:
				if "#" not in line:
					printed = "no"
					cluster = line.split("\t")[0].replace(" ","")
					typus = line.split("\t")[1].replace(" ","")
					position = line.split("\t")[2].replace(" ","")
					if "Input" in typus: 
						mutat = line.split("\t")[3].replace(" ","")
						try:
							cosmiccnt = max(extra_dictionary["COSMIC"][mutat])
						except:
							cosmiccnt = 0
						try:
							hetcnt = extra_dictionary["Hereditary"][mutat][0].split("-")[1]
							homcnt = extra_dictionary["Hereditary"][mutat][0].split("-")[2]
						except:
							hetcnt = 0
							homcnt = 0	
						combined_score	= 0

						try:
							bayes_score = float(asdrubael[str(position)][0])
						except:
							bayes_score = 0.0
						try:
							rf_Score = float(asdrubael[str(position)][1])
						except:
							rf_Score = 0.0
					
						### crafting the combined score
						if bayes_score >= 15:
							combined_score += 0.375
						elif bayes_score < 15:
							if bayes_score >= 12.5:
								combined_score += 0.25
							else:
								if bayes_score >= 10:
									combined_score += 0.125
						if rf_Score >= 0.9:
							combined_score += 0.375
						elif rf_Score < 0.9:
							if rf_Score >= 0.7:
								combined_score += 0.25
							else:
								if rf_Score >= 0.5:
									combined_score += 0.125
						### evaluating the combined score
						if combined_score <= 0.125:
							verdict = "No_Impact(high)"	
						elif combined_score > 0.125:	 										
							if combined_score <= 0.25:
								verdict = "No_Impact(medium)"
							if combined_score > 0.25:
								if combined_score < 0.35:
									verdict = "Ambiguous"	### this should cover the area of low confidence
							if combined_score >= 0.35:
								if combined_score <= 0.5:
									verdict = "Impact(medium)"
							if combined_score > 0.5:
								verdict = "Impact(high)"

						if databank != "Uniprot":
							verdict = "Not_scored"	

						if truevariants != "":
							positioninquestion = nums.search(str(mutat)).group(0)
							for item in truevariants.split(","):
								try:
									positiontocheck = nums.search(str(item)).group(0)
								except:
									positiontocheck = "-100"
								if positioninquestion == positiontocheck:
									printed = "yes"
									newstringler = "CAUTION: Probably wrong input ("+str(mutat)+"). We found "+str(item)+" instead. Please check."
									print line.replace("\n","").replace(mutat,newstringler),"\t",str(bayes_score),"\t",cosmiccnt,"\t",hetcnt,"\t",homcnt,"\t",str(rf_Score),"\t",verdict
							if printed == "no":
								print line.replace("\n",""),"\t",str(bayes_score),"\t",cosmiccnt,"\t",hetcnt,"\t",homcnt,"\t",str(rf_Score),"\t",verdict
						else:
							print line.replace("\n",""),"\t",str(bayes_score),"\t",cosmiccnt,"\t",hetcnt,"\t",homcnt,"\t",str(rf_Score),"\t",verdict

					else:
						mutat = line.split("\t")[2].replace(" ","")
						try:
							cosmiccnt = max(extra_dictionary["COSMIC"][mutat])
						except:
							cosmiccnt = "-"
						try:
							hetcnt = extra_dictionary["Hereditary"][mutat][0].split("-")[1]
							homcnt = extra_dictionary["Hereditary"][mutat][0].split("-")[2]
						except:
							hetcnt = "-"
							homcnt = "-"	
						print line.replace("\n",""),"\t","-","\t",cosmiccnt,"\t",hetcnt,"\t",homcnt,"\t","-","\t","-"
			else:
 				print "Clusternumber","\t","Data_Source","\t","Position","\t","Functional_Information","\t","Mechismo_Predictions","\t","SeqIdent%","\t","TypeIdent%","\t","BAYES","\t","COSMIC","\t","gnomAD_Het","\t","gnomAD_Hom","\t","Predictor","\t","Verdict"
except:
	logging.exception("message")
	pass

sys.stdout = orig_stdout
if project == "R_Submissions":
	if r_counter != "none":
		if erledigt_two == "FALSE":
			old_dir = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+test_unip+"-"+test_gn+"/"
			newdir = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+str(r_counter)+"-"+test_unip+"-"+test_gn+"/"


		else:
			old_dir = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+test_unip+"-"+test_gn+"-Custom/"
			newdir = "/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/PDB_Structures/"+project+"/"+str(r_counter)+"-"+test_unip+"-"+test_gn+"/"

		os.rename(old_dir,newdir)


end_time_15 = time.time()
differ_15 = end_time_15 - start_time_15
print "Step 15 time: ",differ_15
########################################################### 		####	#   #	###			###########################################################
###########################################################		#	##  #	#  #	  		###########################################################
###########################################################		###	# # #	#   #			###########################################################
###########################################################		#       #  ##	#  #	 		###########################################################
###########################################################		####    #   #	###	  		###########################################################
end_time = time.time()
print(end_time - start_time)	
