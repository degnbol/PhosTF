#!/usr/bin/env python

"""
NetworKIN(tm), (C) 2005,2006,2007,2013.
Drs Rune Linding, Lars Juhl Jensen, Heiko Horn & Jinho Kim

Usage: ./networkin.py Organism FastaFile SitesFile

If no sites file is given NetworKIN will predict on all T/S/Y residues 
in the given sequences.
"""

### Coding strategy & TODO
# FIX: filtering when no sites given
# 2) add code for prediction of downstream recognition site (14-3-3 first)
# 3) add code for predicting downstream binding module (14-3-3 first)
# 4) finnish up output format and CLI options
# 5) add code for doing 'pathway' analysis
# 6) predict phenotypes
###

# Changelog
# 21.10.06: tested blast -a 8 option, no differences on proteome
# 21.10.06: Filtercode broken
#			-fixed
# 28.01.07: Testing autophosphorylation (self == 1 in update script)
# 30.07.07: Working on v1.5 milestone
# 07.05.08: Initiated 2.5 w. scaling factor
# 11.07.13: New scoring scheme (Bayesian)
#

import sys, os, subprocess, fpformat, re, tempfile, dircache, random, operator, glob
import thread, threading
from optparse import OptionParser
from string import *
from likelihood import *

#debugging
import time

# Weighting parameter, 0=only motif, 1=only STRING
# estimated while benchmarking an then hardcoded here
# feel free to play with it, but it is your own responsibility
ALPHAS = {"9606": 0.85, "4932": 0.65}
dSpeciesName = {"9606": "human", "4932": "yeast"}
dPenalty = {"9606": {"hub penalty": 100, "length penalty": 800}, "4932": {"hub penalty": 170, "length penalty": 1000}}

NETWORKIN_SITE_FILE = 1
PROTEOME_DISCOVERER_SITE_FILE = 2
MAX_QUANT_DIRECT_OUTPUT_FILE = 3
	

global options
# Temporary directory to store the files
#tempfile.tempdir= '/tmp'

# Location of NetworKIN input files (string network, alias file etc.)
#datadir = sys.argv[0].rsplit("/", 1)[0]+'/data'
#global DATADIR
#DATADIR = ""

# Number of threads used for NetPhorest and BLAST
# setting it to a high number can also help in case NetworKIN uses too much memory
# as NetPhorest result files will be read one after the other to save memory
# -> the more files, the less memory usage
#global NUMBER_OF_PROCESSES
#NUMBER_OF_PROCESSES = "1";

# should temporary files of netphorest be zipped?
# saves some diskspace and makes it sometimes faster (depending on the CPU/harddisk speed)
#global SAVE_DISKSPACE
#SAVE_DISKSPACE = "";

# Should the analysis be limited to specific trees?
# just comment it out if you want to predict on all
#limitTrees = ['SH2']
#limitTrees = ['KIN', "SH2", "PTP"]

################################################################################
#                                                                              #
#                             the code starts here                             #
#                                                                              #
################################################################################

class CSheet(list):
    pass

# Run system binary
def myPopen(cmd):
	try:
		pipe = subprocess.Popen(cmd, shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout = pipe.stdout.readlines()
	except:
		sys.stderr.write('ERROR executing: '+`cmd`+'\n')
		sys.exit()
	else:
		return stdout 

# Read sequences from fasta file
def readFasta(fastafile):
	id_seq = {}
	aminoacids = re.compile('[^ACDEFGHIKLMNPQRSTVWYXB]')
	data = fastafile.readlines()
	fastafile.close()
	seq = ''
	for line in data:
		if line[0] <> ';':
			if line[0] == '>':
				line = line.strip()
				id = line[1:].split(' ', 1)[0]
#				id = line[1:-1].split(' ', 1)[0]
				seq = ''
			else:
				seq += aminoacids.sub('', line)
				if len(seq) > 0:
					id_seq[id] = seq
		else:
			pass
	return id_seq

def CheckInputType(sitesfile):
	f = open(sitesfile, 'rU')
	line = f.readline()
	f.close()
	
	tokens = line.split()
	if len(tokens) == 3:
		return NETWORKIN_SITE_FILE
	elif len(tokens) == 2:
		return PROTEOME_DISCOVERER_SITE_FILE
	elif tokens[0] == "Proteins" and tokens[4] == "Leading":
		return MAX_QUANT_DIRECT_OUTPUT_FILE
	else:
		sys.stderr("Unknown format of site file")
		sys.exit()

# Read phosphorylation sites from tsv file
# id -> position -> residue
def readPhosphoSites(sitesfile):
	id_pos_res = {}
	f = open(sitesfile, 'rU')
	if f:
		data = f.readlines()
		f.close()
		for line in data:
			tokens = line.split('\t')
			id = tokens[0]
			try:
				pos = int(tokens[1])
			except:
				sys.stderr.write(line)
				raise
			try:
				res = tokens[2][:-1]
			except:
				res = ""
			if id in id_pos_res:
				id_pos_res[id][pos] = res
			else:
				id_pos_res[id] = {pos: res}
	else:
		sys.stderr.write("Could not open site file: %s" % sitesfile)
		sys.exit()
		
	return id_pos_res


def readPhosphoSitesProteomeDiscoverer(fastafile, sitesfile):
	# This function takes a tab-separated list of ProteinID \t Peptide, in ProteomeDiscoverer format (i.e. phosphosites listed as
	# lowercase), maps the peptide to the full-length protein sequence and conducts the NetworKIN search.
	#
	
	# store all fasta sequences
	fasta = open(fastafile).readlines()
	fastadict = {}
	for line in fasta:
	    if line[0]==">":                                         #use this as each new entry in FASTA file starts with '>'
		ensp=line.split()[0].strip(">")[0:15]                     #strip the '>' symbol from the line before copying value into dictionary
		#print ensp
		fastadict[ensp]=""
	#retrieve the AA sequence belonging to individual ESPN identifiers and store as value for ESPN key       
	    else:
		seq=line.strip()
		fastadict[ensp]=seq
		
	# store all peptides in dictionary, with parent protein as key
	peptides = open(sitesfile).readlines()
	peptidedict = {}
	for line in peptides:
	    line = line.strip()
	    tokens = line.split("\t")
	    protID = tokens[0]
	    peptide = tokens[1]
	    if protID in peptidedict.iterkeys():
		if peptide in peptidedict[protID]:
		    pass
		else:
		    peptidedict[protID][peptide] = ""
	    else:
		peptidedict[protID] = { peptide : "" }
	
	
	# now map peptide onto full length sequence, then get absolute phosphosite locations
	id_pos_res = {}
	for protID in peptidedict.iterkeys():
	    for peptide in peptidedict[protID]:
		
		# first make peptide upper case to match to FASTA sequence
		UPPERpeptide = peptide.upper()
		#print UPPERpeptide
		
		# get parent protein sequence
		sequence = fastadict[protID]
		peptideindex = sequence.index(UPPERpeptide)
		#print peptideindex
		
		x = 0
		for letter in peptide:
		    if letter.islower():
			if letter == "s" or letter =="t" or letter == "y":
			    
			    phoslocation = peptideindex + x + 1
			    phosresidue = sequence[phoslocation - 1]
			    
			    if protID in id_pos_res.iterkeys():
				id_pos_res[protID][phoslocation] = phosresidue
			    else:
				id_pos_res[protID] = { phoslocation : phosresidue }
				
			else:
			    # non phosho modification
			    pass
		    else:
			pass
		    
		    x += 1
	return id_pos_res
	    

# return a list of dictionaries
# Usage l[index][column name]
def ReadSheet(fname, offset = 0):
    l = CSheet()

    f = open(fname, 'rU')
    
    for i in range(offset):
        f.readline()

    columns = f.readline().strip().split('\t')

    l.columns = columns
    for line in f.readlines():
        instance = {}
        l.append(instance)
        fields = line.strip().split('\t')
        
        if len(fields) > len(columns):
            sys.stderr.write("Error in data file")
            sys.stderr.write(columns)
            sys.stderr.write(fields)
            raise
        
        for i in range(len(columns)):
            try:
                instance[columns[i]] = fields[i]
            except IndexError:
                instance[columns[i]] = ''
    
    #sys.stderr.write("No. of entries in file input: %d" % len(l))
    f.close()
    
    return l

def readPhosphoSitesMaxQuant(fname, only_leading = False):
    id_pos_res = {}
    
    phosphosites = ReadSheet(fname)
    
    for site in phosphosites:
        Ids = site["Proteins"].split(';')
        positions = map(lambda(x):int(x), site["Positions within proteins"].split(';'))
        aa = site["Amino acid"]
    
        leading_protein_ids = site["Leading proteins"].split(';')
    
        for i in range(len(Ids)):
            Id = Ids[i]
            pos = positions[i]

            if only_leading and not Id in leading_protein_ids:
                continue

            if Id in id_pos_res:
                id_pos_res[Id][pos] = aa
            else:
                id_pos_res[Id] = {pos: aa}

    return id_pos_res


#Alias hashes
def readAliasFiles(organism, datadir):
	alias_hash = {}
	desc_hash = {}

	# Read alias db
	try:
		alias_db = myPopen('gzip -cd %s/%s.alias_best.tsv.gz'%(datadir, organism))
		for line in alias_db:
			(taxID, seqID, alias) = line.strip().split('\t')[:3]
			alias_hash[seqID] = alias
	except:
		sys.stderr.write("No aliases available for organism: '%s'\n"%organism)

	# Read desc db
	try:
		desc_db = myPopen('gzip -cd %s/%s.text_best.tsv.gz'%(datadir, organism))
		for line in desc_db:
			(taxID, seqID, desc) = line.split('\t')[:3]
			desc_hash[seqID] = desc
	except:
		sys.stderr.write("No descriptions available for organism: '%s'\n"%organism)

	return alias_hash, desc_hash

# Run Netphorest
def runNetPhorest(id_seq, id_pos_res, save_diskspace, number_of_processes, number_of_active_processes = 1, fast = False, leave_intermediates = False):
	if number_of_active_processes < 2:
		raise "Number of maximum threads is less than 2"
	id_pos_tree_pred = {}

	#check how many sequences we actually have
	number_of_sequences = 0
	if id_pos_res == {}:
		number_of_sequences = len(id_seq)
	else:
		number_of_sequences = len(id_pos_res)

	if number_of_sequences < number_of_processes:
		number_of_processes = number_of_sequences

	# use multiple instances of netphorest
	file_in = range(number_of_processes)
	file_out = range(number_of_processes)
	class CDummy:
		def __init__(self, name):
			self.name = name
			
	# create filehandles
	for i in range(number_of_processes):
		file_in[i] = tempfile.NamedTemporaryFile()
		
		if fast or leave_intermediates:
			if save_diskspace:
				file_out[i] = CDummy("%s.%s.gz" % (fn_netphorest_output, str(i)))      # jhkim
				
			else:
				file_out[i] = CDummy("%s.%s.txt" % (fn_netphorest_output, str(i)))
		else:
			file_out[i] = tempfile.NamedTemporaryFile()

	# distribute data into different files
	line_counter = 0;
	if id_pos_res == {}:
		for id in id_seq:
			line_counter = line_counter + 1
			number = line_counter % number_of_processes
			file_in[number].write(">%s\n%s\n"%(id, id_seq[id]))
		#file = open(args[1], 'r')
		#for line in file:
		#	if( re.search('^>',line) ):
		#		line_counter = line_counter + 1
		#	number = line_counter % number_of_processes
		#	file_in[number].write(line)
	else:
		for id in id_seq:
			if id in id_pos_res:
				line_counter = line_counter + 1
				number = line_counter % number_of_processes
				file_in[number].write(">%s\n%s\n"%(id, id_seq[id]))
	for i in range(number_of_processes):
		file_in[i].flush()

	# run NetPhorest for each file in parallel
	if options.verbose:
		sys.stderr.write("\n")
		
	if options.verbose:
		sys.stderr.write("Running on %s sequences\n"%line_counter)
		
	for i in range(number_of_processes):
		if fast and os.path.isfile(file_out[i].name):
			if options.verbose:
				sys.stderr.write("%s is already exist.\n" % file_out[i].name)
			continue
				
		while threading.activeCount() >= number_of_active_processes:
			sys.stderr.write('.')
			time.sleep(5)
			
		if options.verbose:
			sys.stderr.write(netphorest_bin+' < '+file_in[i].name+ ' > '+file_out[i].name+'\n');
		if save_diskspace:
			arg = (netphorest_bin+' < '+file_in[i].name+ '| gzip -9 > '+file_out[i].name,)
		else:
			arg = (netphorest_bin+' < '+file_in[i].name+ ' > '+file_out[i].name,)
		
		threading.Thread(target=myPopen, args=arg).start()
	
	# wait for threads to finish
	while (threading.activeCount() > 1):
		sys.stderr.write('.')
		time.sleep(5)

	return file_out

# Parse the NetPhorest output
# expected format:
# Name  Position        Residue Peptide Method    Orgnism    Tree    Classifier      Posterior       Prior
# O00151  2       T       ----MtTQQID     nn    human    KIN     CDK2_CDK3_group    0.040050    0.028284
def parseNetphorestFile(filename, id_pos_res, save_diskspace):
	if save_diskspace:
		command = "gzip -cd %s"%(filename)
	else:
		command = "cat %s"%(filename)

	if options.verbose:
		sys.stderr.write('Parsing NetPhorest result file: "' + filename + '"\n');
	try:
		netphorest_results = myPopen(command)
	except:
		sys.stderr.write("Going to sleep for 1 hour to give you time to debug:\nCrashed with '%s'\n"%command)
		time.sleep(3600)

	id_pos_tree_pred = {}

	for line in netphorest_results:
		if(re.match('^#',line)):
			continue
		tokens = line.split('\t')
		id = tokens[0]
		pos = int(tokens[1])
		# NetPhorest2 introduces organism column
		try:
			(res, peptide, method, organism, tree, pred) = tokens[2:8]
		except ValueError:
			sys.stderr.write(tokens)
			raise
		try:
			if tree not in limitTrees:
				continue
		except:
			pass
		try:
			score = float(tokens[8])
		except:
			continue
		if (id in id_pos_res and pos in id_pos_res[id]) or id_pos_res == {}:
			if id in id_pos_tree_pred:
				if pos in id_pos_tree_pred[id]:
					if tree in id_pos_tree_pred[id][pos]:
						id_pos_tree_pred[id][pos][tree][pred] = (res,peptide,score)
					else:
						id_pos_tree_pred[id][pos][tree] = { pred: (res,peptide,score) }
				else:
					id_pos_tree_pred[id][pos]= { tree: { pred: (res,peptide,score) } }					
			else:
				id_pos_tree_pred[id]= { pos: { tree: { pred: (res,peptide,score) } } }
		else:
			pass
	else:
		pass
	return id_pos_tree_pred

def ReadLines(fname):
	f = open(fname)
	lines = f.readlines()
	f.close()
	return lines

def WriteString(fname, s):
    f = open(fname, 'w')
    f.write(s)
    f.close()
    
# Map incoming peptides to STRING sequences
def mapPeptides2STRING(blastDir, organism, fastafilename, id_pos_res, id_seq, number_of_processes, datadir, fast = False, leave_intermediates = False):
	sys.stderr.write("Mapping using blast\n")
	incoming2string = {}
	string2incoming = {}
	# 
	# # Speedup, only blast sequences with site specified
	# blast_tmpfile = tempfile.NamedTemporaryFile()
	# if id_pos_res == {}:
	# 	for id in id_seq:
	# 		blast_tmpfile.write('>'+id+'\n'+id_seq[id]+'\n')
	# else:
	# 	for id in id_pos_res:
	# 		try:
	# 			blast_tmpfile.write('>'+id+'\n'+id_seq[id]+'\n')
	# 		except:
	# 			sys.stderr.write("No sequence available for '%s'\n"%id)
	# blast_tmpfile.flush()
	# 
	# blastDB = os.path.join(datadir, "%s.protein.sequences.fa"%(organism)).replace(' ', '\\ ')
	# 
	# # Check if blast database is actually initialized, if not: do it
	# if not os.path.isfile(blastDB+'.pin'):
	# 	command = "%s/formatdb -i %s"%(blastDir.rsplit("/", 1)[0], blastDB)
	# 	sys.stderr.write("Looks like blast database is not initialized, trying to run:\n%s\n"%command)
	# 	myPopen(command)
	# 
	# command = "%s -a %s -p blastp -e 1e-10 -m 8 -d %s -i %s | sort -k12nr"%(blastDir, number_of_processes, blastDB, blast_tmpfile.name)
	# 
	# to save time - jhkim
	# if fast and os.path.isfile(fn_blast_output):
	blast_out = ReadLines(fn_blast_output)
	# else:
	# 	blast_out = myPopen(command)
	# 	if leave_intermediates:
	# 		WriteString(fn_blast_output, "".join(blast_out))
		
	for line in blast_out:
		tokens = line.split('\t')
		incoming = tokens[0]
		if incoming not in incoming2string:
			# get rid of organism prefix
			string = tokens[1].replace("%s."%organism, "")
			identity = float(tokens[2])
			evalue = float(tokens[-2])

			if string in string2incoming:
				sys.stderr.write('Best hit for '+incoming+' is not reciprocal\n')
			if identity < 90:
				sys.stderr.write('Best hit for '+incoming+' has only '+fpformat.fix(identity, 2)+' %identity\n')
			if evalue > 1e-40:
				sys.stderr.write('Best hit for '+incoming+' has high E-value '+fpformat.sci(evalue, 2)+' \n')
			if incoming in incoming2string:
				incoming2string[incoming][string] = True
			else:
				incoming2string[incoming] = { string: True }
			if string in string2incoming:
				string2incoming[string][incoming] = True
			else:
				string2incoming[string] = { incoming: True }
		else:
			pass
	return incoming2string, string2incoming

# Random mapping of incoming identifiers to string
def mapRandom(id_seq):
	sys.stderr.write("Mapping random\n")
	import random

	incoming2string = {}
	string2incoming = {}

	stringIDs = []
	file = open('%s/%s.protein.sequences.fa'%(datadir,organism), 'r')
	data = file.readlines()
	file.close()

	for line in data:
		if line[0] == '>':
			name = line[1:-1]
			stringIDs.append(name)
			string2incoming[name] = {}
	max = len(stringIDs) - 1

	for incoming in id_seq:
		if incoming not in incoming2string:
			int = random.randint(0, max)
			string = stringIDs[int]
			incoming2string[incoming] = { string: True }
			if string in string2incoming:
				string2incoming[string][incoming] = True
			else:
				string2incoming[string] = { incoming: True }
	return incoming2string, string2incoming

# In case we run NetworKIN on the same sequence set we use in STRING, we can skip the mapping by blast
def mapOne2one(id_seq):
	sys.stderr.write("Mapping one2one\n")
	incoming2string = {}
	string2incoming = {}
	for incoming in id_seq:
		incoming2string[incoming] = {incoming: True}
		string2incoming[incoming] = {incoming: True}
	return incoming2string, string2incoming

# In case we have a better mapping then we can expect from blasting, we can use an external file to do so
# file format:
# incoming ID -> STRING ID
def mapFromFile(filename):
	sys.stderr.write("Mapping using external mapping file\n")
	incoming2string = {}
	string2incoming = {}

	command = "cat %s"%(filename)
	try:
		mappingFile = myPopen(command)
	except:
		sys.stderr.write("Going to sleep, crashed with '%s'\n"%command)
		time.sleep(3600)

	for line in mappingFile:
		if(re.match('^#',line)):
			continue
	
		line = line.strip()
		tokens = line.split('\t')

		incoming = tokens[0]
		string = tokens[1]

		if incoming in incoming2string:
			incoming2string[incoming][string] = True
		else:
			incoming2string[incoming] = { string: True }
		if string in string2incoming:
			string2incoming[string][incoming] = True
		else:
			string2incoming[string] = { incoming: True }

	return incoming2string, string2incoming

# Load the precalculated STRING network file
def loadSTRINGdata(string2incoming, datadir, number_of_processes):
	#command = 'gzip -cd %s/%s.bestpath.tsv.gz'%(datadir, organism)
	
	fn_bestpath = "%s/%s.string_000_%04d_%04d.tsv.gz" % (os.path.join(datadir, "bestpath"), organism, dPenalty[organism]["hub penalty"], dPenalty[organism]["length penalty"])
	if not os.path.isfile(fn_bestpath):
		sys.stderr.write("Best path file does not exist: %s" % fn_bestpath)
				 
	'''	    
	command = "gzip -cd %s" % fn_bestpath
	
	try:
		data = myPopen(command)
	except:
		sys.stderr.write("Error loading STRING data using '%s', sleeping fo 1h.\n"%command)
		time.sleep(3600)
	'''	    

	tree_pred_string_data = {}

		
	# for memory efficiency
	import gzip
	f = gzip.open(fn_bestpath)
	#for line in data:
	while (True):
		line = f.readline()
		if line == '':
			break
		line = line.strip()
		tokens = line.split('\t')
		if len(tokens) == 8:
			(tree, group, name, string1, string2, stringscore, stringscore_indirect, path) = tokens
		elif len(tokens) == 7:
			(tree, group, name, string1, string2, stringscore, stringscore_indirect) = tokens
			path = ""	# path to itself,  we will miss the path information
		elif len(tokens) == 6:
			(name, string1, string2, stringscore, stringscore_indirect, path) = tokens
		elif len(tokens) == 5:
			(name, string1, string2, stringscore, stringscore_indirect) = tokens
			path = ""	# path to itself,  we will miss the path information


		if string2 in string2incoming:
			if string2 in tree_pred_string_data:
				tree_pred_string_data[string2][string1] = {"_name": name}
			else:
				tree_pred_string_data[string2] = {string1: {"_name": name}}

			if options.path == "direct":
				tree_pred_string_data[string2][string1]["_score"] = float(stringscore)
			elif options.path == "indirect":
				tree_pred_string_data[string2][string1]["_score"] = float(stringscore_indirect)	# Use indirect path
			else:
				raise "Path information should be either direct or indirect."
			tree_pred_string_data[string2][string1]["_path"] = path
		else:
			pass

	f.close()
	return tree_pred_string_data

def InsertValueIntoMultiLevelDict(d, keys, value):
    for i in range(len(keys)-1):
        if not d.has_key(keys[i]):
            d[keys[i]] = {}
        d = d[keys[i]]

    if not d.has_key(keys[-1]):
        d[keys[-1]] = []
    d[keys[-1]].append(value)
    
def ReadGroup2DomainMap(path_group2domain_map):
    map_group2domain = {}   # KIN   group   name
    
    f = open(path_group2domain_map, "rU")

    for line in f.readlines():
        tokens = line.split()
        InsertValueIntoMultiLevelDict(map_group2domain, tokens[:2], tokens[2])

    f.close()
    
    return map_group2domain

def SetValueIntoMultiLevelDict(d, keys, value):
    for i in range(len(keys)-1):
        if not d.has_key(keys[i]):
            d[keys[i]] = {}
        d = d[keys[i]]
    if d.has_key(keys[-1]) and type(d[keys[-1]]) != type(value):
        sys.stderr.write("Causion: multi-dict already has value and try to assign a value of different type")
        pass
    if d.has_key(keys[-1]):
        if d[keys[-1]] != value:
            sys.stderr.write("This operation replaces a value (%s)" % " ".join(map(lambda(x):str(x), keys)))
    d[keys[-1]] = value
    
def printResult(id_pos_tree_pred, tree_pred_string_data, incoming2string, string_alias, string_desc, organism, mode, dir_likelihood_conversion_tbl, map_group2domain):
	ALPHA = ALPHAS[organism]
	species = dSpeciesName[organism]

	dLRConvTbl = {}
	for fname in glob.glob(os.path.join(dir_likelihood_conversion_tbl, "conversion_tbl_*_smooth*")):
		netphorest_or_string, species_of_conversion_table, tree, player_name = re.findall("conversion_tbl_([a-z]+)_smooth_([a-z]+)_([A-Z0-9]+)_([a-zA-Z0-9_/-]+)", os.path.basename(os.path.splitext(fname)[0]))[0]
		#species, tree, player_name = os.path.basename(os.path.splitext(fname)[0]).rsplit('_', 3)[1:]
		
		if species_of_conversion_table != species:
			continue
		
		conversion_tbl = ReadConversionTableBin(fname)
		SetValueIntoMultiLevelDict(dLRConvTbl, [species_of_conversion_table, tree, player_name, netphorest_or_string], conversion_tbl)
    
		if options.verbose:
			sys.stderr.write("Conversion table %s %s %s %s\n" % (species_of_conversion_table, tree, player_name, netphorest_or_string) )
    
	# For each ID in NetPhorest
	for id in id_pos_tree_pred:
		# We have a mapping to STRING
		if id in incoming2string:
			# For each predicted position
			for pos in id_pos_tree_pred[id]:
				# For each of the trees (KIN, SH@ etc.)
				for tree in id_pos_tree_pred[id][pos]:
					score_results = {}
					# For each single classifier
					for pred in id_pos_tree_pred[id][pos][tree]:
						# For each mapped sequence
						for string1 in incoming2string[id]:
							if string1 in string_alias:
								bestName1 = string_alias[string1]
							else:
								bestName1 = ''
							if string1 in string_desc:
								desc1 = string_desc[string1]
							else:
								desc1 = ''
							if string1 in tree_pred_string_data:
								(res, peptide, netphorestScore) = id_pos_tree_pred[id][pos][tree][pred]
								for string2 in tree_pred_string_data[string1]:
									
									if string2 in string_alias:
										bestName2 = string_alias[string2]
									else:
										bestName2 = ''
									if string2 in string_desc:
										desc2 = string_desc[string2]
									else:
										desc2 = ''
									stringScore = tree_pred_string_data[string1][string2]["_score"]
									path = tree_pred_string_data[string1][string2]["_path"]
									name = tree_pred_string_data[string1][string2]["_name"]	# string2 = kinase
									
									#sys.stderr.write("%s %s %s\n" % (tree, pred, name))
									if not map_group2domain.has_key(tree) or not map_group2domain[tree].has_key(pred) or not name in map_group2domain[tree][pred]:
										if options.string_for_uncovered:
											if species == "human":
												if tree in ["1433", "BRCT", "WW", "PTB", "WD40", "FHA"]:
													conversion_tbl_string = dLRConvTbl[species]["SH2"]["general"]["string"]
												else:
													conversion_tbl_string = dLRConvTbl[species][tree]["general"]["string"]
											elif species == "yeast":
												conversion_tbl_string = dLRConvTbl[species][tree]["general"]["string"]
											else:
												raise "This species is not supported"
											
											likelihood_netphorest = 1
											likelihood_string = ConvertScore2L(stringScore, conversion_tbl_string)
											unified_likelihood = likelihood_netphorest * likelihood_string
											networkinScore = unified_likelihood
		
											# NetworKIN result
											result = id+'\t'+res+str(pos)+'\t'+tree+'\t'+pred+'\t'+name+'\t'+ \
											fpformat.fix(networkinScore,4)+'\t'+"N/A"+'\t'+fpformat.fix(stringScore,4)+'\t'+ \
											string1+'\t'+string2+'\t'+bestName1+'\t'+bestName2+'\t'+desc1+'\t'+desc2+'\t'+peptide+'\t'+path+'\n'
										else:
											continue
									else:
										#sys.stderr.write("%s %s\n" % (string1, string2))
	
										if species == "human":
											if tree in ["1433", "BRCT", "WW", "PTB", "WD40", "FHA"]:
												conversion_tbl_netphorest = dLRConvTbl[species]["SH2"]["general"]["netphorest"]
												conversion_tbl_string = dLRConvTbl[species]["SH2"]["general"]["string"]
											else:
												if dLRConvTbl[species][tree].has_key(name):
													conversion_tbl_netphorest = dLRConvTbl[species][tree][name]["netphorest"]
													conversion_tbl_string = dLRConvTbl[species][tree][name]["string"]
												else:
													conversion_tbl_netphorest = dLRConvTbl[species][tree]["general"]["netphorest"]
													conversion_tbl_string = dLRConvTbl[species][tree]["general"]["string"]
										elif species == "yeast":
											if dLRConvTbl[species][tree].has_key(name):
												conversion_tbl_netphorest = dLRConvTbl[species][tree][name]["netphorest"]
												conversion_tbl_string = dLRConvTbl[species][tree][name]["string"]
											else:
												conversion_tbl_netphorest = dLRConvTbl[species][tree]["general"]["netphorest"]
												conversion_tbl_string = dLRConvTbl[species][tree]["general"]["string"]
										else:
											raise "This species is not supported"
										
										likelihood_netphorest = ConvertScore2L(netphorestScore, conversion_tbl_netphorest)
										likelihood_string = ConvertScore2L(stringScore, conversion_tbl_string)
										unified_likelihood = likelihood_netphorest * likelihood_string
										networkinScore = unified_likelihood
										#networkinScore = pow(stringScore, ALPHA)*pow(netphorestScore, 1-ALPHA)
	
										# NetworKIN result
										result = id+'\t'+res+str(pos)+'\t'+tree+'\t'+pred+'\t'+name+'\t'+ \
										fpformat.fix(networkinScore,4)+'\t'+fpformat.fix(netphorestScore,4)+'\t'+fpformat.fix(stringScore,4)+'\t'+ \
										string1+'\t'+string2+'\t'+bestName1+'\t'+bestName2+'\t'+desc1+'\t'+desc2+'\t'+peptide+'\t'+path+'\n'

									if networkinScore not in score_results:
										score_results[networkinScore] = []

									score_results[networkinScore].append(result)
					if mode == 'network':
						highestScore = sorted(score_results.keys(), reverse=True)[0]

						if len(score_results[highestScore]) > 1:
							index = random.randint(0,len(score_results[highestScore])-1)
							sys.stdout.write((score_results[highestScore][index]))
						else:
							sys.stdout.write((score_results[highestScore][0]))
						pass
					else:
						for score in sorted(score_results.keys(), reverse=True):
							sys.stdout.write("".join(score_results[score]))
		else:
			pass
	return

#MAIN
def Main():
	sys.stderr.write("Reading fasta input file\n")
	id_seq = readFasta(fastafile)
	if options.verbose:
		sys.stderr.write("%s sequences loaded\n"%len(id_seq.keys()))

	if sitesfile:
		sys.stderr.write("Reading phosphosite file\n")
		input_type = CheckInputType(sitesfile)
		if input_type == NETWORKIN_SITE_FILE:
			id_pos_res = readPhosphoSites(sitesfile)
		elif input_type == PROTEOME_DISCOVERER_SITE_FILE:
			id_pos_res = readPhosphoSitesProteomeDiscoverer(fn_fasta, sitesfile)
		elif input_type == MAX_QUANT_DIRECT_OUTPUT_FILE:
			id_pos_res = readPhosphoSitesMaxQuant(sitesfile)
	else:
		id_pos_res = {}

	sys.stderr.write("Loading aliases and descriptions\n")
	(string_alias, string_desc) = readAliasFiles(args[0], options.datadir);

	if organism == "9606":
		path_group2domain_map = os.path.join(options.datadir, "group_human_protein_name_map.tsv")
	elif organism == "4932":
		path_group2domain_map = os.path.join(options.datadir, "group_yeast_KIN.tsv")
	

	map_group2domain = ReadGroup2DomainMap(path_group2domain_map)
	
	# Default way of mapping using BLAST
	incoming2string, string2incoming = mapPeptides2STRING(blastDir, organism, fastafile.name, id_pos_res, id_seq, options.threads, options.datadir, options.fast, options.leave)

	# Hack for random mapping to proteins
	#incoming2string, string2incoming = mapRandom(id_seq)

	# Use if a mapping file for the input can be provided
	#incoming2string, string2incoming = mapFromFile("/home/red1/hhorn/projects2/2012_03_22_Jesper/ensmusp_ensp.tsv")

	# Used if only ensembl of the same version used
	#incoming2string, string2incoming = mapOne2one(id_seq)

	# Load the STRING network data
	sys.stderr.write("Loading STRING network\n")
	tree_pred_string_data = loadSTRINGdata(string2incoming, options.datadir, options.threads)

	# Run NetPhorest
	sys.stderr.write("Running NetPhorest")
	netphorestTmpFiles = runNetPhorest(id_seq, id_pos_res, options.compress, options.threads, options.active_threads, options.fast, options.leave)
	sys.stderr.write('\n')

	# Writing result to STDOUT
	sys.stderr.write("Writing results\n")
	sys.stdout.write("#Name\tPosition\tTree\tNetPhorest Group\tKinase/Phosphatase/Phospho-binding domain\tNetworKIN score\tNetPhorest probability\tSTRING score\tTarget STRING ID\tKinase/Phosphatase/Phospho-binding domain STRING ID\tTarget description\tKinase/Phosphatase/Phospho-binding domain description\tTarget Name\tKinase/Phosphatase/Phospho-binding domain Name\tPeptide sequence window\tIntermediate nodes\n")
	for i in range(len(netphorestTmpFiles)):
		id_pos_tree_pred = parseNetphorestFile(netphorestTmpFiles[i].name, id_pos_res, options.compress)
		if options.path == "direct":
			dir_likelihood_conversion_tbl = os.path.join(options.datadir, "likelihood_conversion_table_direct")
		elif options.path == "indirect":
			dir_likelihood_conversion_tbl = os.path.join(options.datadir, "likelihood_conversion_table_indirect")
		else:
			raise "Path information should be either direct or indirect."
		printResult(id_pos_tree_pred, tree_pred_string_data, incoming2string, string_alias, string_desc, args[0], options.mode, dir_likelihood_conversion_tbl, map_group2domain)

	return

if __name__ == '__main__':
	#BLAST
	try:
		blastDir = os.environ['BLAST_PATH']
	except:
		blastDir=""
	#NETPHOREST
	try:
		netphorest_bin = os.environ['NETPHOREST_PATH']
	except:
		netphorest_bin=""

	usage = "usage: %prog [options] organism FASTA-file [sites-file]"
	parser = OptionParser(usage=usage, version="%prog 3.0")
	parser.add_option("-n", "--netphorest", dest="netphorest_bin", default=netphorest_bin,
										help="set the location of the NetPhorest binary, overwrites the 'NETPHOREST_PATH' environmental variable. [ENV: %default]")
	parser.add_option("-b", "--blast", dest="blast", default=blastDir,
										help="set the directory for the BLAST binaries (formatdb and blastall), overwrites the 'BLAST_PATH' environmental variable. [ENV: %default]")
	parser.add_option("-m", "--mode", dest="mode", default=False,
										help="if set to 'network', gives only one best scoring result for each site. In case of multiple candidate kinases with the same core, the selection hapens randomly. [default: %default]")
	parser.add_option("-p", "--path", dest="path", default="direct",
										help="NetworKIN uses both direct and indirect paths. Otherwise, it uses only indirect paths. [default: %default]")
	parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
										help="print out everything [default: %default]")
	parser.add_option("-f", "--fast", dest="fast", default=False, action="store_true",
										help="Speed up by using the intermediate files of previous run [default: %default]")
	parser.add_option("-l", "--leave", dest="leave", default=False, action="store_true",
										help="leave intermediate files [default: %default]")
	parser.add_option("-u", "--uncovered", dest="string_for_uncovered", default=False, action="store_true",
										help="Use STRING likelihood for uncovered Kinases [default: %default]")
	parser.add_option("-t", "--threads", dest="threads", default=1, type="int",
										help="number of available threads/CPUs. Also leads to less memory usage, as result files are read sequentially [default: %default]")
	parser.add_option("--nt", dest="active_threads", default=2, type="int",
										help="number of active threads at a time")
	
	parser.add_option("-c", "--compress", dest="compress", default=True,
										help="compress temporary result files, saves discspace [default: %default]")
	parser.add_option("-d", "--data", dest="datadir", default=os.path.join(os.path.split(os.path.realpath(sys.argv[0]))[0], 'data'),
										help="location for the additional files like the pre-computed STRING network, STRING sequence database etc. [default: %default]")
	parser.add_option("--blast-result")
	
	if os.environ.has_key("TMPDIR"):
		parser.add_option("--tmp", dest="tmpdir", default=os.environ["TMPDIR"],
										help="location for the temporary files [default: %default]")
	else:
		print >> sys.stderr, "TMPDIR environmental variable is not defined. Please define this variable or specify tmp directory by using --tmp command line option"
		sys.exit()


	global options
	(options, args) = parser.parse_args()

	if options.active_threads < 2:
		parser.error("Number of active thread (--nt) is less than 2")
	tempfile.tempdir= options.tmpdir

	#ORGANISM
	try:
		organism = args[0]
	except:
		parser.error("Organism not defined!")

	#SEQUENCE FILE
	try:
		fn_fasta = args[1]
		if options.blast_result is not None:
			fn_blast_output = options.blast_result
		else:
			fn_blast_output = "%s.%s.blast.out" % (fn_fasta, organism)
			
		fn_netphorest_output = "%s.%s.netphorest.out" % (fn_fasta, organism)
		fastafile = open(fn_fasta, 'rU')
	except:
		sys.stderr.write("%s"%args)
		parser.error("FASTA-file not defined!")

	#SITES FILE
	try:
		sitesfile = args[2]
	except:
		sitesfile = False

	#BLAST
	if options.blast:
		blastDir = options.blast
	
	#NETPHOREST
	if options.netphorest_bin:
		netphorest_bin = options.netphorest_bin
		
	# Show runtime parameters
	if(options.verbose):
		sys.stderr.write('\nPredicting using parameters as follows:\nOrganism:\t%s\nFastaFile:\t%s\n'%(organism, fn_fasta))
		sys.stderr.write('Threads:\t%s\nCompress:\t%s\n'%(options.threads, options.compress))
		if options.string_for_uncovered:
			sys.stderr.write("Use STRING likelihood when the kinases is not covered NetPhorest.\n")
		if sitesfile:
			sys.stderr.write('Sitesfile:\t%s\n'%sitesfile)
		else:
			sys.stderr.write('No sites-file given, predicting on all S,T,Y residues.\n')
		if options.mode:
			sys.stderr.write('Mode:\t\t%s\n'%options.mode)
		sys.stderr.write("Blast dir: %s" % blastDir)
		sys.stderr.write("NetPhorest binary: %s" % netphorest_bin)
		sys.stderr.write('\n')
	Main()
