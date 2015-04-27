#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File   		: BlobCollection.py
Version 	: 0.1
Author 		: Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Bugs 		: ?
To do 		: 
"""

from __future__ import division
from MiscFunctions import n50, keyWithMaxVal
import numpy
import time
import subprocess
import commands
import re
import sys

class Blob():
	'''
	Blobs contain the information of a contig in an assembly.
	Each line that eventually gets printed into the blobplot file is a blob. 
	'''

	def __init__(self, header, seq):
		self.name = header
		self.seq = seq.upper()
		self.length = len(seq)
		self.corrected_length = len(self.seq) - self.seq.count('N') 
		self.gc = float((self.seq.count('G') + self.seq.count('C') ) / self.corrected_length ) if self.corrected_length > 0 else 0.0
		self.covs = dict()
		self.tax = dict() 

class BlobCollection():
	'''
	The BlobCollection class contains the blobs and all methods needed 
	for parsing the input files and doing the necessary computations 
	for arriving at the final results. 
	'''
	def __init__(self, target_ranks, rank_level):
		self.contigs = dict()
		self.index = list()
		self.outfile = str()
		self.stats = dict()
		self.cov_libs = list()
		self.blast_libs = list()
		self.target_ranks = target_ranks # this one could be removed, i think
		self.rank_level = rank_level
		self.blast_order = list()

	def addBlob(self, blob):
		'''
		Adds a blob to the BlobCollection
			- puts it in the contigs dict with header as the key
			- puts it header in the index list so that it can be found by index (this is for CAS files and printing blobs in order)
		'''
		if not blob.name in self.contigs: 
			self.contigs[blob.name] = blob
		else: 
			sys.exit("[ERROR] - Sequence header {} occurs more than once".format(blob.name))
		self.index.append(blob.name)

	def getBlobsFromAssembly(self, assembly_file, assembly_type, exclude_assembly_cov):
		print "[STATUS] - Parsing assembly %s" % (assembly_file)
		header, seq = '', ''
		with open(assembly_file) as fh:
			for line in fh:
				line_data = line.rstrip("\n")
				if line.startswith('>'):
					if (seq): 
						blob = Blob(header, seq) 
						self.addBlob(blob)
						seq = ''
						if assembly_type == 'unknown' or exclude_assembly_cov:
							pass
						else:
							cov = float(self.parseCovFromHeader(header, assembly_type))
							self.addBlobCov(header, assembly_type, cov)
					header = line.rstrip("\n").lstrip(">")
				else:
					seq += line.rstrip("\n").upper() 
			blob = Blob(header, seq) 
			self.addBlob(blob)
			if assembly_type == 'unknown' or exclude_assembly_cov:
				pass
			else:
				cov = float(self.parseCovFromHeader(header, assembly_type))
				self.addBlobCov(header, assembly_type, cov)
		if not assembly_type == 'unknown' and not exclude_assembly_cov:
			self.cov_libs.append(assembly_type)

	def addBlobCov(self, header, mapping_lib, cov):
		'''
		Adds coverage to the covs dict: key = cov_lib name, value = coverage
		'''
		if mapping_lib in self.contigs[header].covs:
			sys.exit("[ERROR] - Contig {} received more than one coverage from the {} assembly file".format( header, mapping_lib))
		else:
			self.contigs[header].covs[mapping_lib] = float(cov)

	def parseCovFromHeader(self, header, assembly_type):
		''' 
		Returns the coverage from the header of a FASTA 
		sequence depending on the assembly type
		'''
		if assembly_type == 'spades':
			return header.split("_")[-3]
		elif assembly_type == 'velvet':
			return header.split("_")[-1]
		elif assembly_type == 'abyss':
			temp = header.split(" ")
			return temp[2]/(temp[1]+1-75)
		else:
			sys.exit("[ERROR] - Coverage could not be parsed from header {} of type {}".format( header, assembly_type ))

	def getCovForBlobs(self, mapping_files):
		'''
		- Coordinates parsing of coverages from different mapping files
		- If more than one cov_lib is specified then a new cov_lib is created for each blob which contain the sum of the other cov_libs 
		'''
		for cov_lib, mapping_file in mapping_files.items():
			self.cov_libs.append(cov_lib)
			print "[STATUS] - Parsing coverage from cov_lib \"{}\" from file {}".format(cov_lib, mapping_file)
			if cov_lib.startswith("CAS"):
				self.parseCovFromCasFile(cov_lib, mapping_file)
			elif cov_lib.startswith("BAM"):
				pass
			elif cov_lib.startswith("SAM"):
				pass
			elif cov_lib.startswith("COV"):
				self.parseCovFromCovFile(cov_lib, mapping_file)
			else:
				sys.exit("[ERROR] - Unknown cov_lib type {} in {}".format(cov_lib, mapping_file))	
	
		if (len(self.cov_libs) > 1):
			for contig_name, blob in self.contigs.items():
				cov_sum = 0.0
				for cov_lib in self.cov_libs:
					cov_sum += blob.covs[cov_lib]	
				self.contigs[contig_name].covs['SUM']=cov_sum		
			self.cov_libs.append("SUM")

	def parseCovFromCasFile(self, lib_name, cas_file):
		'''
		Parse coverage from CAS file
		'''
		error, message = commands.getstatusoutput("clc_mapping_info -s -f " + cas_file)
		if (error):
			sys.exit("[ERROR] - Please add clc_mapping_info to you PATH variable.") 
		p = subprocess.Popen("clc_mapping_info -n " + cas_file , stdout=subprocess.PIPE, bufsize=1, shell=True)
		read_counter = 1
		cas_line_re = re.compile(r"\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+.\d{2})\s+(\d+)\s+(\d+.\d{2})")
		for line in iter(p.stdout.readline, b''):
			match = cas_line_re.search(line)
			if match:
				contig_index = int(match.group(1)) - 1 # -1 because index of contig list starts with zero 
				contig_cov = float(match.group(4))
				contig_id = self.index[contig_index]
				self.addBlobCov(contig_id, lib_name, contig_cov)

	def parseCovFromCovFile(self, lib_name, cov_file):
		'''
		Parse coverage from COV file
		'''
		with open(cov_file) as fh:
			for line in fh:
				cov_line_re = re.compile(r"^(\S+)\t(\d+\.\d+)")
				match = cov_line_re.search(line)
				if match:
					contig_id, contig_cov = match.group(1), float(match.group(2))
					self.addBlobCov(contig_id, lib_name, contig_cov)

	def getTaxForBlobs(self, blast_files, taxdb):
		'''
		- Parses NCBI nodes, names and all BLAST files
		- gets taxonomy for each BLAST hit taxid (getTaxonomy)
		- for each BLAST file it sums up the bitscore by rank_level (e.g. phylum) of hit taxonomy and stores them in blob.tax[blast_lib]
		- then it goes through all blobs again and sets those without hits to 'no-hit'
		'''
		nodes_dict = self.parse_taxdb_nodes(taxdb['nodes'])
		names_dict = self.parse_taxdb_names(taxdb['names'])

		for blast_lib, blast_file in blast_files.items():
			self.blast_libs.append(blast_lib)
			with open(blast_file) as fh:
				for line in fh:
					line_data = line.rstrip("\n").split("\t")
					blast_line_re = re.compile(r"^(\S+)\t(\S+)\t\s*(\S+)") # turns out that blastn output is not tab delimited but tab/(+space) delimited
					match = blast_line_re.search(line)
					if match:
						qseqid, taxid, bitscore = match.group(1), int(match.group(2).split(";")[0]), int(match.group(3)) # if more than one taxid is found ... first one will be used
						if qseqid in self.contigs: 
							taxonomy = {}
							taxonomy = self.getTaxonomy(taxid, taxonomy, nodes_dict, names_dict) # infers taxonomy based on 
							if not blast_lib in self.contigs[qseqid].tax:
								self.contigs[qseqid].tax[blast_lib] = {}
							taxonomic_group = taxonomy[self.rank_level]
							self.contigs[qseqid].tax[blast_lib][taxonomic_group] = self.contigs[qseqid].tax[blast_lib].get(taxonomic_group, 0) + bitscore							
						else:
							sys.exit("[ERROR] - {} in {} does not seem to be part of the assembly.".format(qseqid, blast_file))
					else:
						sys.exit("[ERROR] - BLAST results in {} do not seem to be in the right format. Please run BLAST with the following option : -outfmt '6 qseqid staxids bitscore'".format(blast_file))

			for contig_name in self.contigs:
				if not blast_lib in self.contigs[contig_name].tax:
					self.contigs[contig_name].tax[blast_lib] = {} 
					self.contigs[contig_name].tax[blast_lib]['no-hit'] = 0

	def getTaxonomy(self, taxid, taxonomy, nodes_dict, names_dict):
		'''
		gets target_ranks from nodes/names based on taxid and returns taxonomy
		'''
		if taxid in nodes_dict: 
			rank = nodes_dict[int(taxid)]['rank']
			parent = int(nodes_dict[int(taxid)]['parent'])
			if taxid == 1 or rank == 'superkingdom':
				# finish if taxid == 1 (root) or rank == superkingdom
				taxonomy[rank] = names_dict[int(taxid)] 
				for rank in self.target_ranks:
					# set rank to undef if ranks in target ranks was not populated 
					taxonomy[rank] = taxonomy.get(rank, "undef")
				return taxonomy
			else:
				if rank in self.target_ranks:
					taxonomy[rank] = names_dict[int(taxid)] 
				self.getTaxonomy(parent, taxonomy, nodes_dict, names_dict)
		else:
			# if taxid is not in nodes_dict then print warning and set all ranks to undef
			print "[WARN] - Taxid %s not found in TAXDB. This is probably because you have BLASTed against a newer version of NCBI nt. You should update your TAXDB." % taxid
			taxonomy = {rank: "undef" for rank in self.target_ranks}
		return taxonomy

	def parse_taxdb_names(self, infile):
		print "[STATUS] - Parsing names.dmp from NCBI taxdb" 
		names_dict = {}
		with open(infile) as fh:
			for line in fh:
				fields = line.split("\t")
				if fields[6] == "scientific name":
					names_dict[int(fields[0])] = fields[2]
		return names_dict

	def parse_taxdb_nodes(self, infile):
		print "[STATUS] - Parsing nodes.dmp from NCBI taxdb" 
		nodes_dict = {}
		with open(infile) as fh:
			for line in fh:
				fields = line.split("\t")
				nodes_dict[int(fields[0])]={ 'parent' : fields[2] , 'rank' : fields[4] }
		return nodes_dict

	
	def getConsensusTaxForBlobs(self, taxrule, blast_order):
		''' 
		- Based on taxrule ("A" or "B") and the blast_order (list in order in which blast files where specified) 
		it calculates the consensus taxonomy for each blob 
		- if taxrule == A:
			- it puts all taxonomic groups in a dict with their summed scores as values
			- if a taxonomic group occurs in hits of more than one BLAST file, the highest score is used
		- if taxrule == B:
			- taxonomic groups are put in the dict with their summed scores as values IF they come from the first BLAST file
			- If there was no hit then take the taxonomic groups from the next one  	
		- The highest scoring taxonomic group is selected as consensus taxonomy for each blob
		'''
		for contig_name in self.contigs:
			dict_for_tax_merging = {}
			for blast_lib in blast_order:
				for tax, score in sorted(self.contigs[contig_name].tax[blast_lib].items(), key=lambda x: x[1], reverse=True):
					# loops through tax/score with decreasing score
					if taxrule == 'A':
						if not tax in dict_for_tax_merging:
							dict_for_tax_merging[tax] = score
						else:
							if score > dict_for_tax_merging[tax]:
								dict_for_tax_merging[tax] = score 
					elif taxrule == 'B':
						if blast_lib == blast_order[0]:
							# First blast_lib
							dict_for_tax_merging[tax] = score
						else:
							if len(dict_for_tax_merging) <= 1 and ('no-hit' in dict_for_tax_merging):
								dict_for_tax_merging[tax] = score	

			tax = keyWithMaxVal(dict_for_tax_merging)
			self.contigs[contig_name].tax['tax'] = {}
			self.contigs[contig_name].tax['tax'][tax]=dict_for_tax_merging[tax]
		self.blast_libs.append('tax')

	def getStats(self):
		'''
		Calculates all different kind of stats for each taxonomic group based on each BLAST file and the consensus... 
		'''
		self.stats['count'] = {}
		self.stats['span']= {}
		self.stats['n50']= {}
		self.stats['lengths']= {}
		self.stats['gc']= {}
		self.stats['cov'] = {}
		self.stats['total_count'] = 0
		self.stats['total_span'] = 0
		self.stats['total_n50'] = 0
		self.stats['total_lengths'] = []
		self.stats['total_cov'] = {}
		self.stats['total_gc'] = {'raw' : [], 'mean' : 0.0, 'stdev' : 0.0}
		self.stats['cov_libs'] = []

		for contig_name in self.contigs:
			blob = self.contigs[contig_name]
			self.stats['total_count'] += 1
			self.stats['total_span'] += blob.length
			self.stats['total_lengths'].append(blob.length)
			self.stats['total_gc']['raw'].append(blob.gc)

			for blast_lib in self.blast_libs:
				
				bestTax = keyWithMaxVal(blob.tax[blast_lib])
				if not blast_lib in self.stats['count']:
					self.stats['count'][blast_lib] = {}
					self.stats['span'][blast_lib] = {}
					self.stats['lengths'][blast_lib] = {}
					self.stats['gc'][blast_lib] = {}
					self.stats['cov'][blast_lib] = {}
				self.stats['count'][blast_lib][bestTax] = self.stats['count'][blast_lib].get(bestTax, 0) + 1	
				self.stats['span'][blast_lib][bestTax] = self.stats['span'][blast_lib].get(bestTax, 0) + blob.length

				if not bestTax in self.stats['gc'][blast_lib]:
					self.stats['gc'][blast_lib][bestTax] = {'raw' : [], 'mean' : 0.0, 'stdev' : 0.0}
					self.stats['lengths'][blast_lib][bestTax] = []	
				self.stats['gc'][blast_lib][bestTax]['raw'].append(blob.gc) 
				self.stats['lengths'][blast_lib][bestTax].append(blob.length)
				
				for cov_lib, cov in blob.covs.items():
					if not cov_lib in self.stats['cov'][blast_lib]:
						self.stats['total_cov'][cov_lib] = {'raw' : [], 'mean' : 0.0, 'stdev' : 0.0}
						self.stats['cov'][blast_lib][cov_lib]={}
					if not bestTax in self.stats['cov'][blast_lib][cov_lib]:
						self.stats['cov'][blast_lib][cov_lib][bestTax] = {'raw' : [], 'mean' : 0.0, 'stdev' : 0.0} 
					self.stats['cov'][blast_lib][cov_lib][bestTax]['raw'].append(cov)
			
			for cov_lib, cov in blob.covs.items():
				self.stats['total_cov'][cov_lib]['raw'].append(cov)

		for blast_lib in self.blast_libs:
			# calculate N50
			for tax, list_of_lengths in self.stats['lengths'][blast_lib].items():
				if not blast_lib in self.stats['n50']:
					self.stats['n50'][blast_lib] = {}
				self.stats['n50'][blast_lib][tax] = n50(list_of_lengths)
			self.stats['total_n50'] = n50(self.stats['total_lengths'])

			# calculate total gc mean/stdev
			for tax in self.stats['gc'][blast_lib]:
				self.stats['gc'][blast_lib][tax]['mean'] = "{0:.2f}".format(numpy.mean(self.stats['gc'][blast_lib][tax]['raw']))
				self.stats['gc'][blast_lib][tax]['stdev'] = "{0:.2f}".format(numpy.std(self.stats['gc'][blast_lib][tax]['raw']))

			# calculate total cov mean/stdev
			for cov_lib in self.stats['cov'][blast_lib]:
				self.stats['total_cov'][cov_lib]['mean'] = "{0:.2f}".format(numpy.mean(self.stats['total_cov'][cov_lib]['raw']))
				self.stats['total_cov'][cov_lib]['stdev'] = "{0:.2f}".format(numpy.std(self.stats['total_cov'][cov_lib]['raw']))
				
				# calculate tax-specific cov mean/stdev
				for tax in self.stats['cov'][blast_lib][cov_lib]:
					self.stats['cov'][blast_lib][cov_lib][tax]['mean'] = "{0:.2f}".format(numpy.mean(self.stats['cov'][blast_lib][cov_lib][tax]['raw']))
					self.stats['cov'][blast_lib][cov_lib][tax]['stdev'] = "{0:.2f}".format(numpy.std(self.stats['cov'][blast_lib][cov_lib][tax]['raw']))
		
		self.stats['total_gc']['mean'] = "{0:.2f}".format(numpy.mean(self.stats['total_gc']['raw']))
		self.stats['total_gc']['stdev'] = "{0:.2f}".format(numpy.std(self.stats['total_gc']['raw']))

	def writeOutput(self, version):
		
		'''
		Writes outputfiles:
		- stats.txt which contains tables with the stats calculated through getStats()
		- blobs.txt which contains the blobs
		'''

		header = '# makeblobs.py v{}\n# {} {}\n# {}\n'.format(version, time.strftime("%Y-%m-%d"), time.strftime("%H:%M:%S"), " ".join(sys.argv))
		
		# Writing of stats.txt
		stats_fh = open(self.outfiles['stats'], 'w')
		stats_fh.write(header)
		stats_string = ''
		for blast_lib in self.blast_libs:
			stats_string += "\tTAX:{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:^10}".format(blast_lib, "contigs", "span", "N50", "GC")
			for cov_lib in self.cov_libs:
				stats_string += "\t{:<10}".format(cov_lib)
			stats_string += "\n\t" + (("-" * 10) + "\t") * (5 + len(self.cov_libs)) + "\n"
			for tax, score in sorted(self.stats['span'][blast_lib].items(), key=lambda x: x[1], reverse = True):
				stats_string += "\t{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:>10}".format(tax, self.stats['count'][blast_lib][tax], self.stats['span'][blast_lib][tax], self.stats['n50'][blast_lib][tax],  str(self.stats['gc'][blast_lib][tax]['mean']) + " SD:" + str(self.stats['gc'][blast_lib][tax]['stdev']))
				for cov_lib in self.cov_libs:
					stats_string += "\t{:>10}".format( str(self.stats['cov'][blast_lib][cov_lib][tax]['mean']) + " SD:" + str(self.stats['cov'][blast_lib][cov_lib][tax]['stdev']))
				stats_string += "\n"
			stats_string += "\t{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:>10}".format("Total", self.stats['total_count'], self.stats['total_span'], self.stats['total_n50'],  str(self.stats['total_gc']['mean']) + " SD:" + str(self.stats['total_gc']['stdev']))
			for cov_lib in self.cov_libs:
				stats_string += "\t{:>10}".format( str(self.stats['total_cov'][cov_lib]['mean']) + " SD:" + str(self.stats['total_cov'][cov_lib]['stdev']))
			stats_string += "\n\n"
		stats_fh.write(stats_string)
		stats_fh.close()

		# Writing of blobs.txt
		blobs_fh = open(self.outfiles['blobs'], 'w')
		blobs_fh.write(header)
		blobs_fh.write("# contig_id\tlength\tgc\tcov\ttaxonomy\n")
		blobs_string = ''
		for contig_name in self.index:
			blob = self.contigs[contig_name]
			blobs_string += "{}\t{}\t{:.3f}".format(blob.name, blob.length, blob.gc)
			cov_string, tax_string = '\t','\t'
			for cov_lib, cov in blob.covs.items():
				cov_string += "{}={};".format(cov_lib, cov) 
			blobs_string += cov_string[0:-1]
			for blast_lib, tax in blob.tax.items():
				tax_string += "{}=".format(blast_lib)
				for phylum, score in sorted(tax.items(), key=lambda x: x[1], reverse = True):
					tax_string += "{}:{},".format(phylum, score)
				tax_string = tax_string[0:-1] + ";"
			blobs_string += tax_string[0:-1]
			blobs_string += "\n"
		blobs_fh.write(blobs_string)
		blobs_fh.close()