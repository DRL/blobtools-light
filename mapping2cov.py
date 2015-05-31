#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File   		: mapping2cov.py
Version 	: 0.1
Author 		: Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Descrition	: Takes a BAM/SAM/CAS file and writes a COV file
Bugs 		: ?
To do 		: ?
"""

from __future__ import division
import sys 
import argparse
import os
import BlobCollection

def getInput():
	parser = argparse.ArgumentParser(
		prog='makeblobs.py',
		usage = '%(prog)s -a <ASSEMBLY> -cas <CAS> -blast <BLAST> -taxdb <PATH_TO_TAXDB> -o <OUTPUT> [-h]',
		add_help=True)
	# only ONE assembly file
	parser.add_argument('-a', metavar = 'ASSEMBLY_FASTA', default='', help='Assembly file')
	# multiple mapping files
	parser.add_argument('-bam', metavar = 'BAM_FILE', default=[], nargs='+', help='BAM (mapping) file')
	parser.add_argument('-sam', metavar = 'SAM_FILE', default=[], nargs='+', help='SAM (mapping) file')
	parser.add_argument('-cas', metavar = 'CAS_FILE', default=[], nargs='+', help='CAS (mapping) file')
	# output prefix
	parser.add_argument('-o', metavar = 'OUTPUT_PREFIX', default='', help='Output prefix') 
	# Version number
	parser.add_argument('-v', action='version', version='%(prog)s version 0.1')
	args = parser.parse_args()
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	return InputObject(args)

class InputObject():
	def __init__(self, args):
		self.assembly_type, self.assembly_file = self.getAssembly(args)
		self.mapping_files = self.getDictOfMappingFiles(args)

	def printParameters(self):
		print "\n\tAssembly file :\n\t\t- {} (type = '{}')".format(self.assembly_file, self.assembly_type)
		print "\tMapping file(s) :"
		for lib_name, mapping_file in self.mapping_files.items():
			print "\t\t- {} : {}".format(lib_name, mapping_file)
		if not self.assembly_type == 'unknown':
			if self.exclude_assembly_cov:
				print "\t\t- Exclude assembly coverage : {} ".format(self.assembly_file)
			else:
				print "\t\t- Include assembly coverage : {} ".format(self.assembly_file)
		print


	def getAssembly(self, args):
		assembly = {'unknown' : args.a}
		assembly = {k: v for k, v in assembly.items() if v}	
		try:
			# get assembly type, assembly file and check whether file exists
			[(assembly_type, assembly_file)] = assembly.items()
			if os.path.exists(assembly_file):
				pass
			else:
				sys.exit("[ERROR] : Assembly file {} does not exist.".format(assembly_file))		
		except ValueError, e:
			# throw ValueError if there are too many (>1) elements to unpack from assembly   
			sys.exit("[ERROR] : Please specify an assembly file.")	
		if assembly_type == 'unknown':
			self.exclude_assembly_cov = True
		return assembly_type, assembly_file

	def getDictOfMappingFiles(self, args):
		mappings = {}
		files = {'CAS' : args.cas, 'BAM' : args.bam, 'SAM' : args.sam}
		files = {k: v for k, v in files.items() if v}	
		count = {}
		set_of_files = set()
		for mapping_type, mapping_files in files.items():
			for mapping_file in mapping_files:
				if os.path.exists(mapping_file):
					if mapping_file in set_of_files:
						sys.exit("[ERROR] : Mapping file {} listed more than once.".format(mapping_file))
					set_of_files.add(mapping_file)
					count[mapping_type] = count.get(mapping_type, 0) + 1
					mappings[mapping_type + "_" + str(count[mapping_type])] = mapping_file
				else:
					sys.exit("[ERROR] : Mapping file %s does not exist." %mapping_file)	
		if not mappings and (self.exclude_assembly_cov or self.assembly_type == 'unknown'):
			sys.exit("[ERROR] : Please specify at least one mapping file.")	
		return mappings

def parseInput(parameters):
	'''
	1. Create a BlobCollection Object
	'''
	data = BlobCollection.BlobCollection(RANKS, None)
	'''
	2. Parse contigs of assembly into Blob objects in BlobCollection object
	'''
	data.getBlobsFromAssembly(parameters.assembly_file, parameters.assembly_type, parameters.exclude_assembly_cov)
	'''
	3a. Parse coverage library files into Blob objects in BlobCollection object
	'''
	data.getCovForBlobs(parameters.mapping_files)
	'''
	3b. Print COVs to files if mappings are BAM/SAM/CAS
	'''
	data.printCOVToFiles(parameters.mapping_files)
	'''
	4. Return BlobCollection object
	'''
	return data

if __name__ == "__main__":
	
	__version__ = 0.1

	RANKS = [] # LEGACY

	# Put all parameters of the run into InputObject
	parameters = getInput()
	
	# Print InputObject so that user gets feedback about whats happening
	parameters.printParameters()
	
	# Parse files according to InputObject into BlobCollection object 
	data = parseInput(parameters)
