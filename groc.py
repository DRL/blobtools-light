#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File    	: groc.py
Author  	: Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Version 	: 0.8
Description : groc.py ("Get-Reads-Of-Casfile") takes a CLC mapping file in *.cas format and filters reads based on contig lists
Requirements: *.fa and .1, .2 ... mapping with '-q -i' 
Beware of 	: Non-unique read names (checking for duplicates would make the whole thing much much slower)
To do 		: Check that exclude/include file contains only one column
"""

from __future__ import division
import os, re, sys, argparse, subprocess, itertools, commands
from itertools import izip

def get_input():
	'''Gets, checks and returns input'''
	parser = argparse.ArgumentParser(
		prog='groc.py',
		usage = '%(prog)s -c -e [-dry -h]',
		add_help=True)
	parser.add_argument('-c', metavar = 'cas', default='', help='*.cas file which will be parsed')
	parser.add_argument('-i', metavar = 'include', default='', help='NOT IMPLEMENTED YET.List of contigs which determines which reads should be INCLUDED (Incompatible with -e)')
	parser.add_argument('-e', metavar = 'exclude', default='', help='List of contigs which determines which reads should be EXCLUDED (Incompatible with -i)') 
	parser.add_argument('-o', metavar = 'out', default='', help='Add output prefix to filtered read files') 
	parser.add_argument('-u', action='store_true' , help='Exclude unmapped reads.')
	parser.add_argument('-dry', action='store_true' , help='Set flag for optional dry-run') # dry is for only running until end of cas parsing to output stats and then prompt to continue 


	args = parser.parse_args()

	cas_file, i_file, e_file, dry, exclude_unmapped, out_prefix = args.c, args.i, args.e, args.dry, args.u, args.o

	if not cas_file:
		sys.exit("ERROR: Please specify a CLC mapping file")
	if not os.path.isfile(cas_file):
		sys.exit("ERROR: " + cas_file + " is not a file")
	if i_file and e_file:
		sys.exit("ERROR: Please specify either -i OR -e")
	if not i_file and not e_file:
		sys.exit("ERROR: Please specify -i or -e")
	if i_file and not os.path.isfile(i_file):
		sys.exit("ERROR: " + i_file + " is not a file")
	#for e_file in e_files:
	if e_file and not os.path.isfile(e_file):
		sys.exit("ERROR: " + e_file + " is not a file")

	return cas_file, i_file, e_file, dry, exclude_unmapped, out_prefix

def get_read_info(cas_file):

	"Uses clc_mapping_info. Returns read files and number of reads."
	read_file_type = ''
	reads_re = re.compile(r" -q -i (\S+) (\S+)")
	assembly_re = re.compile(r"-d (\S+)")

	number_of_reads_re = re.compile(r"\s+Reads\s+(\d+)")
	error, message = commands.getstatusoutput("clc_mapping_info -s -f " + cas_file)
	if (error):
		sys.exit("ERROR: Please load the CLC module ('module load clc')") 
	mapping_info = subprocess.check_output("clc_mapping_info -s -f " + cas_file, stderr=subprocess.STDOUT, shell=True)
	# if not re.match(q_i_re, str(mapping_info)):
	#	sys.exit("ERROR: The mapping was not carried out with the -i parameter")
	read_files = re.search(reads_re, str(mapping_info)).group(1,2)
	assembly_file = re.search(assembly_re, str(mapping_info)).group(1)
	number_of_reads = re.search(number_of_reads_re, str(mapping_info)).group(1)
	print "The *.cas file indicates " + str(int(int(number_of_reads)/2)) + " read pairs"
	print

	if (read_files[0][-2:] == "fa" and read_files[1][-2:] == "fa") or (read_files[0][-5:] == "fasta" and read_files[1][-5:] == "fasta"):
		read_file_type = "fa" 
	elif (read_files[0][-2:] == "fq" and read_files[1][-2:] == "fq") or (read_files[0][-5:] == "fastq" and read_files[1][-5:] == "fastq"):
		read_file_type = "fq"
	else:
		sys.exit("ERROR: Both read have to be *.fa/*.fasta or *.fq/*.fastq.")

	return number_of_reads, read_files, read_file_type, assembly_file

def parse_contig_list(contig_file):

	'''Parses contig list file and returns the set "contigs".'''
	contigs = set()
	with open(contig_file) as fh:
		for line in fh:
			if line.startswith("#"):
				pass
			else:
				contig = line.lstrip(">").rstrip("\n")
				contigs.add(contig)
	if len(contigs) == 0:
		sys.exit("ERROR: " + contig_file + " is empty. Nothing to exclude.")
	return contigs

def parse_assembly_file(assembly_file):

	'''Parses assembly file and returns the names of the contigs.'''

	assembly_contigs = []
	with open(assembly_file) as fh:
		for line in fh:
			if line.startswith(">"):
				assembly_contigs.append(line.lstrip(">").rstrip("\n"))
	return assembly_contigs

def parse_cas(cas_file, assembly_contigs, contigs_to_exclude, number_of_reads):

	'''Uses clc_mapping_table. Parses *.cas file and unless both reads map to a contig in the set "contigs" they are added to the set "reads". 
	Returns the set "reads".'''

	print "Started parsing *.cas file : " + cas_file

	included_readset = set()
	excluded_readset = set()
	unmapped_readset = set() #

	i = 1
	line_counter = 1
	pair_counter = 0
	both_excluded = 0
	both_unmapped = 0
	both_included = 0
	one_excluded_one_unmapped = 0
	one_excluded_one_included = 0
	one_included_one_unmapped = 0

	excluded_reads = []
	included_reads = []
	unmapped_reads = [] #

	p = subprocess.Popen("clc_mapping_table -n -p " + cas_file, stdout=subprocess.PIPE, bufsize=1, shell=True)
	
	for line in iter(p.stdout.readline, b''):

		if line_counter % 5000 == 0:		
			sys.stdout.write('\r')
			progress = int(line_counter)/int(number_of_reads)
			print "Progress:\t" + format(float(progress),'.2%'),
			sys.stdout.flush()

		column = line.rstrip("\n").split()

		contig_index = column[5]

		if contig_index == '-1':
			unmapped_reads.append(column)
		elif assembly_contigs[int(contig_index)] in contigs_to_exclude:  
			excluded_reads.append(column)
		elif (assembly_contigs[int(contig_index)] not in contigs_to_exclude):
			included_reads.append(column)
		else:
			sys.exit("ERROR: " + column)
			
		
		if not i % 2 == 0: # First read
			i += 1
		else:  		   # Second read
			i = 1
			pair_counter += 1
			if len(excluded_reads) == 2: # both reads excluded 
				both_excluded += 1
				column_1 = excluded_reads[0]
				column_2 = excluded_reads[1]
				excluded_readset.add(column_1[1])
				excluded_readset.add(column_2[1])
			elif len(unmapped_reads) == 2: # both reads do not map #
				both_unmapped += 1 #
				column_1 = unmapped_reads[0] #
				column_2 = unmapped_reads[1] #
				unmapped_readset.add(column_1[1]) #
				unmapped_readset.add(column_2[1]) #
				if not (exclude_unmapped):
					included_readset.add(column_1[1]) #
					included_readset.add(column_2[1]) #
			elif len(included_reads) == 2: # both reads included
				both_included += 1
				column_1 = included_reads[0]
				column_2 = included_reads[1]
				included_readset.add(column_1[1])
				included_readset.add(column_2[1])
			elif (len(excluded_reads) == 1) and (len(unmapped_reads) == 1): # one read excluded, one unmapped
				one_excluded_one_unmapped += 1
				column_1 = excluded_reads[0]
				column_2 = unmapped_reads[0]
				excluded_readset.add(column_1[1])
				excluded_readset.add(column_2[1])
			elif (len(excluded_reads) == 1) and (len(unmapped_reads) == 0): # one read excluded, one included
				one_excluded_one_included += 1
				column_1 = excluded_reads[0]
				column_2 = included_reads[0]
				included_readset.add(column_1[1])
				included_readset.add(column_2[1])
			elif (len(unmapped_reads) == 1) and (len(included_reads) == 1): # one read included, one unmapped
				one_included_one_unmapped += 1
				column_1 = included_reads[0]
				column_2 = unmapped_reads[0]
				included_readset.add(column_1[1])
				included_readset.add(column_2[1])
			else:
				print included_reads
				print excluded_reads
				print unmapped_reads
				sys.exit("Error : something broke ...")

			included_reads = []
   			excluded_reads = []
   			unmapped_reads = [] #
   		line_counter += 1

   	number_of_pairs = int(int(number_of_reads)/2)

   
	print "Finished parsing *.cas file"
   	print
   	print "    Total pairs     : \t\t" + format(pair_counter,',d').rjust(len(str(pair_counter)) + 2, ' ') + "\t(" + format(float(pair_counter/number_of_pairs),'.2%').rjust(6, ' ') + ")"
   	print "-------------------- --------------------------------------------------- "
   	print "Excluded | Excluded : \t\t" + format(both_excluded,',d').rjust(len(str(pair_counter)) + 2, ' ') + "\t(" + format(float(both_excluded/number_of_pairs),'.2%').rjust(6, ' ') + ")\t[EXCLUDED]"
   	print "Excluded | Unmapped : \t\t" + format(one_excluded_one_unmapped,',d').rjust(len(str(pair_counter)) + 2, ' ') + "\t(" + format(float(one_excluded_one_unmapped/number_of_pairs),'.2%').rjust(6, ' ') + ")\t[EXCLUDED]"
   	print "Unmapped | Unmapped : \t\t" + format(both_unmapped,',d').rjust(len(str(pair_counter)) + 2, ' ') + "\t(" + format(float(both_unmapped/number_of_pairs),'.2%').rjust(6, ' ') + ")\t",
   	if (exclude_unmapped):
   		print "[UNMAPPED]"
   	else:
   		print "[UNMAPPED] + [INCLUDED]"
   	print "Included | Included : \t\t" + format(both_included,',d').rjust(len(str(pair_counter)) + 2, ' ') + "\t(" + format(float(both_included/number_of_pairs),'.2%').rjust(6, ' ') + ")\t[INCLUDED]"
   	print "Included | Unmapped : \t\t" + format(one_included_one_unmapped,',d').rjust(len(str(pair_counter)) + 2, ' ') + "\t(" + format(float(one_included_one_unmapped/number_of_pairs),'.2%').rjust(6, ' ') + ")\t[INCLUDED]"
   	print "Included | Excluded : \t\t" + format(one_excluded_one_included,',d').rjust(len(str(pair_counter)) + 2, ' ') + "\t(" + format(float(one_excluded_one_included/number_of_pairs),'.2%').rjust(6, ' ') + ")\t[INCLUDED]"
   	print "-------------------- --------------------------------------------------- "
   	print "Excluded read pairs : \t\t" + format(int(len(excluded_readset)/2),',d').rjust(len(str(pair_counter)) + 2, ' ') + "\t(" + format(float(len(excluded_readset)/int(number_of_reads)),'.2%').rjust(6, ' ') + ")"
   	print "Unmapped read pairs : \t\t" + format(int(len(unmapped_readset)/2),',d').rjust(len(str(pair_counter)) + 2, ' ') + "\t(" + format(float(len(unmapped_readset)/int(number_of_reads)),'.2%').rjust(6, ' ') + ")"
   	print "Included read pairs : \t\t" + format(int(len(included_readset)/2),',d').rjust(len(str(pair_counter)) + 2, ' ') + "\t(" + format(float(len(included_readset)/int(number_of_reads)),'.2%').rjust(6, ' ') + ")"
   	print 

   	user_input = ''

   	if (dry):
   		while (user_input.strip() not in ['Y','y','N','n'] ):
   			user_input = raw_input('Do you want to filter the read files [y/n] : ')
   			user_input.split()
   		if user_input.strip() in ['N','n']:
   			sys.exit("\nTerminated")

   	if (len(excluded_readset) or len(unmapped_readset)):
   		return included_readset, excluded_readset, unmapped_readset 
   	else:
   		sys.exit("No reads to filter ... ")
   	
def parse_read_files(read_files, included_readset, excluded_readset, unmapped_readset, read_file_type, number_of_reads, out_prefix):
	
	'''Parses both read files and if they are contained in the set "reads" they are written to file.'''

	print "Started filtering read files : " + read_files[0] + " and " + read_files[1]
	
	included_out_1 = open(read_files[0] + ".included." + out_prefix + "." + read_file_type, 'w')
	included_out_2 = open(read_files[1] + ".included." + out_prefix + "." + read_file_type, 'w')
	excluded_out_1 = open(read_files[0] + ".excluded." + out_prefix + "." + read_file_type, 'w')
	excluded_out_2 = open(read_files[1] + ".excluded." + out_prefix + "." + read_file_type, 'w')
	unmapped_out_1 = open(read_files[0] + ".unmapped." + out_prefix + "." + read_file_type, 'w') #
	unmapped_out_2 = open(read_files[1] + ".unmapped." + out_prefix + "." + read_file_type, 'w') #

	this_seq_good = 0
	this_seq_unmapped = 0 #
	lines_parsed = 0
	reads_filtered = 0
	chunk_size = 0

	if (read_file_type == "fa"):
		chunk_size = 2
	else:
		chunk_size = 4

	with open(read_files[0]) as infile_A, open(read_files[1]) as infile_B: 

		line_number = chunk_size
		unmapped_span = 0
		included_span = 0
		excluded_span = 0
		total_span = 0
		number_of_unmapped_reads = 0
		number_of_included_reads = 0
		number_of_excluded_reads = 0
		for line_A, line_B in izip(infile_A, infile_B):

			line_A = line_A.rstrip("\n")
			line_B = line_B.rstrip("\n")

			if lines_parsed % 5000 == 0:
				sys.stdout.write('\r')
				progress = int(reads_filtered)/int(number_of_reads)
				print "Progress:\t" + format(float(progress),'.2%'),
				sys.stdout.flush()

			if line_number == chunk_size:
				
				lines_parsed += 1

				# Here is where the logic happens ...

				if (line_A[1:] in unmapped_readset) and (line_B[1:] in unmapped_readset):
					unmapped_out_1.write(line_A + "\n")
					unmapped_out_2.write(line_B + "\n")
					this_seq_unmapped = 1
					number_of_unmapped_reads += 2
					if not (exclude_unmapped):
						included_out_1.write(line_A + "\n")
						included_out_2.write(line_B + "\n")
						number_of_included_reads += 2
						this_seq_good = 1
						reads_filtered += 2

				elif (line_A[1:] in included_readset) and (line_B[1:] in included_readset):
					included_out_1.write(line_A + "\n")
					included_out_2.write(line_B + "\n")						
					this_seq_good = 1
					number_of_included_reads += 2
					reads_filtered += 2

				elif (line_A[1:] in excluded_readset) and (line_B[1:] in excluded_readset):
					excluded_out_1.write(line_A + "\n")
					excluded_out_2.write(line_B + "\n")
					number_of_excluded_reads += 2
					this_seq_good = 0

				else:
					sys.exit("ERROR : " + line_A + " and/or " + line_B + " were not encountered when reading the *.cas file")
				line_number -= 1

			elif (line_number < chunk_size) and (line_number):
				lines_parsed += 1
				if (this_seq_unmapped):
					if (line_number == (chunk_size - 1)): # nucleotides
						unmapped_span += len(line_A) + len(line_B)
						total_span += len(line_A) + len(line_B)
					unmapped_out_1.write(line_A + "\n")
					unmapped_out_2.write(line_B + "\n")
					if not (exclude_unmapped) and (this_seq_good):
						if (line_number == (chunk_size - 1)): # nucleotides
							included_span += len(line_A) + len(line_B)
						included_out_1.write(line_A + "\n")
						included_out_2.write(line_B + "\n")

				elif (this_seq_good):
					if (line_number == (chunk_size - 1)): # nucleotides
						included_span += len(line_A) + len(line_B)
						total_span += len(line_A) + len(line_B)
					included_out_1.write(line_A + "\n")
					included_out_2.write(line_B + "\n")
				
				else:
					if (line_number == (chunk_size - 1)): # nucleotides
						excluded_span += len(line_A) + len(line_B)
						total_span += len(line_A) + len(line_B)
					excluded_out_1.write(line_A + "\n")
					excluded_out_2.write(line_B + "\n")
				line_number -= 1
			else:
				pass

			if line_number == 0:
				this_seq_good = 0
				this_seq_unmapped = 0 #
				line_number = chunk_size
	
	sys.stdout.write('\r')
	progress = int(lines_parsed)/int(number_of_reads)
	print "Progress:\t100.00%"
	print "Finished filtering read files"
	print
	print "    Total span      :\t" + format(total_span,',d').rjust(len(str(total_span)) + 2, ' ') + "nt\t(" + format((total_span/total_span),'.2%').rjust(6, ' ') + ")\tMean : " + ("%.2f" % (total_span/int(number_of_reads))).rjust(6, ' ')+ "nt"
	print "-------------------- --------------------------------------------------- "
	print "   Included span    :\t" + format(included_span,',d').rjust(len(str(total_span)) + 2, ' ') + "nt\t(" + format((included_span/total_span),'.2%').rjust(6, ' ') + ")\tMean : " + ("%.2f" % (included_span/int(number_of_included_reads))).rjust(6, ' ')+ "nt"
	print "   Unmapped span    :\t" + format(unmapped_span,',d').rjust(len(str(total_span)) + 2, ' ') + "nt\t(" + format((unmapped_span/total_span),'.2%').rjust(6, ' ') + ")\tMean : " + ("%.2f" % (unmapped_span/int(number_of_unmapped_reads))).rjust(6, ' ')+ "nt"
	print "   Excluded span    :\t" + format(excluded_span,',d').rjust(len(str(total_span)) + 2, ' ') + "nt\t(" + format((excluded_span/total_span),'.2%').rjust(6, ' ') + ")\tMean : " + ("%.2f" % (excluded_span/int(number_of_excluded_reads))).rjust(6, ' ')+ "nt"
	print "-------------------- --------------------------------------------------- "
	print
	print "Cleaning up ... \n"

	included_out_1.close()
	included_out_2.close()
	excluded_out_1.close()
	excluded_out_2.close()
	unmapped_out_1.close()
	unmapped_out_2.close()

	print "Done."

if __name__ == "__main__":

	cas_file, i_file, e_file, dry, exclude_unmapped, out_prefix = get_input()

 	read_files = []
 	number_of_reads, read_files, read_file_type, assembly_file = get_read_info(cas_file)
 	
 	contigs_to_include = set()
 	if i_file:
 		contigs_to_include = parse_contig_list(i_file)
 		sys.exit("ERROR: This functions is not yet implemented.")

 	contigs_to_exclude = set()
 	if e_file:
 		contigs_to_exclude = parse_contig_list(e_file)

 	assembly_contigs = parse_assembly_file(assembly_file)
	included_readset, excluded_readset, unmapped_readset = parse_cas(cas_file, assembly_contigs, contigs_to_exclude, number_of_reads)
	parse_read_files(read_files, included_readset, excluded_readset, unmapped_readset, read_file_type, number_of_reads, out_prefix)