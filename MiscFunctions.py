#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File   		: MiscFunctions.py
Version 	: 0.1
Author 		: Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Bugs 		: None.
To do 		: 
"""

def keyWithMaxVal(d):
	"""
    http://stackoverflow.com/a/12343826 
    a) create a list of the dict's keys and values; 
	b) return the key with the max value
	"""  
	v=list(d.values())
	k=list(d.keys())
	return k[v.index(max(v))]

def n50(list_of_lengths):
	total_span = 0
	sorted_list_of_lengths=sorted(list_of_lengths, reverse=True)
	for contig_length in sorted_list_of_lengths:
		total_span += contig_length
	teoN50 = total_span/2.0
	running_sum = 0
	N50 = 0
	for contig_length in sorted_list_of_lengths:
		running_sum += contig_length
		if teoN50 <= running_sum:
			N50 = contig_length
			break
	return N50