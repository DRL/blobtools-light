#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File   	: plotblobs.py
Version : 0.1
Author 	: Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Bugs 	: ?
To do 	: 	- Add read proportion histogram 
			- legend placement at bottom
			- total reads mapped per phylum
		
"""

# # # # # 
# MODULES										
# # # # # 

from __future__ import division

import numpy as np
import math as math
import matplotlib as mat
import matplotlib.pyplot as plt
import sys
import argparse
import os

from matplotlib import cm
from matplotlib.ticker import NullFormatter
from matplotlib.lines import Line2D

from MiscFunctions import n50

mat.rcParams.update({'font.size': 30})
mat.rcParams['xtick.major.pad']='8'
mat.rcParams['ytick.major.pad']='8'
mat.rcParams['lines.antialiased']=True

def set_canvas():
	left, width = 0.1, 0.60
	bottom, height = 0.1, 0.60
	bottom_h = left_h = left+width+0.02
	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]
	rect_legend = [left_h, bottom_h, 0.2, 0.2]
	return rect_scatter, rect_histx, rect_histy, rect_legend

def set_format_scatterplot(axScatter):
	axScatter.set_xlabel("GC proportion", fontsize=35)
	axScatter.set_ylabel("Coverage", fontsize=35)
	axScatter.grid(True, which="major", lw=2., color=white, linestyle='-') 
	axScatter.set_axisbelow(True)
	axScatter.set_xlim( (0, 1) )
	axScatter.set_ylim( (0.01, max_cov+1000) ) # This sets the max-Coverage so that all libraries + sum are at the same scale
	axScatter.xaxis.labelpad = 20
	axScatter.xaxis.labelpad = 20
	return axScatter

def set_format_hist_x(axHistx, axScatter):
	axHistx.set_xlim( axScatter.get_xlim() )
	axHistx.grid(True, which="major", lw=2., color= white, linestyle='-')
	axHistx.xaxis.set_major_formatter(nullfmt) # no labels since redundant
	axHistx.set_axisbelow(True)
	axHistx.yaxis.labelpad = 20
	return axHistx

def set_format_hist_y(axHisty, axScatter):
	axHisty.set_yscale('log')
	axHisty.yaxis.set_major_formatter(nullfmt) # no labels since redundant
	axHisty.set_ylim( axScatter.get_ylim() )
	axHisty.grid(True, which="major", lw=2., color= white, linestyle='-')
	axHisty.set_axisbelow(True)
	axHisty.xaxis.labelpad = 20
	return axHisty

def plot_ref_legend(axScatter):
	s = 15
	# markersize in scatter is in "points^2", markersize in Line2D is in "points" ... that's why we need math.sqrt()
	ref_1 = (Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=math.sqrt(1000/15),  markerfacecolor=grey))
	ref_2 = (Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=math.sqrt(5000/15), markerfacecolor=grey))
	ref_3 = (Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=math.sqrt(10000/15), markerfacecolor=grey))
	axScatter.legend([ref_1,ref_2,ref_3], ["1,000nt", "5,000nt", "10,000nt"], numpoints=1, loc = 4, fontsize=fontsize)

def plot(data, cov_data, outfile):
	""" Plotting function which gets masked data and plots to outfile"""

	rect_scatter, rect_histx, rect_histy, rect_legend = set_canvas()
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Setting up plots and axes
	plt.figure(1, figsize=(35,35), dpi=400)
	axScatter = plt.axes(rect_scatter, axisbg=background_grey, yscale = 'log')
	axScatter = set_format_scatterplot(axScatter)
	axHistx = plt.axes(rect_histx, axisbg=background_grey)
	axHistx = set_format_hist_x(axHistx, axScatter)
	axHisty = plt.axes(rect_histy, axisbg=background_grey)
	axHisty = set_format_hist_y(axHisty, axScatter)
	axScatter.yaxis.get_major_ticks()[0].label1.set_visible(False)
	axScatter.yaxis.get_major_ticks()[1].label1.set_visible(False)
	#plt.suptitle(out_file, fontsize=25, verticalalignment='bottom')
	
	axLegend = plt.axes(rect_legend, axisbg=white)
	axLegend.xaxis.set_major_locator(plt.NullLocator())
	axLegend.xaxis.set_major_formatter(nullfmt)
	axLegend.yaxis.set_major_locator(plt.NullLocator())
	axLegend.yaxis.set_major_formatter(nullfmt)
	#
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Setting bins for histograms
	top_bins = np.arange(0, 1, 0.01)
	right_bins = np.logspace(-2, (int(math.log(max_cov)) + 1), 200, base=10.0)

	# empty handles for big legend
	legend_handles = []
	legend_labels = []

	# change file name if span (should be in input parsing function)
	if hist_span:
		outfile += ".hist_span"
	else:
		outfile += ".hist_count"

	# counter necessary for multiplot so that PNGs are in order when sorted by name
	i = 0

	# initiate variables for plotting
	s, lw, alpha, color = 0, 0, 0, ''
	
	# Maybe make an STDOUT printing func?
	print "[STATUS] Plotting : " + outfile

	# for each phylum ... they are ordered
	for tax in tax_list:

		i += 1

		# get indices for those rows in data where the phylum == tax
		index_for_tax = np.where(data[:,3].astype(str) == tax, True, False)
		# count of contigs ... total number of contigs comes from previous step?
		number_of_contigs_for_tax = np.sum(index_for_tax)
		
		# uses number_of_contigs for checking whether plotting should be carried out ... maybe there is a better place for this ...
		if number_of_contigs_for_tax == 0:
			pass
		else:
			# sums span for phylum in mb and not ... do we need both?
			span_of_contigs_for_tax = np.sum(data[index_for_tax][:,1].astype(int))
			span_of_contigs_for_tax_in_mb = span_of_contigs_for_tax/1000000

			# create np_arrays for length, gc and cov for all contigs in phylum 
			len_array = data[index_for_tax][:,1].astype(int)
			gc_array = data[index_for_tax][:,2].astype(float)
			cov_array = cov_data[index_for_tax].astype(float)
			# generates label ... this should be turned into a table ...
			label = tax + " (" + "{:,}".format(number_of_contigs_for_tax) + "; " + "%.2f" % round(span_of_contigs_for_tax_in_mb,2) + "MB; " + "{:,}".format(n50(len_array)) + "nt)"

			# another status message
			print "\t" + label 
			s_array = []
			# ignore contig length ... maybe do this in input and set these params for plotting there ...
			if (ignore_contig_len):
				if tax == 'no-hit':
					s, lw, alpha, color = 15, 0.5, 0.5, grey
				else:
					s, lw, alpha, color = 65, 0.5, 1, color_dict[tax]
				s_array = [s for contig_length in len_array]
			else:
				if tax == 'no-hit':
					s, lw, alpha, color = 15, 0.5, 0.5, grey
				else:
					s, lw, alpha, color = 15, 0.5, 1, color_dict[tax]
				# these are the sizes for plotting with contig sizes
				s_array = [contig_length/s for contig_length in len_array]
			
			# making copies of gc/cov_array
			gc_hist_array = gc_array
			cov_hist_array = cov_array

			#######
			# if hist span ... 
			# 	make a new array ...
			# 	add to the array : (gc * len/1000) - 1
			# substitute old array with new array

			# set histogram labels depending on type ... can be set before ... 
			weights_array = len_array/1000 
			if (hist_span):
				axHistx.set_ylabel("Span (kb)")
				axHisty.set_xlabel("Span (kb)", rotation='horizontal')
			else:
				axHistx.set_ylabel("Count")
				axHisty.set_xlabel("Count", rotation='horizontal')
		
			# this should be set before ... or after? ... but only once 
			for xtick in axHisty.get_xticklabels(): # rotate text for ticks in cov histogram 
				xtick.set_rotation(270)

			# add text to legend ... label was build before ... could be a function
			legend_handles.append(Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=24, markerfacecolor=color))
			legend_labels.append(label)
			
			if (number_of_contigs_for_tax):
				if (hist_span):
					axHistx.hist(gc_hist_array, weights=weights_array , color = color, bins = top_bins, histtype='step', lw = 3)
					axHisty.hist(cov_hist_array, weights=weights_array , color = color, bins = right_bins, histtype='step', orientation='horizontal', lw = 3)
				else:			
					axHistx.hist(gc_hist_array, color = color, bins = top_bins, histtype='step', lw = 3)
					axHisty.hist(cov_hist_array , color = color, bins = right_bins, histtype='step', orientation='horizontal', lw = 3)
		
			axScatter.scatter(gc_array, cov_array, color = color, s = s_array, lw = lw, alpha=alpha, edgecolor=black, label=label)
		
			axLegend.axis('off')

			if (multi_plot): # MULTI-PLOT!!!
				axLegend.legend(legend_handles, legend_labels, loc=6, numpoints=1, fontsize=fontsize, frameon=True)
				plot_ref_legend(axScatter)
				#plt.savefig(outfile + "." + str(i) + "_"+ tax.replace("/","") + "." + fig_format, format=fig_format)
				plt.savefig(outfile + "." + str(i) + "_"+ tax.replace("/","") + "." + fig_format, format=fig_format)
	
	if ignore_contig_len:
		pass
	else: # print scale-legend
		plot_ref_legend(axScatter)

	axLegend.legend(legend_handles, legend_labels, numpoints=1, fontsize=fontsize, frameon=True, loc=6 )		
	sys.stdout.write("Saving file " + outfile)
	plt.savefig(outfile + "." + fig_format, format=fig_format)
	plt.close()
	print " [Done]\n" 

def getInput():

	parser = argparse.ArgumentParser(
		prog='tagc_plot.py',
		usage = '%(prog)s infile [-p] [-f] [-t] [-e] [-n] [-o] [-l] [-c] [-s] [-h]',
		add_help=True)
	parser.add_argument('i', metavar = 'infile', help='Input file (blobplot.txt)')
	parser.add_argument('-p', metavar = 'max_taxa_plot', default=7, type = int, help='Maximum number of phyla to plot (Default = 7)')
	#parser.add_argument('-t', metavar = 'tax_level', default=2, type = int, help='Taxonomic level on which to plot.  Species = 0, Order = 1, Phylum = 2, Superkingdom = 3 (Default = 2)')
	#parser.add_argument('-e', metavar = 'eval_cutoffs' , default=[1.0], type = float, nargs='+', help='Set maximal e-value(s) (Default = 1.0)') 
	parser.add_argument('-c', metavar = 'len_cutoffs' , default=[100], type = int, nargs='+', help='Set minium contig length(s) (Default = 100)') 
	parser.add_argument('-s', action='store_true' , help='Ignore contig length for plotting.') 
	parser.add_argument('-n', action='store_true', help='Hides "no-hit" contigs') 
	parser.add_argument('-o', metavar ='out_prefix', default='' , help='Set output file prefix.') 
	parser.add_argument('-m', action='store_true' , help='Multi-plot. Print PNG after each tax-addition.') 
	parser.add_argument('-sort', action='store_false' , help='Sort by number of contigs per phylum (Default: Sort by span by tax)') 
	parser.add_argument('-span', action='store_true' , help='Make histograms based on assembled span by tax (Expensive!)') 
	parser.add_argument('-v', action='version', version='%(prog)s version 0.1')
	args = parser.parse_args()
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	out_prefix, sort_by_span, multi_plot = args.o, args.sort, args.m

	len_cutoffs, ignore_contig_len, hide_not_annotated, hist_span = args.c, args.s, args.n, args.span

	infile = args.i

	if not os.path.exists(infile):
		parser.exit("[ERROR] : Input file {} does not exist!".format(infile))

	max_taxa_plot = args.p
	
	if (max_taxa_plot < 1):
		parser.exit("[ERROR] : 'max_taxa_plot' must be a positive integer!")
	
	return infile, out_prefix, max_taxa_plot, len_cutoffs, ignore_contig_len, sort_by_span, multi_plot, hist_span  

def parseInfile(data):
	
	contig_data_list = []
	cov_data_dict = {}
	tax_dict = {}
	with open(infile) as fh:
		for line in fh:
			if line.startswith('#'):
				pass
			else:
				line_data = line.rstrip("\n").split("\t")
				contig_id = line_data[0]
				length = int(line_data[1])
				gc = float(line_data[2])
				cov_dict = dict(string.split('=') for string in line_data[3].split(";"))
				cov_dict = {k: (float(v) if float(v) >= 0.1 else 0.1) for k, v in cov_dict.items()} # set coverages below 0.1 to 0.1
				blast_dict = dict(string.split('=') for string in line_data[4].split(";"))
				blast_dict = {k: v.split(":")[0] for k, v in blast_dict.items()} # removing bitscores, since not needed anymore
				tax = blast_dict['tax']
				
				contig_data_list.append([(contig_id), (length), (gc), (tax)])
				for cov_lib, cov in cov_dict.items():
					if not cov_lib in cov_data_dict:
						cov_data_dict[cov_lib] = []
					cov_data_dict[cov_lib].append(cov)

				if not tax in tax_dict:
					tax_dict[tax] = {}
				tax_dict[tax]['count'] = tax_dict[tax].get('count', 0) + 1
				tax_dict[tax]['span'] = tax_dict[tax].get('span', 0) + int(length)

	contig_data_array = np.array(contig_data_list)
	cov_data_array = {}
	for cov_lib, cov_data in cov_data_dict.items():
		cov_array = np.array(cov_data, dtype=np.float)
		cov_data_array[cov_lib]=cov_array

	return contig_data_array, tax_dict, cov_data_array

def getMasks(data, len_cutoffs):
	""" Returns dict with parameters as keys and mask-array as values. """
	mask_dict = {}
	for len_cutoff in len_cutoffs:
		key = str(len_cutoff) #+ "_" + str(eval_cutoff)
		mask = np.where(data[:,1].astype(int) >= len_cutoff) #& (data[:,4].astype(float) <= eval_cutoff))	
		mask_dict[key]=mask
	return mask_dict

def getSortedTax(tax_dict, sort_by_span, max_taxa_plot):
	""" Returns list of tax of size max_taxa_plot sorted by span or count. """
	tax_list = []
	if (sort_by_span):
		tax_list = sorted(tax_dict, key = lambda x : tax_dict[x]['span'], reverse=True)
	else:
		tax_list = sorted(tax_dict, key = lambda x : tax_dict[x]['count'], reverse=True)

	tax_list = tax_list[0:max_taxa_plot] # only return those that will be plotted	
	return tax_list

def getColorDict(tax_list, colormap):
	""" Returns colour dict, Not annotated is always grey. """
	colors = cm.get_cmap(name=colormap)
	color_index = 1
	color_dict = {}
	for tax in tax_list:
		if tax == "no-hit":
			color_dict[tax] = grey	
		else:
			color_dict[tax] = mat.colors.rgb2hex(colors(1.0 * (color_index/len(tax_list))))
			color_index += 1
	return color_dict

def getMinMaxCov(cov_dict):
	max_cov, min_cov = 100.00, 100.00
	for lib in cov_dict:
		lib_max_cov, lib_min_cov = np.amax(cov_dict[lib].astype(float)), np.amin(cov_dict[lib].astype(float))
		if lib_max_cov > max_cov:
			max_cov = lib_max_cov
		if lib_min_cov < min_cov:
			min_cov = lib_min_cov
	print "[STATUS] - Max.cov = " + str(max_cov) + " / Min.cov = " + str(min_cov)
	return max_cov, min_cov

if __name__ == "__main__":
	fontsize = 24
	colormap = "Set2" # "Paired"
	black, grey, background_grey, white = '#262626', '#d3d3d3', '#F0F0F5', '#ffffff'
	fig_format = 'png'
	nullfmt = NullFormatter()         # no labels on axes

	infile, out_prefix, max_taxa_plot, len_cutoffs, ignore_contig_len, sort_by_span, multi_plot, hist_span = getInput()

	data, tax_dict, cov_dict = parseInfile(infile) # data is a numpy array, tax_dict is a dict, cov_dict is a dict of numpy arrays

	mask_dict = getMasks(data, len_cutoffs) # allows filtering of blobs by length 

	tax_list = getSortedTax(tax_dict, sort_by_span, max_taxa_plot)
	
	color_dict = getColorDict(tax_list, colormap)

	max_cov, min_cov = getMinMaxCov(cov_dict)

	for lib in sorted(cov_dict):
	# sanitation of names occurs in cov_dict 
		for key in mask_dict:
			mask = mask_dict[key]
			cov = cov_dict[lib]
			#print lib
			#print cov[mask]
			if out_prefix: 
				outfile = out_prefix + "." + lib + "." + key 
			else: 
				outfile =  infile + "." + lib + "." + key 
			plot(data[mask], cov[mask], outfile)
		#