#!/usr/bin/env python
import os
import sys
import re
import time
import argparse
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as mpl

# absolute path to this script
SIM_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1])
sys.path.append(SIM_PATH+'/py/')

from sysFunc import *
from samFunc import *
from strFunc import *
from graphFunc import *

from msa import msa_consensus
from build_report import gen_plots

##### EXECUTABLES
####BLASTALL = '/Volumes/Epoch9/tools/blast-2.2.26/bin/blastall'
####FORMATDB = '/Volumes/Epoch9/tools/blast-2.2.26/bin/formatdb'
####SAMTOOLS = '/Volumes/Epoch9/tools/samtools'
####
##### FILTERING PARAMETERS
####MIN_CLIP_LEN = 10
####RLEN_MIN     = 400
####RLEN_BUFF    = 8000
####
##### ALGORITHM PARAMETERS
####MIN_ALN_LEN          = 8	# blast alignment must be at least this long to be considered
####MAX_PRE_CLUSTER_DIST = 5	# when intiatally constructing graph, connect alignment edges at most this far apart
########                        #     --- MAX_PRE_CLUSTER_DIST must be < MIN_ALN_LEN, or resulting graph could possibly have cycles
####MAX_NOVEL_INS        = 500	# maximum length of novel insertion graph edges
####MAX_AMBIGUOUS_DEL    = 500	# maximum length for ambiguous deletion graph edges
####MAX_AMBIGUOUS_INS    = 1000	# maximum length for ambiguous insertion graph edges
####MIN_AMBIGUOUS_NOVEL  = 50	# the receiving node of an ambiguous deletion graph edge must span at least this many novel read positions
########
####NODE_PRUNE_THRESH    = 50	# if node spans less read positions than this, prune it if only connected edges are novel_ins
########
####PATH_MIN_EDGE_LEN    = 100	# the alignments on the flanking end of the best path must be at least this long, prune if lower
####MULTIPLICITY_FILTER  = 2	# an observed junction must occur in at least this many reads to be included in report
####
##### VARIOUS CONSTANTS
####JUNC_TO_STRING       = [chr(65)]
####PLOT_STUFF           = True


# Lets be civilized and read values in from a config file...
CORGI_CFG = SIM_PATH+'/corgi.cfg'

# EXECUTABLES
BLASTALL = getConfigVal(CORGI_CFG,'BLASTALL','str')
FORMATDB = getConfigVal(CORGI_CFG,'FORMATDB','str')
SAMTOOLS = getConfigVal(CORGI_CFG,'SAMTOOLS','str')
IGV_PATH = getConfigVal(CORGI_CFG,'IGV','str')
# FILTERING PARAMETERS
MIN_CLIP_LEN = getConfigVal(CORGI_CFG,'MIN_CLIP_LEN','int')
RLEN_MIN     = getConfigVal(CORGI_CFG,'RLEN_MIN','int')
RLEN_BUFF    = getConfigVal(CORGI_CFG,'RLEN_BUFF','int')
# ALGORITHM PARAMETERS
MIN_ALN_LEN          = getConfigVal(CORGI_CFG,'MIN_ALN_LEN','int')
MAX_PRE_CLUSTER_DIST = getConfigVal(CORGI_CFG,'MAX_PRE_CLUSTER_DIST','int')
MAX_NOVEL_INS        = getConfigVal(CORGI_CFG,'MAX_NOVEL_INS','int')
MAX_AMBIGUOUS_DEL    = getConfigVal(CORGI_CFG,'MAX_AMBIGUOUS_DEL','int')
MAX_AMBIGUOUS_INS    = getConfigVal(CORGI_CFG,'MAX_AMBIGUOUS_INS','int')
MIN_AMBIGUOUS_NOVEL  = getConfigVal(CORGI_CFG,'MIN_AMBIGUOUS_NOVEL','int')
NODE_PRUNE_THRESH    = getConfigVal(CORGI_CFG,'NODE_PRUNE_THRESH','int')
PATH_MIN_EDGE_LEN    = getConfigVal(CORGI_CFG,'PATH_MIN_EDGE_LEN','int')
MULTIPLICITY_FILTER  = getConfigVal(CORGI_CFG,'MULTIPLICITY_FILTER','int')
# VARIOUS CONSTANTS
JUNC_TO_STRING       = [chr(getConfigVal(CORGI_CFG,'JUNC_TO_STRING','int'))]
PLOT_STUFF           = getConfigVal(CORGI_CFG,'PLOT_STUFF','bool')

def processAlignment(f_in,report_dir):
	f = open(f_in,'r')
	byQuery = {}
	sLen = None
	sChr = None
	sPos = None
	sEnd = None
	for line in f:
		(q,l,qs,qf,ss,sf,ev,bs) = readBlastLine(line)
		if sLen == None:
			splt2 = line.split('\t')[1].split('_')
			sChr  = splt2[1]
			sLen  = int(splt2[3]) - int(splt2[2])
			sPos  = int(splt2[2])
			sEnd  = int(splt2[3])
		if q not in byQuery:
			byQuery[q] = []
		byQuery[q].append((qs,qf,ss,sf,bs))
	f.close()

	# for each read in the dataset...
	####all_paths      = []
	full_path_data = []
	current_ind    = 1
	total_ind      = len(byQuery.keys())
	for q in sorted(byQuery.keys()):
		print q, '('+str(current_ind)+'/'+str(len(byQuery))+')'

		#
		# CONNECT ALL ALIGNMENT EDGES THAT ARE WITHIN A SPECIFIED DISTANCE (ALONG THE READ, OF COURSE)
		#
		aln_graph        = {}	# aln_graph[(q_min, q_max, node_ind)] = set((connected_node,dist), ...)
		node_dat         = {}
		node_ind         = 0
		#dist_hist = {}
		for dat in byQuery[q]:
			(qs,qf,ss,sf,bs) = dat
			if abs(qf-qs) < MIN_ALN_LEN:
				continue
			#if abs(qf-qs) not in dist_hist:
			#	dist_hist[abs(qf-qs)] = 0
			#dist_hist[abs(qf-qs)] += 1
			node_dat[node_ind]   = (qs,qf,ss,sf,bs)		# ind corresponding to qs
			[q_min, q_max]       = sorted([qs,qf])
			aln_graph[(q_min,q_max,node_ind)] = set()
			node_ind += 1

		for node in aln_graph.keys():			# enumerate all connections node -> node2
			for node2 in aln_graph.keys():
				if node2 != node:
					wasnt_ambigDel_but_tried = 0
					# forbid connections into position 1 of reference (most likely erroneous connection because read starts before ref start)
					if node_dat[node2[2]][2] == 1 or node_dat[node2[2]][3] == 1:
						continue
					# similarly, forbid connections that go beyond the final reference position
					if node_dat[node[2]][2] == sLen or node_dat[node[2]][3] == sLen:
						continue
					# normal connections
					if abs(node2[0]-node[1]) <= MAX_PRE_CLUSTER_DIST:
						aln_graph[node].add((node2,abs(node2[0]-node[1]),'normal'))
					# novel insertions
					elif node2[0]-node[1] > MAX_PRE_CLUSTER_DIST and node2[0]-node[1] <= MAX_NOVEL_INS:
						aln_graph[node].add((node2,node2[0]-node[1],'novel_ins'))
					# ambiguous deletion (e.g. copy number loss)
					elif node2[0]-node[1] < -MAX_PRE_CLUSTER_DIST and node2[0]-node[1] >= -MAX_AMBIGUOUS_DEL:
						wasnt_ambigDel_but_tried = 1
						(ambig_del_forward, ambig_del_within, ambig_del_novel) = (False, False, False)
						# must make forward progress along reference
						if node_dat[node2[2]][2] > node_dat[node[2]][3]:
							ambig_del_forward = True
						# is the landing point for the connection within the read-position span of where we're coming from?
						if node2[0] <= node[1] and node2[0] >= node[0]:
							ambig_del_within = True
						# does landing point involve a span of read sequence sufficiently unique from where we're coming from?
						if node2[1]-node[1] >= MIN_AMBIGUOUS_NOVEL:
							ambig_del_novel = True
						# all conditions met?
						if all([ambig_del_forward, ambig_del_within, ambig_del_novel]):
							aln_graph[node].add((node2,abs(node2[0]-node[1]),'ambiguous_del'))
							wasnt_ambigDel_but_tried = 2
					# ambiguous insertion (e.g. copy number gain)
					if wasnt_ambigDel_but_tried == 1 and node2[0]-node[1] < -MAX_PRE_CLUSTER_DIST and node2[0]-node[1] >= -MAX_AMBIGUOUS_INS:
						(ambig_ins_backward, ambig_ins_within, ambig_ins_novel) = (False, False, False)
						# must make backward progress along the reference
						if node_dat[node2[2]][2] < node_dat[node[2]][3]:
							ambig_ins_backward = True
						# is the landing point for the connection within the read-position span of where we're coming from?
						if node2[0] <= node[1] and node2[0] >= node[0]:
							ambig_ins_within = True
						# does landing point involve a span of read sequence sufficiently unique from where we're coming from?
						if node2[1]-node[1] >= MIN_AMBIGUOUS_NOVEL:
							ambig_ins_novel = True
						# all conditions met?
						if all([ambig_ins_backward, ambig_ins_within, ambig_ins_novel]):
							aln_graph[node].add((node2,abs(node2[0]-node[1]),'ambiguous_ins'))

		#print 'len(aln_graph):', len(aln_graph), len(aln_graph)*len(aln_graph), '-->',

		# prune nodes to accelerate traversal
		# --- remove node if node_prize is below a threshold and all edges to/from node are novel_ins (or ambiguous_del/ins)
		for node in aln_graph.keys():
			if node[1]-node[0] < NODE_PRUNE_THRESH:
				only_nov_ins = True
				for node2 in aln_graph[node]:
					if node2[2] == 'normal':
						only_nov_ins = False
				for node2 in aln_graph.keys():
					for node3 in aln_graph[node2]:
						if node3[0] == node and node3[2] == 'normal':
							only_nov_ins = False
				if only_nov_ins:
					aln_graph[node] = set()
					for node2 in aln_graph.keys():
						aln_graph[node2] = set([n for n in aln_graph[node2] if n[0] != node])
		##### --- remove nodes of degree zero
		####myCount = 0
		####node_rename_dict = {}
		####for node in aln_graph.keys():
		####	if len(aln_graph[node]) == 0:
		####		notFoundElsewhere = True
		####		for node2 in aln_graph.keys():
		####			if node != node2:
		####				if node in [n[0] for n in aln_graph[node2]]:
		####					notFoundElsewhere = False
		####					break
		####		if notFoundElsewhere:
		####			del aln_graph[node]
		####		else:
		####			node_rename_dict[node[0]] = myCount
		####			myCount += 1
		####	else:
		####		node_rename_dict[node[0]] = myCount
		####		myCount += 1
		####print len(aln_graph), sum([len(aln_graph[k]) for k in aln_graph.keys()])
		##### --- renumber nodes to reflect reduced size
		#####print node_rename_dict
		####continue

		# convert to a more sensible graph structure for traversal
		score_graph = {}
		for node in aln_graph.keys():
			score_graph[(node[2],node[1]-node[0])] = set()
			for node2 in aln_graph[node]:
				score_graph[(node[2],node[1]-node[0])].add(((node2[0][2],node2[0][1]-node2[0][0]),node2[1]))

		#
		# create DAG from score_graph and sort topologically
		#
		####
		####	prize(node) = bp spanned (in read) by this alignment
		####
		####	weight(node1,node2) = [distance (in read) between edges of alignment] - 2 * [bp shared (in read) by alignments]
		####
		graph_unsorted = []
		graph_weights  = [[None for m in xrange(len(score_graph))] for n in xrange(len(score_graph))]
		graph_prize    = {n[0]:n[1] for n in score_graph.keys()}
		for node in score_graph.keys():
			tempList = []
			for node2 in score_graph[node]:
				(qs1,qf1,qs2,qf2) = (node_dat[node[0]][0],node_dat[node[0]][1],node_dat[node2[0][0]][0],node_dat[node2[0][0]][1])
				bp_overlap_in_read = 0
				for i in xrange(qs1,qf1+1):
					if i >= qs2 and i <= qf2:
						bp_overlap_in_read += 1
				graph_weights[node[0]][node2[0][0]] = node2[0][1] - node2[1] - 2*bp_overlap_in_read
				tempList.append(node2[0][0])
			graph_unsorted.append((node[0],[n for n in tempList]))
		# topological sort
		graph_sorted = []
		graph_unsorted = dict(graph_unsorted)
		while graph_unsorted:
			acyclic = False
			for node, edges in graph_unsorted.items():
				for edge in edges:
					if edge in graph_unsorted:
						break
				else:
					acyclic = True
					del graph_unsorted[node]
					graph_sorted.append((node, edges))
			if not acyclic:
				raise RuntimeError("A cyclic dependency occurred")

		#
		# DAG_longestPath
		#
		tt = time.time()
		[bestScore,bestPath] = exhaustive_DAG(graph_sorted,graph_weights,graph_prize)
		#print 'BEST_SCORE:',bestScore
		#print 'BEST_PATH: ',bestPath
		#print time.time()-tt,'(sec)'
		#exit(1)

		# prune flanking alignments of best path, if they are too small
		if len(bestPath) > 1:
			delList = []
			for i in xrange(len(bestPath)):
				myLen = node_dat[bestPath[i]][1] - node_dat[bestPath[i]][0]
				if abs(myLen) >= PATH_MIN_EDGE_LEN:
					break
				else:
					delList.append(i)
			for i in xrange(len(bestPath)-1,-1,-1):
				myLen = node_dat[bestPath[i]][1] - node_dat[bestPath[i]][0]
				if abs(myLen) >= PATH_MIN_EDGE_LEN:
					break
				else:
					delList.append(i)
			delList = sorted(list(set(delList)),reverse=True)
			for n in delList:
				del bestPath[n]

		# from the best path, create symbol string for each detected junction
		if len(bestPath) > 1:
			####myConn = []
			####for i in xrange(len(bestPath)-1):
			####	myConn.append((node_dat[bestPath[i]][3],'-->',node_dat[bestPath[i+1]][2]))
			####all_paths.append([n for n in myConn])
			full_path_data.append([node_dat[n] for n in bestPath])

			##### deal with ambiguous deletions
			####ambig_del_byNodeKey = {}
			####ambigDel = [(full_path_data[-1][n+1][0] - full_path_data[-1][n][1] >= -MAX_AMBIGUOUS_DEL and full_path_data[-1][n+1][0] - full_path_data[-1][n][1] < -MAX_PRE_CLUSTER_DIST) for n in xrange(len(full_path_data[-1])-1)]
			####while any(ambigDel):
			####	i = ambigDel.index(True)
			####	# TEMPORARY: REDUCE AMBIGUOUS DELETIONS A SINGLE DELETION
			####	####print full_path_data[-1][i], full_path_data[-1][i+1]
			####	####delta = full_path_data[-1][i][1] - full_path_data[-1][i+1][0] + 1
			####	####isForward = (full_path_data[-1][i+1][3] > full_path_data[-1][i+1][2])
			####	####(d1,d2,d3,d4,d5) = full_path_data[-1][i+1]
			####	####if isForward:
			####	####	full_path_data[-1][i+1] = (d1+delta,d2,d3+delta,d4,d5)
			####	####else:
			####	####	full_path_data[-1][i+1] = (d1+delta,d2,d3,d4+delta,d5)
			####
			####	# hypothetically, an alternative would be to say that reference region has shrunk, ambiguously:
			####	delta = full_path_data[-1][i][1] - full_path_data[-1][i+1][0]
			####	(d1,d2,d3,d4,d5) = full_path_data[-1][i]
			####	isForward        = (d4 > d3)
			####	if isForward:
			####		full_path_data[-1][i] = (d1,d2-delta,d3,d4-delta,d5)
			####	else:
			####		full_path_data[-1][i] = (d1,d2-delta,d3,d4+delta,d5)
			####	(d1,d2,d3,d4,d5) = full_path_data[-1][i+1]
			####	isForward        = (d4 > d3)
			####	if isForward:
			####		full_path_data[-1][i+1] = (d1+delta,d2,d3+delta,d4,d5)
			####	else:
			####		full_path_data[-1][i+1] = (d1+delta,d2,d3-delta,d4,d5)
			####	ambig_del_byNodeKey[(full_path_data[-1][i],full_path_data[-1][i+1])] = True
			####
			####	ambigDel = [(full_path_data[-1][n+1][0] - full_path_data[-1][n][1] >= -MAX_AMBIGUOUS_DEL and full_path_data[-1][n+1][0] - full_path_data[-1][n][1]) < -MAX_PRE_CLUSTER_DIST for n in xrange(len(full_path_data[-1])-1)]


		if PLOT_STUFF:
			PLOT_DIR = report_dir+'alignments/'
			makedir(PLOT_DIR)
			fig = mpl.figure(0,figsize=(14,6))
			mpl.subplot(121)
			#ax = mpl.axes()
			plot_me_later = []
			for node in aln_graph.keys():
				(qs,qf,ss,sf,bs) = node_dat[node[2]]
				x = [ss,sf]
				y = [qs,qf]
				mpl.plot(x,y,'-k',linewidth=2)
				for node2 in aln_graph[node]:
					x = [sf,node_dat[node2[0][2]][2]]
					y = [node[1],node2[0][0]]
					if node2[2] == 'normal':
						plot_me_later.append([[x[0],x[1]],[y[0],y[1]]])
					elif node2[2] == 'novel_ins':
						mpl.plot(x,y,'--r',linewidth=1)
					elif node2[2] == 'ambiguous_del':
						mpl.plot(x,y,'--m',linewidth=2)
					elif node2[2] == 'ambiguous_ins':
						mpl.plot(x,y,'--b',linewidth=2)
			for n in plot_me_later:
				mpl.plot(n[0],n[1],'-g',linewidth=2)
			mpl.grid()
			mpl.title('read '+str(current_ind)+'/'+str(total_ind))
			mpl.xlabel('ref position')
			mpl.ylabel('read position')
			current_ind += 1

			mpl.subplot(122)
			for i in xrange(len(bestPath)):
				(qs,qf,ss,sf,bs) = node_dat[bestPath[i]]
				x = [ss,sf]
				y = [qs,qf]
				mpl.xlabel('ref position')
				mpl.ylabel('read position')
				mpl.plot(x,y,'-k',linewidth=2)
			mpl.xlim([0,sLen])
			#mpl.show()
			mpl.savefig(PLOT_DIR+'read_'+str(current_ind-1)+'.pdf')
			#mpl.clf()
			mpl.close(fig)

	##### prune regions that don't pass multiplicity filter
	####if len(full_path_data):
	####	passedFilter = {}
	####	keyList = []
	####	for i in xrange(len(full_path_data)):
	####		key = [(full_path_data[i][j][3],full_path_data[i][j+1][2]) for j in xrange(len(full_path_data[i])-1)]
	####		keyList.append(tuple(key))
	####		if keyList[-1] in passedFilter:
	####			passedFilter[keyList[-1]] = True
	####		else:
	####			passedFilter[keyList[-1]] = False
	####	for i in xrange(len(full_path_data)-1,-1,-1):
	####		if not passedFilter[keyList[i]]:
	####			del full_path_data[i]

	#
	# infer chunked regions
	#
	if len(full_path_data):
		region_map = {}
		for i in xrange(len(full_path_data)):
			print '---', full_path_data[i]
			for j in xrange(len(full_path_data[i])-1):
				key = (full_path_data[i][j][3],full_path_data[i][j+1][2])
				if key not in region_map:
					region_map[key] = 0
				region_map[key] += 1
		# apply multiplicity filter
		region_map = {k:region_map[k] for k in region_map.keys() if region_map[k] >= MULTIPLICITY_FILTER}
		#print 'region_map:',region_map
		boundaries = sorted(list(set([n for n2 in region_map.keys()+[(1,sLen)] for n in n2])))
		print 'boundaries:',boundaries
		region_map = [('anchor1',boundaries[0],boundaries[1]-1)]
		for i in xrange(1,len(boundaries)-2):
			region_map.append((chr(64+i),boundaries[i],boundaries[i+1]-1))
		region_map.append(('anchor2',boundaries[-2],boundaries[-1]))
		#rev_region_map = {region_map[i][0]:i for i in xrange(len(region_map))}
		#for i in xrange(1,len(rev_region_map)-1):
		#	rev_region_map[region_map[i][0]+'*'] = i
		#print 'region_map:',region_map
		#print 'rev_region_map:',rev_region_map

		sorted_ps = [n for n in get_path_strings(full_path_data,region_map,MAX_PRE_CLUSTER_DIST) if n[0] >= MULTIPLICITY_FILTER]
		print '-sorted_ps_1-'
		for n in sorted_ps:
			print n

		# redo region chunking and string assignment based on patterns that passed multiplicity filter
		boundaries = [1,sLen]
		for i in xrange(len(sorted_ps)):
			for j in xrange(len(sorted_ps[i][1])):
				#myRegion = [n[0] for n in region_map].index(sorted_ps[i][1][j][2])
				myRegion = sorted_ps[i][1][j][0]
				toAdd = [region_map[myRegion][1]]
				if j < len(sorted_ps[i][1])-1 and region_map[myRegion][2]+1 < sLen:
					toAdd.append(region_map[myRegion][2]+1)
				boundaries.extend(toAdd)
		boundaries = sorted(list(set(boundaries)))
		print 'boundaries (again):',boundaries
		region_map = [('anchor1',boundaries[0],boundaries[1]-1)]
		for i in xrange(1,len(boundaries)-2):
			region_map.append((chr(64+i),boundaries[i],boundaries[i+1]-1))
		region_map.append(('anchor2',boundaries[-2],boundaries[-1]))
		rev_region_map = {region_map[i][0]:i for i in xrange(len(region_map))}
		for i in xrange(1,len(rev_region_map)-1):
			rev_region_map[region_map[i][0]+'*'] = i
		print 'region_map (again):',region_map
		# remove any full_path_data entries that contain breakpoints not in our new cleaned boundaries list
		for i in xrange(len(full_path_data)-1,-1,-1):
			anyBad = False
			for j in xrange(len(full_path_data[i])-1):
				if full_path_data[i][j][3] not in boundaries or full_path_data[i][j+1][2] not in boundaries:
					anyBad = True
			if anyBad:
				del full_path_data[i]
				
		print 'full_path_data:'
		for n in full_path_data:
			print n
		sorted_ps = [n for n in get_path_strings(full_path_data,region_map,MAX_PRE_CLUSTER_DIST) if n[0] >= MULTIPLICITY_FILTER]
		print '-sorted_ps_2-'
		for n in sorted_ps:
			print n

		#
		# string annotation
		#
		# construct reference string to feed into interpretation function
		ref_int = [n[0] for n in region_map]
		# novel region names, indexed by (prev_region, novel_len)
		novel_dict = {}
		# ambig gap info, indexed by (region, prev_region)
		ambig_dict = {}
		if len(region_map) > 2:
			novel_val = chr(max([ord(n[0]) for n in region_map[1:-1]])+1)
		else:
			novel_val = JUNC_TO_STRING[0]
		# for each observed rearrangement...
		sample_strings = []
		for i in xrange(len(sorted_ps)):
			# construct sample string to feed into interpretation function
			samp_int   = []
			for j in xrange(len(sorted_ps[i][1])):
				samp_int.append(sorted_ps[i][1][j][2])
				if sorted_ps[i][1][j][3][1]:
					novel_key = (samp_int[-1],sorted_ps[i][1][j][3][0])
					if novel_key not in novel_dict:
						novel_dict[novel_key] = novel_val
						novel_val = chr(ord(novel_val)+1)
					samp_int.append(novel_dict[novel_key])
				if sorted_ps[i][1][j][3][2]:
					if j < len(sorted_ps[i][1])-1:
						ambig_key = (sorted_ps[i][1][j][2],sorted_ps[i][1][j+1][2])
						ambig_dict[ambig_key] = sorted_ps[i][1][j][3][0]
			sample_strings.append([n for n in samp_int])

			####novel_list = [n for n in samp_int if n in novel_dict.values()]
			####
			####print sorted_ps[i]
			####print '===',ref_int,samp_int,novel_list
			#####print sorted_ps[i][0], [m[2] for m in sorted_ps[i][1]]
			#####print sorted_ps[i]
			####
			##### if event is enclosed by anchors, interpret whole event normally
			####if samp_int[0] == 'anchor1' and samp_int[-1] == 'anchor2':
			####	(path_int,clust_int,str_int) = interpret_rearrangement(ref_int,samp_int,novel_list)
			##### otherwise we need to deal with the fact that we're only partially covering the event
			####else:
			####	str_int = 'incomplete'
			####
			####interpretation_list.append(str_int)
		rev_novel_dict = {novel_dict[k]:k for k in novel_dict.keys()}
		ambig_plotting = {}
		for k in ambig_dict.keys():
			ambig_plotting[k[0]] = (ambig_dict[k],True)
			ambig_plotting[k[1]] = (ambig_dict[k],False)
		print 'sample_strings:',sample_strings

		# rename inverted regions so that all string entries are a single letter (to facilitate MSA)
		inv_dict = {}
		strings_for_msa = []
		for i in xrange(len(sample_strings)):
			msaStr = []
			for m in sample_strings[i]:
				if m == 'anchor1' or m == 'anchor2' or m[-1] == '*':
					if m not in inv_dict:
						inv_dict[m] = novel_val
						novel_val = chr(ord(novel_val)+1)
					msaStr.append(inv_dict[m])
				else:
					inv_dict[m] = m
					msaStr.append(m)
			strings_for_msa.append(''.join(msaStr))
		rev_inv_dict = {inv_dict[k]:k for k in inv_dict.keys()}
		print 'inv_dict:',inv_dict
		print 'rev_inv_dict:',rev_inv_dict
		print 'strings_for_msa:',strings_for_msa,'-->',
		# remove strings that are subsets of other strings, for efficiency..
		if len(strings_for_msa) > 1:
			delList = []
			for i in xrange(len(strings_for_msa)):
				for j in xrange(len(strings_for_msa)):
					if i == j:
						continue
					if strings_for_msa[i] in strings_for_msa[j] and len(strings_for_msa[i]) < len(strings_for_msa[j]):
						#print '*',(i,j), strings_for_msa[i], strings_for_msa[j]
						delList.append(i)
			delList = sorted(list(set(delList)),reverse=True)
			for i in delList:
				del strings_for_msa[i]
		print strings_for_msa

		# if only 1 string is present, no need to do any fancy MSA steps
		sorted_events_final = []
		interpretation_list = []
		if len(strings_for_msa) == 1:
			multiplicity_list = [sorted_ps[0][0]]
			samp_int   = [rev_inv_dict[n] for n in strings_for_msa[0]]
			novel_list = [n for n in samp_int if n in rev_novel_dict]
			s_ps_entry = []
			for i in xrange(len(samp_int)):
				if samp_int[i] not in novel_list:
					if i < len(samp_int)-1 and samp_int[i+1] in novel_list:
						novelDat = (rev_novel_dict[samp_int[i+1]][1],True,False)
					else:
						novelDat = (0,False,False)	# WE'RE LOSING INFO HERE. TODO: FIND A WAY TO PRESERVE UNCERTAINTY GAP IF NOT NOVEL INS
					# check for ambig
					if i < len(samp_int)-1 and (samp_int[i],samp_int[i+1]) in ambig_dict:
						novelDat = (ambig_dict[(samp_int[i],samp_int[i+1])],False,True)
					s_ps_entry.append((rev_region_map[samp_int[i]],samp_int[i][-1]!='*',samp_int[i],novelDat))
			sorted_events_final.append((str(sum(multiplicity_list))+'*'*(len(multiplicity_list)>1),tuple(s_ps_entry)))
			(path_int,clust_int,str_int) = interpret_rearrangement(ref_int,samp_int,novel_list)
			interpretation_list.append(str_int)
		# otherwise, cluster msa strings by shared junctions and determine consensus sequences...
		else:
			msa_adj = [[0 for n in strings_for_msa] for m in strings_for_msa]
			for i in xrange(len(msa_adj)):
				for j in xrange(len(msa_adj)):
					if i != j:
						# I'm a substring of you
						if strings_for_msa[i] in strings_for_msa[j]:
							msa_adj[i][j] = 1
							msa_adj[j][i] = 1
						# my suffix is your prefix
						# (test all possible suffix lengths)
						for suffixLen in xrange(2,min(len(strings_for_msa[i]),len(strings_for_msa[j]))+1):
							if strings_for_msa[i][-suffixLen:] == strings_for_msa[j][:suffixLen]:
								msa_adj[i][j] = 1
								msa_adj[j][i] = 1
						#for k in xrange(0,len(strings_for_msa[i])-1):
						#	if strings_for_msa[i][k:k+2] in strings_for_msa[j]:
						#		msa_adj[i][j] = 1
			msa_clusters = get_connected_subgraphs(msa_adj)
			print msa_clusters
			# perform msa
			for i in xrange(len(msa_clusters)):
				multiplicity_list = [sorted_ps[n][0] for n in msa_clusters[i]]
				consensus_seq     = msa_consensus([strings_for_msa[n] for n in msa_clusters[i]])
				# trim gap characters from ends of consensus, if present
				consensus_seq     = [n for n in consensus_seq if n[0] != '-']
				print 'consensus_seq:',consensus_seq
				samp_int   = [rev_inv_dict[n[0]] for n in consensus_seq]
				print 'samp_int (msa):',samp_int
				# code reuse, yay I'm bad at programming!
				novel_list = [n for n in samp_int if n in rev_novel_dict]
				s_ps_entry = []
				for i in xrange(len(samp_int)):
					if samp_int[i] not in novel_list:
						if i < len(samp_int)-1 and samp_int[i+1] in novel_list:
							novelDat = (rev_novel_dict[samp_int[i+1]][1],True,False)
						else:
							novelDat = (0,False,False)	# WE'RE LOSING INFO HERE. TODO: FIND A WAY TO PRESERVE UNCERTAINTY GAP IF NOT NOVEL INS
						# check for ambig
						if i < len(samp_int)-1 and (samp_int[i],samp_int[i+1]) in ambig_dict:
							novelDat = (ambig_dict[(samp_int[i],samp_int[i+1])],False,True)
						s_ps_entry.append((rev_region_map[samp_int[i]],samp_int[i][-1]!='*',samp_int[i],novelDat))
				sorted_events_final.append((str(sum(multiplicity_list))+'*'*(len(multiplicity_list)>1),tuple(s_ps_entry)))
				print 'interpret_rearrangement input:',ref_int, samp_int, novel_list
				(path_int,clust_int,str_int) = interpret_rearrangement(ref_int,samp_int,novel_list)
				interpretation_list.append(str_int)

		#
		# GENERATE HTML REPORT
		#
		print sorted_events_final
		print interpretation_list
		gen_plots((sChr,sPos,sEnd),region_map,sorted_events_final,interpretation_list,ambig_plotting,report_dir)

		#
		# GENERATE OUTPUT FILE
		#
		ouf = report_dir+'results.p'
		pickle.dump([(sChr,sPos,sEnd),region_map,sorted_events_final,interpretation_list], open(ouf,"wb"))


""" ***************************************************************************
********                                                                 ******
*******                             MAIN()                              *******
******                                                                 ********
*************************************************************************** """

def main():

	#	Parse input arguments
	#
	parser = argparse.ArgumentParser(description='plotLongRead5.py')
	parser.add_argument('-r', type=str, required=True,  metavar='<str>',                    help="* ref.fa")
	parser.add_argument('-b', type=str, required=True,  metavar='<str>',                    help="* input.bam")
	parser.add_argument('-c', type=str, required=True,  metavar='<str>',                    help="* chromosome")
	parser.add_argument('-p', type=int, required=True,  metavar=('<int>','<int>'), nargs=2, help='* start & end coordinates')
	parser.add_argument('-o', type=str, required=True,  metavar='<str>',                    help="* outputDir/")
	parser.add_argument('--skip-igv',   required=False, action='store_true', default=False, help='skip IGV screenshot')
	parser.add_argument('--skip-bam',   required=False, action='store_true', default=False, help='skip read extraction')
	args = parser.parse_args()

	(REF_FILE, INPUT_BAM, CHROM, POS_S, POS_E, OUT_DIR) = (args.r, args.b, args.c, args.p[0], args.p[1], args.o)
	(PROCESS_BAM, PROCESS_IGV) = (not(args.skip_bam), not(args.skip_igv))

	if OUT_DIR[-1] != '/':
		OUT_DIR += '/'
	makedir(OUT_DIR)
	OUT_DIR   += 'corgi_'+CHROM+'_'+str(POS_S)+'_'+str(POS_E)+'/'
	TEMP_DIR  = OUT_DIR
	makedir(TEMP_DIR)

	if PROCESS_BAM:
		print 'grabbing reads from input bam...'
		#
		refChunk = TEMP_DIR+'refChunk.fa'
		writeRefSection(REF_FILE, CHROM, POS_S, POS_E, refChunk)
		exe(FORMATDB + ' -i ' + refChunk + ' -l ' + TEMP_DIR + 'formatdb.log -p F')

		#
		tempF1   = TEMP_DIR+'allReads.sam'
		samPos_s = max([0,POS_S - RLEN_BUFF])
		samPos_e = POS_E + RLEN_BUFF
		exe(SAMTOOLS + ' view ' + '-o ' + tempF1 + ' ' + INPUT_BAM + ' ' + CHROM + ':' + str(samPos_s) + '-' + str(samPos_e+RLEN_BUFF))
		# select only reads with sufficient clipped content at either end, also by mapping position
		tempF1_f1 = TEMP_DIR+'allReads_withClip.sam'
		tempF1_f2 = TEMP_DIR+'allReads_filtered.sam'
		tempF1_f3 = TEMP_DIR+'allReads_filtered2.sam'
		filterSamByClipContent(tempF1,tempF1_f1,MIN_CLIP_LEN,indelLen=MIN_ALN_LEN)
		filterSamByRefOverlap(tempF1_f1,tempF1_f2,POS_S,POS_E)
		filterSamByReadLen(tempF1_f2,tempF1_f3,rMin=RLEN_MIN)

		#
		tempF2 = TEMP_DIR+'allReads_filtered.fa'
		sam2fa(tempF1_f3,tempF2)

		#
		tempF3 = TEMP_DIR+'blastOut.txt'
		exe(BLASTALL + ' -p blastn -d ' + refChunk + ' -i ' + tempF2 + ' -m8 -W5 -a4 -FF > ' + tempF3)
	else:
		tempF3 = TEMP_DIR+'blastOut.txt'

	#
	if PROCESS_IGV:
		print 'grabbing alignment screenshot...'
		view_IGV(IGV_PATH,REF_FILE,INPUT_BAM,CHROM,POS_S,POS_E,screenShotDir=OUT_DIR)
	#exit(1)
	print 'analyzing reads...'
	RESOURCE_DIR = OUT_DIR+'resources/'
	BOKEH_DIR    = '/'.join(os.path.realpath(__file__).split('/')[:-1])+'/plotting_resources/'
	BOKEH_FILES  = ["bokeh-0.12.6.min.css","bokeh-0.12.6.min.js"]
	makedir(RESOURCE_DIR)
	for n in BOKEH_FILES:
		exe('cp '+BOKEH_DIR+n+' '+RESOURCE_DIR+n)
	processAlignment(tempF3,OUT_DIR)


if __name__ == '__main__':
	main()