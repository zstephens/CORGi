import copy


def getConfigVal(cfg,varName,varType):
	myItem = None
	f = open(cfg,'r')
	for line in f:
		if line[0] != '#' and len(line.strip()):
			splt = [n for n in line.strip().split(' ') if len(n)]
			if splt[0] == varName:
				myItem = splt[1]
				break
	f.close()
	if myItem == None:
		print '\nError: Config variable "'+varName+'" not found.\n'
		exit(1)
	else:
		if varType == 'str':
			return myItem
		elif varType == 'int':
			return int(myItem)
		elif varType == 'float':
			return float(myItem)
		elif varType == 'bool':
			if myItem.upper() == 'TRUE':
				return True
			else:
				return False
		else:
			print '\nError: Unknown config variable type: "'+varName+'"\n'
			exit(1)

# if a read encompasses all the ref positions plus or minus this amount around the
# beginning of a reference region (at the end of a long read alignment), lets append
# it to the end of the path string
REF_REGION_EXTEND_SPAN = 20

def get_path_strings(full_path_data, region_map, MAX_PRE_CLUSTER_DIST):
	path_string_count = {}
	for i in xrange(len(full_path_data)):
		path_string = []
		for j in xrange(len(full_path_data[i])):
			# region is reverse-complemented if ref coordinates are inverted
			isForward = (full_path_data[i][j][3] > full_path_data[i][j][2])
			myInd     = min([full_path_data[i][j][3], full_path_data[i][j][2]])
			myR_start = None
			for k in xrange(len(region_map)):
				if myInd >= region_map[k][1] and myInd <= region_map[k][2]:
					myR_start = k
					break
			# find out what regions we span
			myInd     = max([full_path_data[i][j][3], full_path_data[i][j][2]]) - 1
			myR_span  = None
			for k in xrange(len(region_map)):
				#print myInd, region_map[k]
				if myInd > region_map[k][1] and myInd <= region_map[k][2]:
					myR_span = k
					break
			rInds = range(myR_start,myR_span+1)
			#print 'rInds:',myR_start,myR_span+1,rInds
			if isForward == False:
				rInds.reverse()

			myGap = 0
			if j < len(full_path_data[i])-1:
				myGap = full_path_data[i][j+1][0] - full_path_data[i][j][1]
			# path_string[i] = (region_index, isForward, region_str, (connection_gap, isNovel, isAmbig))
			for myRegion in rInds[:-1]:
				path_string.append((myRegion,isForward,region_map[myRegion][0]+'*'*(~isForward+2),(0,False,False)))
			myRegion = rInds[-1]
			path_string.append((myRegion,isForward,region_map[myRegion][0]+'*'*(~isForward+2),(myGap,(myGap>MAX_PRE_CLUSTER_DIST),(myGap<-MAX_PRE_CLUSTER_DIST))))
			#print 'path_string:',path_string

		##### see if the end of a long read spans additional reference regions
		##### -- at the moment I'm only interested in forward reigons that might span to anchor2...
		####if path_string[-1][1]:
		####	myRefSpan = (full_path_data[i][-1][2],full_path_data[i][-1][3])
		####	#print 'myRefSpan:',myRefSpan
		####	#print 'all possible regions after:',region_map[path_string[-1][0]][0],'...'
		####	for j in xrange(path_string[-1][0]+1,len(region_map)):
		####		#print region_map[j]
		####		# allow us to begin in small regions
		####		tl = max([myRefSpan[0],region_map[j][1]-REF_REGION_EXTEND_SPAN])
		####		# but be more strict in the larger ones
		####		tr = region_map[j][1]+REF_REGION_EXTEND_SPAN
		####		if (tl >= myRefSpan[0] and tl < myRefSpan[1]) and (tr >= myRefSpan[0] and tr < myRefSpan[1]):
		####			path_string.append((j,True,region_map[j][0],(0,(False))))

		# count up the junction patterns we encounter
		if tuple(path_string) not in path_string_count:
			path_string_count[tuple(path_string)] = 0
		path_string_count[tuple(path_string)] += 1

	return sorted([(path_string_count[n],n) for n in path_string_count.keys()],reverse=True)



def str_inv(s):
	if s[-1] != '*':
		return s+'*'
	elif len(s) >= 2:
		return s[:-1]
	else:
		print 'Error: Bad string input in str_inv()'
		exit(1)

def ref_dist(r,s):
	allCounts = {}
	for n in r:
		if n not in allCounts:
			allCounts[n] = [0,0]
		allCounts[n][0] += 1
	for n in s:
		if n not in allCounts:
			allCounts[n] = [0,0]
		allCounts[n][1] += 1
	dist = sum([abs(v[0]-v[1]) for v in allCounts.values()])
	return dist

TOO_MANY_REGIONS = 8
TOO_MANY_NOVEL   = 4

def interpret_rearrangement(ref,samp,novel=[]):
	# via exhaustive search, find the simplest combination of
	# insertions/deletions/inversions that describe the aggregate event
	#
	# input example: ref  = ['anchor','A','B','C','anchor']
	#                samp = ['anchor','B*','anchor']
	#
	#                novel = list of regions not originating in ref
	#

	if len(novel) >= TOO_MANY_NOVEL:
		print 'warning: too many novel regions, skipping interpretation...'
		return ([],[],'[unknown]')
	regionCount = list(set([n for n in ref+samp if 'anchor' not in n]))
	if len(regionCount) >= TOO_MANY_REGIONS:
		print 'warning: too many reference regions, skipping interpretation...'
		return ([],[],'[unknown]')

	if samp[0] != 'anchor1':
		print 'warning: event does not start in anchor sequence'
		# append regions that precede the region we start at
		# THIS MAY OR MAY OR MAY NOT BE A GOOD IDEA. BETTER STRATEGY WOULD BE TO DETERMINE
		# PHASING INFORMATION AND TRY TO BRIDGE READS TOGETHER --> HIGHLY DIFFICULT
		query = samp[0]
		if query[-1] == '*':
			query = query[:-1]
		samp = ref[:ref.index(query)] + samp
	if samp[-1] != 'anchor2':
		print 'warning: event does not end in anchor sequence'
		# repeat similarly if samp does not end as expected
		query = samp[-1]
		if query[-1] == '*':
			query = query[:-1]
		samp = samp + ref[ref.index(query)+1:]

	MAX_DEPTH = 100
	MAX_REGION_LEN = len(ref)-2
	queue    = [[[n for n in ref],[]]]	# [current_str,path]
	allPaths = []
	#visited  = {}
	while queue:
		[myRef,myPath] = queue.pop(0)
		#if tuple(myRef) in visited:
		#	continue
		#visited[tuple(myRef)] = True
		current_dist = ref_dist(myRef,samp)
		#print (len(myPath),myPath,current_dist)

		for regionLen in xrange(1,MAX_REGION_LEN+1):
			for i in xrange(1,len(myRef)-regionLen):
				#
				# apply deletion (if current string is bigger (OR EQUAL LEN) than samp)
				#
				if len(myRef) >= len(samp):
					myCopy = [n for n in myRef]
					del myCopy[i:i+regionLen]
					if ref_dist(myCopy,samp) < current_dist:
						newPath = myPath + [('del',(i,i+regionLen))]
						if myCopy == samp:
							queue = [n for n in queue if len(n[1]) < MAX_DEPTH]
							allPaths.append([n for n in newPath])
							MAX_DEPTH = len(newPath)
						if len(newPath) < MAX_DEPTH:
							queue.append([[n for n in myCopy],[n for n in newPath]])

				#
				# apply inversion (only allow a single region at a time be inverted)
				#
				if regionLen == 1:
					myCopy = [n for n in myRef]
					myCopy[i:i+regionLen] = [str_inv(n) for n in myCopy[i:i+regionLen][::-1]]
					if ref_dist(myCopy,samp) < current_dist:
						newPath = myPath + [('inv',(i,i+regionLen))]
						if myCopy == samp:
							queue = [n for n in queue if len(n[1]) < MAX_DEPTH]
							allPaths.append([n for n in newPath])
							MAX_DEPTH = len(newPath)
						if len(newPath) < MAX_DEPTH:
							queue.append([[n for n in myCopy],[n for n in newPath]])
				#
				# apply insertions at each possible position (if current string is smaller than samp)
				#
				if len(myRef) < len(samp):
					myRegions = [myRef[i:i+regionLen]]
					# for regions of length one, also allow inverted duplications
					if regionLen == 1:
						myRegions.append([str_inv(myRegions[0][0])])
					for myRegion in myRegions:
						for j in xrange(1,len(myRef)):
							myCopy = [n for n in myRef]
							myCopy = myCopy[:j] + myRegion + myCopy[j:]
							dup_prefix = (myRegion[0][-1]=='*')*'inv-'
							if j == i or j == i+regionLen:
								dup_suffix = 'tdup'
							elif j < i or j > i+regionLen:
								dup_suffix = 'ddup'
							else:
								dup_suffix = 'dup'
							if ref_dist(myCopy,samp) < current_dist:
								newPath = myPath + [(dup_prefix+dup_suffix,(j,myRegion))]
								if myCopy == samp:
									queue = [n for n in queue if len(n[1]) < MAX_DEPTH]
									allPaths.append([n for n in newPath])
									MAX_DEPTH = len(newPath)
								if len(newPath) < MAX_DEPTH:
									queue.append([[n for n in myCopy],[n for n in newPath]])

		#
		# apply novel insertions
		#
		for i in xrange(len(novel)):
			myRegion = [novel[i]]
			for j in xrange(1,len(myRef)):
				myCopy = [n for n in myRef]
				myCopy = myCopy[:j] + myRegion + myCopy[j:]
				newPath = myPath + [('ins',(j,myRegion))]
				if myCopy == samp:
					queue = [n for n in queue if len(n[1]) < MAX_DEPTH]
					allPaths.append([n for n in newPath])
					MAX_DEPTH = len(newPath)
				if len(newPath) < MAX_DEPTH:
					queue.append([[n for n in myCopy],[n for n in newPath]])

	# we found solutions?
	if len(allPaths):
		# heirarchally cluster alterations
		all_dependencies  = []
		all_refAlterCount = []
		for path in allPaths:
			isDependentOn = [None for n in path]
			alteredOnTurn = [0 for n in ref]
			for i in xrange(len(path)):
				#print i,path[i],alteredOnTurn
				if path[i][0] == 'del':
					isDependentOn[i] = alteredOnTurn[path[i][1][0]:path[i][1][1]]
					del alteredOnTurn[path[i][1][0]:path[i][1][1]]
				elif path[i][0] == 'inv':
					isDependentOn[i] = alteredOnTurn[path[i][1][0]:path[i][1][1]]
					for j in xrange(path[i][1][0],path[i][1][1]):
						alteredOnTurn[j] = i+1
				elif path[i][0] == 'ins' or 'dup' in path[i][0]:
					isDependentOn[i] = [min(alteredOnTurn[path[i][1][0]-1:path[i][1][0]+1])]
					alteredOnTurn = alteredOnTurn[:path[i][1][0]] + [i+1 for n in path[i][1][1]] + alteredOnTurn[path[i][1][0]:]
			all_refAlterCount.append([item for sublist in isDependentOn for item in sublist].count(0))
			all_dependencies.append([list(set(n)) for n in isDependentOn])

		# choose to report the event which minimizes (in the following order of priority)...
		#
		#    1) number of total alterations
		#    2) number of insertion operations
		#    3) number of inserted elements
		#    4) number of reference regions altered
		#
		#       if there's a tie, just choose the first option, for consistency
		#
		sorted_paths = []
		for i in xrange(len(allPaths)):
			insCount = sum([len(n[1][1]) for n in allPaths[i] if n[0] == 'ins']+[len(n[1][1]) for n in allPaths[i] if n[0] == 'dup'])
			sort_tuple = (len(allPaths[i]), [n[0] for n in allPaths[i]].count('ins') + [n[0] for n in allPaths[i]].count('dup'), insCount, all_refAlterCount[i])
			sorted_paths.append([sort_tuple,i])
		sorted_paths = sorted(sorted_paths)

		ind_to_report = sorted_paths[0][1]

		#for n in sorted_paths:
		#	print n[0], allPaths[n[1]]

		topPath = allPaths[ind_to_report]
		topDeps = all_dependencies[ind_to_report]
		#print 'topPath:',topPath
		#print 'topDeps:',topDeps

		#
		# output list of independent (can be nested) alterations that we decided describe the event,
		# as well as a text description that can be used for plotting.
		#
		clusters = [[i+1] for i in xrange(len(topPath)) if topDeps[i] == [0]]
		toAdd    = [i for i in xrange(len(topPath)) if topDeps[i] != [0]]
		#print 'clusters:',clusters
		#print 'toAdd:   ',toAdd
		prevLen  = len(toAdd)
		while toAdd:
			for i in xrange(len(toAdd)):
				found = False
				for j in xrange(len(clusters)):
					if all([clusters[j].count(n) for n in topDeps[toAdd[i]]]):
						clusters[j].append(toAdd.pop(i)+1)
						found = True
						break
				if found:
					break
			if len(toAdd) == prevLen:
				print 'Error: There was a problem in clustering:'
				print '---',topPath
				print '---',topDeps
				exit(1)
			prevLen = len(toAdd)
		for i in xrange(len(clusters)):
			clusters[i] = [n-1 for n in clusters[i]]
		print 'clusters:',clusters

		# create description string for human interpretation
		string_description = ''
		for i in xrange(len(clusters)):
			if i > 0:
				string_description += ' + '
			string_description += '['
			for j in xrange(len(clusters[i])):
				if j > 0:
					string_description += ' + '
				string_description += topPath[clusters[i][j]][0]
			string_description += ']'
		print string_description

		# return all the data structures needed to print and further process
		return (topPath,clusters,string_description)

	else:
		#print 'Error: Invalid event comparison'
		#exit(1)
		return ([],[],'[unknown]')


if __name__ == '__main__':
	#tt = time.time()
	#interpret_rearrangement(['anchor1','A','B','C','anchor2'],['anchor1','B*','anchor2'])
	#interpret_rearrangement(['anchor1','A','B','C','anchor2'],['anchor1','A','A','B','C','A','B*','C','anchor2'])
	#interpret_rearrangement(['anchor1','A','B','C','anchor2'],['anchor1','A','B','C','A','C','anchor2'])
	#interpret_rearrangement(['anchor1','A','B','C','anchor2'],['anchor1','Z','anchor2'],novel=['Z'])
	#print time.time()-tt,'(sec)'

	# harder cases
	#
	#interpret_rearrangement(['anchor1', 'A', 'B', 'C', 'anchor2'],['anchor1', 'C*', 'B', 'C', 'anchor2'],novel=[])
	#interpret_rearrangement(['anchor1','A','B','C','D','E','F','anchor2'], ['anchor1','A','B','C','D','A','B','G','F','anchor2'], novel=['G'])
	interpret_rearrangement(['anchor1', 'A', 'B', 'C', 'anchor2'], ['anchor1', 'A', 'B', 'C', 'B*', 'D', 'A', 'B', 'C', 'anchor2'], novel=['D'])



