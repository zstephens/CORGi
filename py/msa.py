from operator import mul
import numpy as np

MATCH    = 10
MISMATCH = -1000
GAP      = -10
GAP_CHR  = '-'

def S(xList):
	myScore = 0
	for i in xrange(len(xList)):
		if xList[i] == GAP_CHR:
			myScore += GAP
		else:
			for j in xrange(i+1,len(xList)):
				if xList[j] != GAP:
					if xList[i] == xList[j]:
						myScore += MATCH
					else:
						myScore += MISMATCH
	return myScore

#
#	OPTIMAL MSA ALIGNMENT USING A REALLY CUMBERSOME DYNAMIC PROGRAMMING APPROACH
#
def msa_optimal(sequences):

	# degenerate condition
	if len(sequences) == 1:
		return sequences

	N = [len(n) for n in sequences]
	D = len(sequences)

	score = np.zeros([n+1 for n in N])
	trace = np.zeros([n+1 for n in N],dtype='<i4')

	#
	#	CREATE DELTA VECTORS
	#
	delta_vec = [tuple([-1 for n in N])]
	ranges     = [(0,2) for n in N]
	ops        = reduce(mul,(p[1]-p[0] for p in ranges))-1
	result     = [n[0] for n in ranges]
	pos        = len(ranges)-1
	increments = 0
	while increments < ops:
		if result[pos] == ranges[pos][1]-1:
			result[pos] = ranges[pos][0]
			pos        -= 1
		else:
			result[pos] += 1
			increments  += 1
			pos          = len(ranges)-1 #increment the innermost loop
			if sum(result) < len(result):
				delta_vec.append(tuple([delta_vec[0][n]+result[n] for n in xrange(D)]))

	#
	#	BIG MATRIX TRAVERSAL
	#
	resultList = [[0 for n in xrange(D)]]
	seenBefore = {}
	dirList    = [n for n in delta_vec if sum(n) == -1]
	while resultList:
		areDone = False
		basePos = resultList.pop()
		for dv in dirList:
			result = [basePos[n]-dv[n] for n in xrange(D)]
			if tuple(result) in seenBefore:
				continue
			if len([n for n in result if n < 0]) or len([n for n in xrange(D) if result[n] > N[n]]):
				continue
			else:
				resultList.insert(0,[n for n in result])
				seenBefore[tuple(resultList[0])] = True
		if basePos != [0 for n in xrange(D)]:
			#print '--', basePos
			result = [n for n in basePos]

			scoreList = []
			for dvi in xrange(len(delta_vec)):
				dv = delta_vec[dvi]
				newPos = [result[n]+dv[n] for n in xrange(D)]
				if len([n for n in newPos if n < 0]):
					continue
				#print '-',dv
				xList = []
				for i in xrange(D):
					if dv[i] == -1:
						xList.append(sequences[i][newPos[i]])
					else:
						xList.append(GAP_CHR)
				scoreList.append((score[tuple(newPos)] + S(xList),dvi))
			scoreList = sorted(scoreList,reverse=True)
			score[tuple(result)] = scoreList[0][0]
			trace[tuple(result)] = scoreList[0][1]

	#
	#	TRACEBACK
	#
	pos = [n for n in N]
	msa_output = ['' for n in xrange(D)]
	while True:
		dv = delta_vec[trace[tuple(pos)]]
		#print pos, dv
		for i in xrange(D):
			if dv[i] == 0:
				msa_output[i] += GAP_CHR
			elif dv[i] == -1:
				msa_output[i] += sequences[i][pos[i]-1]
		pos = [pos[n]+dv[n] for n in xrange(D)]
		if pos == [0]*D:
			break
	for i in xrange(D):
		msa_output[i] = msa_output[i][::-1]

	return msa_output

# nw alignment
def needleman_wunsch(mismatch, gap, sequence1, sequence2):
	n = len(sequence1)
	m = len(sequence2)
	subproblems = [[0 for x in range(m+1)] for x in range(n+1)]
	for i in range(n+1):
		subproblems[i][0] = i * gap
	for j in range(m+1):
		subproblems[0][j] = j * gap

	for i in range(1, n+1):
		for j in range(1, m+1):
			case1 = subproblems[i-1][j-1]
			if sequence1[i-1] != sequence2[j-1]:
				case1 += mismatch
			case2 = subproblems[i-1][j] + gap
			case3 = subproblems[i][j-1] + gap
			subproblems[i][j] = min([case1, case2, case3])
	penalty = subproblems[n][m]

	# Backtrace
	alignment1 = ""
	alignment2 = ""
	i = n
	j = m
	while i > 0 or j > 0:
		pos = subproblems[i][j]
		case1_match = subproblems[i-1][j-1]
		case1_mismatch = case1_match + mismatch
		case2 = subproblems[i-1][j] + gap
		case3 = subproblems[i][j-1] + gap
		if i > 0 and pos == case1_match:
			alignment1 = sequence1[i-1] + alignment1
			alignment2 = sequence2[j-1] + alignment2
			i -= 1
			j -= 1
		elif i > 0 and pos == case1_mismatch:
			alignment1 = sequence1[i-1] + alignment1
			alignment2 = sequence2[j-1] + alignment2
			i -= 1
			j -= 1
		elif i > 0 and pos == case2:
			alignment1 = sequence1[i-1] + alignment1
			alignment2 = GAP_CHR + alignment2
			i -= 1
		elif j > 0 and pos == case3:
			alignment1 = GAP_CHR + alignment1
			alignment2 = sequence2[j-1] + alignment2
			j -= 1
	return (penalty, alignment1, alignment2)

##### a much sloppier version of MSA
####def msa_sloppy_consensus(sequences):
####	seq_sort = [m[1] for m in sorted([(len(n),n) for n in sequences],reverse=True)]
####	print seq_sort
####	cons = msa_consensus(seq_sort[0:2])
####	print cons
####	for i in xrange(2,len(seq_sort)):
####		inS  = ''.join([n[0] for n in cons])
####		cons = msa_consensus([inS,seq_sort[i]])
####		print cons

# get consensus sequence via MSA
def msa_consensus(sequences):
	if len(sequences) == 1:
		return sequences[0]
	msa_result = msa_optimal(sequences)
	#print msa_result
	consensus  = []
	for j in xrange(len(msa_result[0])):
		col = [msa_result[i][j] for i in xrange(len(msa_result))]
		cDict = {}
		for n in col:
			if n not in cDict:
				cDict[n] = 0
			cDict[n] += 1
		topVal = max([[cDict[k],k] for k in cDict.keys()])[::-1]
		consensus.append(tuple(topVal))
	return consensus


if __name__ == '__main__':
	seq = ['TCAGGATGAAC',
		   'ATCACGATGAACC',
		   'ATCAGGAATGAATCC',
		   'TCACGATTGAATCGC',
		   'TCAGGAATGAATCGCM']
	seq = seq[0:1]+seq[2:4]

	seq = ['BCDEFCDDEFG','CDDEFG','CDDE','HABCDEFCDDE','HABCDEFCDDEFG']

	msa = msa_optimal(seq)
	for n in msa:
		print n

	print 'optimal:'
	print msa_consensus(seq)
	#print 'sloppy:'
	#print msa_sloppy_consensus(seq)

