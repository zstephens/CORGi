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

# get consensus sequence via MSA
def msa_consensus(sequences):
	if len(sequences) == 1:
		return sequences[0]
	msa_result = msa_optimal(sequences)
	print msa_result
	exit(1)
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
	#seq = seq[2:5]

	msa = msa_optimal(seq)

	for n in msa:
		print n

