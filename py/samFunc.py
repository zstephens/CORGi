import os
import re

from sysFunc import openWriteClose, exe, rm

def indexRef(refPath):
	fn = None
	if os.path.isfile(refPath+'i'):
		#print 'found index '+refPath+'i'
		fn = refPath+'i'
	if os.path.isfile(refPath+'.fai'):
		#print 'found index '+refPath+'.fai'
		fn = refPath+'.fai'

	ref_inds = []
	if fn != None:
		fai = open(fn,'r')
		for line in fai:
			splt = line[:-1].split('\t')
			seqLen = int(splt[1])
			offset = int(splt[2])
			lineLn = int(splt[3])
			nLines = seqLen/lineLn
			if seqLen%lineLn != 0:
				nLines += 1
			ref_inds.append((splt[0],offset,offset+seqLen+nLines,seqLen))
		fai.close()
		return ref_inds

	#sys.stdout.write('index not found, creating one... ')
	#sys.stdout.flush()
	refFile = open(refPath,'r')
	prevR   = None
	prevP   = None
	seqLen  = 0
	while 1:
		data = refFile.readline()
		if not data:
			ref_inds.append( (prevR, prevP, refFile.tell()-len(data), seqLen) )
			break
		if data[0] == '>':
			if prevP != None:
				ref_inds.append( (prevR, prevP, refFile.tell()-len(data), seqLen) )
			seqLen = 0
			prevP  = refFile.tell()
			prevR  = data[1:-1]
		else:
			seqLen += len(data)-1
	refFile.close()
	#print ''
	return ref_inds

def readRef(refPath,ref_inds_i):
	refFile = open(refPath,'r')
	refFile.seek(ref_inds_i[1])
	myDat = ''.join(refFile.read(ref_inds_i[2]-ref_inds_i[1]).split('\n'))
	myDat = bytearray(myDat.upper())
	return myDat

def writeRefSection(refFile, chrom, samPos_s, samPos_e, fn):
	refInds  = indexRef(refFile)
	refIndex = -1
	for i in xrange(len(refInds)):
		if refInds[i][0] == chrom:
			refIndex = i
			break
	if refIndex == -1:
		print '\n\nCOULD NOT FIND SPECIFIED CHROM IN REF FASTA.\n\n'
		exit(1)
	refChunk = readRef(refFile,refInds[refIndex])[samPos_s-1:samPos_e-1]
	f = open(fn,'w')
	f.write('>refSeq_'+chrom+'_'+str(samPos_s)+'_'+str(samPos_e)+'\n'+refChunk+'\n')
	f.close()

def filterSamByClipContent(sam_in,sam_out,minClip,allow_hard=False):
	fi = open(sam_in,'r')
	fo = open(sam_out,'w')
	CLIP_CHAR = 'S'
	if allow_hard:
		CLIP_CHAR += 'H'
	for line in fi:
		s = line.split('\t')[5]
		letters = re.split(r"\d+",s)[1:]
		numbers = [int(n) for n in re.findall(r"\d+",s)]
		if len(letters) == 0 or len(numbers) == 0:
			print 'invalid line?'
			print line
			continue
		#print letters, numbers
		if (letters[0] in CLIP_CHAR and numbers[0] >= minClip) or (letters[-1] in CLIP_CHAR and numbers[-1] >= minClip):
			fo.write(line)
	fo.close()
	fi.close()

def filterSamByRefOverlap(sam_in,sam_out,pos_s,pos_e):
	fi = open(sam_in,'r')
	fo = open(sam_out,'w')
	CLIP_CHAR = 'SH'
	REF_CHAR  = 'MX=D'
	MIN_MAPQ  = 3
	for line in fi:
		s = line.split('\t')[5]
		if int(line.split('\t')[4]) >= MIN_MAPQ:
			letters = re.split(r"\d+",s)[1:]
			numbers = [int(n) for n in re.findall(r"\d+",s)]
			overlap = False
			bOffset = 0
			map_ind = int(line.split('\t')[3])
			for i in xrange(len(letters)):
				if (i==0 and letters[i] in CLIP_CHAR) or (i==len(letters)-1 and letters[i] in CLIP_CHAR):
					pass
				elif letters[i] in REF_CHAR:
					for j in xrange(numbers[i]):
						if map_ind+bOffset >= pos_s and map_ind+bOffset < pos_e:
							overlap = True
							break
						bOffset += 1
			if overlap:
				fo.write(line)
	fo.close()
	fi.close()

def filterSamByReadLen(sam_in,sam_out,rMax=1e12,rMin=-1):
	fi = open(sam_in,'r')
	fo = open(sam_out,'w')
	for line in fi:
		rLen = len(line.split('\t')[9])
		if rLen >= rMin and rLen <= rMax:
			fo.write(line)
	fo.close()
	fi.close()

def sam2fa(f_in,f_out):
	f1 = open(f_in,'r')
	f2 = open(f_out,'w')
	for line in f1:
		if line[0] != '@':
			splt = line.strip().split('\t')
			f2.write('>'+splt[0]+'_'+str(len(splt[9]))+'\n'+splt[9]+'\n')
	f2.close()
	f1.close()

def readBlastLine(line):
	splt = line.strip().split('\t')
	return (splt[0],int(splt[3]),int(splt[6]),int(splt[7]),int(splt[8]),int(splt[9]),float(splt[10]),float(splt[11]))

IGV_EXE      = 'java -Xmx1500m -jar /Volumes/Epoch9/tools/IGV_2.3.57/igv.jar'
TEMP_IGV     = 'temp_igv_batch.txt'
TEMP_IGLOG   = 'temp_igvLog.txt'
def view_IGV(in_bam,chrom,ind_s,ind_e,screenShotDir=None):
	igvStr    = 'load '+in_bam+'\n'
	if screenShotDir != None:
		igvStr   += 'snapshotDirectory '+screenShotDir+'\n'
	igvStr   += 'genome hg19\n'
	igvStr   += 'goto '+chrom+':'+str(ind_s)+'-'+str(ind_e)+'\n'
	#igvStr   += 'viewaspairs\n'
	#igvStr   += 'sort insert\n'
	igvStr   += 'collapse\n'
	if screenShotDir != None:
		igvStr   += 'snapshot\n'
		igvStr   += 'exit\n'
	openWriteClose(TEMP_IGV,igvStr)
	igvCMD = IGV_EXE + ' -b ' + TEMP_IGV + ' > ' + TEMP_IGLOG + ' 2>&1'
	exe(igvCMD)
	rm(TEMP_IGV)
	rm(TEMP_IGLOG)

