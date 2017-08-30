import os

def openReadClose(fn):
	f = open(fn,'r')
	fr = f.read()
	f.close()
	return fr

def openWriteClose(fn,s):
	f = open(fn,'w')
	f.write(s)
	f.close()

def rm(f):
	if os.path.isfile(f):
		os.system('rm '+f)

def makedir(d):
	if not os.path.isdir(d):
		os.system('mkdir '+d)

def exe(c):
	#print c
	os.system(c)