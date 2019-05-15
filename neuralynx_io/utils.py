"""module for input/output utilities for Neuralynx files """
from numpy import *
import _TTfile

def subsetTTfile(filenameIn, filenameOut, time1, time2):
	"""\brief writes a TT file with a subset of spikes in the input file with times from time 1 to time2"""
	inFh = _TTfile.TTFile(filenameIn)
	outFh = _TTfile.TTFile(filenameOut, 'w')
	for (t, wv, istart, iend, wvRaw) in inFh.genTTData():
		ix = where (t > time1) & (t < time2)
		outFh.writeSpikes(t[ix], wv[ix, :])
	outFh.close()