#!/usr/bin/env python

from numpy import *

import sonIO
import neuralynx_io
import klusters_io
from spikeFeatures import *
import spikePreprocess.computeSpikeFeatures as computeSpikeFeatures
import pymatfile


FEATURES = ('AvWaveform: ', 'StdWaveform: ')
N_PROBES = 4

def getAvWaveform(files, clufile, clusterIdx, nProbes, features):
	"""run (as a simple script) all the steps necessary to prepare for running KK"""

	featuresToCompute = []
	for c in range(1, nProbes+1):
		for f in features:
			featuresToCompute.append(f + str(c))
	
	# create data file handler objects 
	spikeFiles = []
	if files[0][-3:] == 'smr':
		fileType = 'son'
	elif files[0][-5:-2] == 'spk':
		fileType = 'klusters'
	else:
		fileType = 'nlynx'
	
	if fileType == 'son':
		for fname in files:
			spikeFiles.append(sonIO.SONfile(fname))
	elif fileType == 'nlynx':
		for fname in files:
			spikeFiles.append(neuralynx_io.TTFile(fname))
	elif fileType  == 'klusters':
		print('klusters_io')
		for fname in files:
			spikeFiles.append(klusters_io.TTFile(fname))
	if fileType == 'son':
		adcMarkChans = spikeFiles[0].channelsOfKind('ADCMarkerChannel')
	elif fileType == 'nlynx':
		adcMarkChans = (0,)
		
	name = files[0][:-4]
	name = name + "_c_" + clusterIdx + 'wv.mat'
	print name
	outfh = pymatfile.pymatfile(name, 'w')
	if fileType == 'son':
		ttFiles = []
		for sp in spikeFiles:
			ttFiles.append(sonIO.TTFile(sp))
	else:
		ttFiles = spikeFiles
	print(ttFiles)
	computer = computeSpikeFeatures.FeatureComputer(ttFiles, spkname=None, params=None, interleave=120)
	computer.setMaskFromCluFile(clufile, clusterIdx)
	computer.computeFeatures(featuresToCompute)
	features = computer.computedFeatures
	a = array([])
	for f in featuresToCompute:
		fg = reshape(features[f], (1, len(features[f])))
		if len(a):
			a = concatenate([a, fg], 0)
		else:
			a = fg
	outfh.writeItem(a, 'waveformstd')
	outfh.close()
		
if __name__ == '__main__':
	import sys, os

	if len(sys.argv) != 4:
		print "error: call as get_av_waveform TTfile cluFile clusterIdx"
		sys.exit(-1)

		
	getAvWaveform(files = sys.argv[1:-2], clufile = sys.argv[-2], clusterIdx = sys.argv[-1], nProbes=N_PROBES, 
		features=FEATURES)
	
