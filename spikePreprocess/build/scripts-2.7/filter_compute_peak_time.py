#!/home/tara/miniconda3/envs/spikesorting/bin/python

from numpy import *
import sonIO
import neuralynx_io
from spikeFeatures import *
import spikePreprocess.computeSpikeFeatures as computeSpikeFeatures



FEATURES_PER_TRODE = ( 'Peak: ',)
FEATURES_TIM = ('Time in File (1/10000 s)', 'Time (1/10000 s)', 'File ID' )
#FEATURES_MASK = {'Peak Index: ': (-10, 500)}
FEATURES_MASK = {}
N_PROBES = 4

def runKKFeaturesFromFiles(name, files, nProbes, featuresPerProbe, featuresTim, featuresMask ):
	"""run (as a simple script) all the steps necessary to prepare for running KK"""
	fetFeatures = []
	maskFeatures = {}
	fM = featuresMask.keys()
	for c in range(1, nProbes+1):
		for f in featuresPerProbe:
			fetFeatures.append(f + str(c))
		for f in fM:
			maskFeatures[f + str(c)] = featuresMask[f]
	fetFeatures.append('Time in File (1/10000 s)')
	fetFeatures = tuple(fetFeatures)
	print fetFeatures	
	featuresToCompute = []
	featuresToCompute.extend(fetFeatures)
	featuresToCompute.extend(featuresTim)
	
	

	#TODO add mean waveforms
	
	# create data file handler objects 
	spikeFiles = []
	if files[-3:] == 'smr':
		fileType = 'son'
	else:
		fileType = 'nlynx'
	
	if fileType == 'son':
		for fname in files:
			spikeFiles.append(sonIO.SONfile(fname))
	elif fileType == 'nlynx':
		for fname in files:
			spikeFiles.append(neuralynx_io.TTFile(fname))
	
	if fileType == 'son':
		adcMarkChans = spikeFiles[0].channelsOfKind('ADCMarkerChannel')
	elif fileType == 'nlynx':
		adcMarkChans = (0,)
		
	for chan in adcMarkChans:
		fname = name + ".peak.1" #principal components file name needed
		# by klusters
		spkname = name + ".spk.1" # spike waveform file, needed by
		# klusters
		parname = name + ".xml"
		timname = name + ".tim.1"
		if fileType == 'son':
			ttFiles = []
			for sp in spikeFiles:
				ttFiles.append(sonIO.TTFile(sp))
		else:
			ttFiles = spikeFiles
		computer = computeSpikeFeatures.FeatureComputer(ttFiles, spkname, params=None, interleave=120)
		computer.setMaskFromFeatures(maskFeatures)
		computer.computeFeatures(featuresToCompute)
		features = computer.computedFeatures
		computeSpikeFeatures.writeFetFile(fname, features, fetFeatures)
		computeSpikeFeatures.writeTimFile(timname, features)
		computeSpikeFeatures.makeXMLParFile(spikeFiles, parname, len(featuresPerProbe))




if __name__ == '__main__':
	import sys, os

	if len(sys.argv) == 1:
		os.chdir('testData')
		runKKFeaturesFromFiles(name = 'tt1', files = ('tt1_test.ntt', ), nProbes=N_PROBES, 
			featuresPerProbe=FEATURES_PER_TRODE, featuresTim=FEATURES_TIM, featuresMask=FEATURES_MASK)
		sys.exit(0)
		
	runKKFeaturesFromFiles(name = sys.argv[1], files = sys.argv[2:], nProbes=N_PROBES, 
		featuresPerProbe=FEATURES_PER_TRODE, featuresTim=FEATURES_TIM, featuresMask=FEATURES_MASK)
	
